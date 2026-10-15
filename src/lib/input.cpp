/*
input.cpp: trumut input readers.
*/
#include "input.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <zlib.h>

#include <condition_variable>
#include <deque>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "third_party/kseq.h"

struct BgzfBlock {
    uint64_t seq;
    std::vector<uint8_t> cdata;
    uint32_t crc, isize;
    std::vector<uint8_t> out;
};

struct BgzfPool {
    int fd;
    int n_threads;
    std::vector<uint8_t> carry;

    std::mutex mu;
    std::condition_variable cv_work, cv_done, cv_space;
    std::deque<BgzfBlock *> work;
    std::map<uint64_t, BgzfBlock *> done;
    size_t in_flight = 0;
    size_t max_in_flight;
    uint64_t next_seq_read = 0;     // sequence number the reader gives out next
    uint64_t next_seq_serve = 0;    // the block the caller wants next
    bool reader_done = false;
    bool stop = false;
    bool failed = false;
    std::string err;
    uint32_t last_isize = 1;

    std::thread reader;
    std::vector<std::thread> workers;

    BgzfBlock *cur = nullptr;       // block the caller is consuming
    size_t cur_pos = 0;
};

static ssize_t pool_read_exact(BgzfPool *p, uint8_t *buf, size_t n, std::string &err)
{
    size_t got = 0;
    if (!p->carry.empty()) {
        got = p->carry.size() < n ? p->carry.size() : n;
        memcpy(buf, p->carry.data(), got);
        p->carry.erase(p->carry.begin(), p->carry.begin() + (long)got);
    }
    while (got < n) {
        ssize_t k = read(p->fd, buf + got, n - got);
        if (k < 0) {
            if (errno == EINTR) continue;
            err = std::string("read failed: ") + strerror(errno);
            return -1;
        }
        if (k == 0) break;
        got += (size_t)k;
    }
    if (got == 0) return 0;
    if (got < n) { err = "input ends inside a BGZF block (file is truncated)"; return -1; }
    return (ssize_t)n;
}

static inline uint16_t le16(const uint8_t *p) { return (uint16_t)(p[0] | (p[1] << 8)); }
static inline uint32_t le32(const uint8_t *p) {
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

/* parse one BGZF block off the descriptor into *out. */
static int pool_read_block(BgzfPool *p, BgzfBlock *b, std::string &err)
{
    uint8_t hdr[12];
    ssize_t got = pool_read_exact(p, hdr, 12, err);
    if (got <= 0) return (int)got;
    if (hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 8 || !(hdr[3] & 4)) {
        err = "not a BGZF block where one was expected";
        return -1;
    }
    uint16_t xlen = le16(hdr + 10);
    std::vector<uint8_t> extra(xlen);
    if (xlen && pool_read_exact(p, extra.data(), xlen, err) <= 0) {
        if (err.empty()) err = "truncated BGZF header";
        return -1;
    }
    int bsize = -1;
    for (size_t i = 0; i + 4 <= extra.size();) {
        uint8_t si1 = extra[i], si2 = extra[i + 1];
        uint16_t slen = le16(&extra[i + 2]);
        if (si1 == 'B' && si2 == 'C' && slen == 2 && i + 6 <= extra.size()) {
            bsize = le16(&extra[i + 4]);
            break;
        }
        i += 4 + slen;
    }
    if (bsize < 0) { err = "BGZF block has no BC size field"; return -1; }
    // total block size = bsize + 1 = 12 + xlen + cdata + 8
    long cdata_len = (long)bsize + 1 - 12 - xlen - 8;
    if (cdata_len < 0) { err = "bad BGZF block size"; return -1; }
    b->cdata.resize((size_t)cdata_len);
    if (cdata_len && pool_read_exact(p, b->cdata.data(), (size_t)cdata_len, err) <= 0) {
        if (err.empty()) err = "truncated BGZF payload";
        return -1;
    }
    uint8_t trailer[8];
    if (pool_read_exact(p, trailer, 8, err) <= 0) {
        if (err.empty()) err = "truncated BGZF trailer";
        return -1;
    }
    b->crc = le32(trailer);
    b->isize = le32(trailer + 4);
    if (b->isize > 65536) { err = "BGZF block claims more than 64 KiB"; return -1; }
    return 1;
}

static void pool_fail(BgzfPool *p, const std::string &msg)
{
    std::lock_guard<std::mutex> lk(p->mu);
    if (!p->failed) { p->failed = true; p->err = msg; }
    p->reader_done = true;
    p->cv_done.notify_all();
    p->cv_work.notify_all();
    p->cv_space.notify_all();
}

static void pool_reader_main(BgzfPool *p)
{
    for (;;) {
        BgzfBlock *b = new BgzfBlock();
        std::string err;
        int r = pool_read_block(p, b, err);
        if (r <= 0) {
            delete b;
            if (r < 0) { pool_fail(p, err); return; }
            std::lock_guard<std::mutex> lk(p->mu);
            p->reader_done = true;
            p->cv_work.notify_all();
            p->cv_done.notify_all();
            return;
        }
        std::unique_lock<std::mutex> lk(p->mu);
        p->cv_space.wait(lk, [&] { return p->stop || p->failed || p->in_flight < p->max_in_flight; });
        if (p->stop || p->failed) { delete b; return; }
        b->seq = p->next_seq_read++;
        p->last_isize = b->isize;
        p->in_flight++;
        p->work.push_back(b);
        p->cv_work.notify_one();
    }
}

static void pool_worker_main(BgzfPool *p)
{
    z_stream zs;
    memset(&zs, 0, sizeof zs);
    if (inflateInit2(&zs, -15) != Z_OK) { pool_fail(p, "inflateInit2 failed"); return; }
    for (;;) {
        BgzfBlock *b;
        {
            std::unique_lock<std::mutex> lk(p->mu);
            p->cv_work.wait(lk, [&] { return p->stop || p->failed || !p->work.empty() || p->reader_done; });
            if (p->stop || p->failed) break;
            if (p->work.empty()) { if (p->reader_done) break; continue; }
            b = p->work.front();
            p->work.pop_front();
        }
        b->out.resize(b->isize);
        if (b->isize > 0) {
            inflateReset(&zs);
            zs.next_in = b->cdata.data();
            zs.avail_in = (uInt)b->cdata.size();
            zs.next_out = b->out.data();
            zs.avail_out = b->isize;
            int ret = inflate(&zs, Z_FINISH);
            if (ret != Z_STREAM_END || zs.avail_out != 0) {
                delete b;
                pool_fail(p, "BGZF block failed to inflate (file is corrupt)");
                break;
            }
            if (crc32(0L, b->out.data(), b->isize) != b->crc) {
                delete b;
                pool_fail(p, "BGZF block CRC mismatch (file is corrupt)");
                break;
            }
        }
        std::lock_guard<std::mutex> lk(p->mu);
        p->done[b->seq] = b;
        p->cv_done.notify_all();
    }
    inflateEnd(&zs);
}

static BgzfBlock *pool_next_block(BgzfPool *p)
{
    std::unique_lock<std::mutex> lk(p->mu);
    for (;;) {
        auto it = p->done.find(p->next_seq_serve);
        if (it != p->done.end()) {
            BgzfBlock *b = it->second;
            p->done.erase(it);
            p->next_seq_serve++;
            p->in_flight--;
            p->cv_space.notify_one();
            return b;
        }
        if (p->failed) return nullptr;
        if (p->reader_done && p->next_seq_serve >= p->next_seq_read) return nullptr;
        p->cv_done.wait(lk);
    }
}

static void pool_start(BgzfPool *p)
{
    p->max_in_flight = (size_t)p->n_threads * 8 + 8;
    for (int i = 0; i < p->n_threads; ++i) p->workers.emplace_back(pool_worker_main, p);
    p->reader = std::thread(pool_reader_main, p);
}

static void pool_stop(BgzfPool *p)
{
    {
        std::lock_guard<std::mutex> lk(p->mu);
        p->stop = true;
        p->cv_work.notify_all();
        p->cv_done.notify_all();
        p->cv_space.notify_all();
    }
    if (p->reader.joinable()) p->reader.join();
    for (auto &t : p->workers) if (t.joinable()) t.join();
    for (auto *b : p->work) delete b;
    for (auto &kv : p->done) delete kv.second;
    delete p->cur;
}

struct seqsrc_s {
    int fd = -1;
    int close_fd = 0;
    pid_t child = 0;
    std::string fetcher;
    int kind = SEQSRC_RAW;
    std::string path;

    std::string err;
    bool eof = false;

    std::vector<uint8_t> in;
    size_t in_pos = 0;
    bool raw_eof = false;

    // gzip stream mode
    z_stream zs;
    bool zs_init = false;
    bool mid_member = false;
    uint64_t cur_member = 0, last_member = 1;

    // bgzf parallel mode
    BgzfPool *pool = nullptr;

    std::vector<uint8_t> ahead;
    size_t ahead_pos = 0;
};

static void src_seterr(seqsrc_t *s, const std::string &m) { if (s->err.empty()) s->err = m; }

static int has_prefix(const char *p, const char *pre) { return strncmp(p, pre, strlen(pre)) == 0; }

int seqsrc_url_kind(const char *path)
{
    if (!path) return 0;
    if (has_prefix(path, "http://") || has_prefix(path, "https://") || has_prefix(path, "ftp://")) return 1;
    if (has_prefix(path, "s3://")) return 2;
    if (has_prefix(path, "gs://")) return 3;
    return 0;
}

static int exec_on_path(const char *prog)
{
    const char *path = getenv("PATH");
    if (!path) path = "/usr/local/bin:/usr/bin:/bin";
    std::string dirs(path);
    size_t a = 0;
    while (a <= dirs.size()) {
        size_t b = dirs.find(':', a);
        if (b == std::string::npos) b = dirs.size();
        std::string cand = dirs.substr(a, b - a);
        if (!cand.empty()) {
            cand += '/'; cand += prog;
            struct stat st;
            if (stat(cand.c_str(), &st) == 0 && S_ISREG(st.st_mode) && access(cand.c_str(), X_OK) == 0) return 1;
        }
        a = b + 1;
    }
    return 0;
}

static int spawn_fetcher(const char *url, int kind, pid_t *pid_out, std::string *prog_out,
                         char *err, size_t err_n)
{
    const char *argv[8]; int n = 0;
    if (kind == 1) {
        argv[n++] = "curl"; argv[n++] = "-sSfL"; argv[n++] = "--retry"; argv[n++] = "3"; argv[n++] = url;
    } else if (kind == 2) {
        argv[n++] = "aws"; argv[n++] = "s3"; argv[n++] = "cp"; argv[n++] = url; argv[n++] = "-";
    } else {
        argv[n++] = "gcloud"; argv[n++] = "storage"; argv[n++] = "cat"; argv[n++] = url;
    }
    argv[n] = 0;
    if (!exec_on_path(argv[0])) {
        if (err) snprintf(err, err_n, "%s is not on PATH, and it is what fetches %s URLs. Install it, "
                          "pass a presigned https URL, or pipe the file in: %s ... | trumut allele - ...",
                          argv[0], kind == 1 ? "http/ftp" : kind == 2 ? "s3://" : "gs://",
                          kind == 2 ? "aws s3 cp URL -" : kind == 3 ? "gcloud storage cat URL" : "curl -sSfL URL");
        return -1;
    }
    int fds[2];
    if (pipe(fds) != 0) { if (err) snprintf(err, err_n, "pipe: %s", strerror(errno)); return -1; }
    pid_t pid = fork();
    if (pid < 0) { if (err) snprintf(err, err_n, "fork: %s", strerror(errno)); close(fds[0]); close(fds[1]); return -1; }
    if (pid == 0) {
        close(fds[0]);
        dup2(fds[1], STDOUT_FILENO);
        close(fds[1]);
        execvp(argv[0], (char *const *)argv);
        _exit(127);
    }
    close(fds[1]);
    *pid_out = pid;
    *prog_out = argv[0];
    return fds[0];
}

static int reap_fetcher(pid_t pid)
{
    int st = 0;
    while (waitpid(pid, &st, 0) < 0) if (errno != EINTR) return -1;
    if (WIFEXITED(st)) return WEXITSTATUS(st);
    if (WIFSIGNALED(st)) return 128 + WTERMSIG(st);
    return -1;
}

static int src_raw_refill(seqsrc_t *s)
{
    if (s->in_pos) { s->in.erase(s->in.begin(), s->in.begin() + (long)s->in_pos); s->in_pos = 0; }
    if (s->raw_eof) return 0;
    size_t old = s->in.size();
    s->in.resize(old + (1u << 20));
    for (;;) {
        ssize_t k = read(s->fd, s->in.data() + old, s->in.size() - old);
        if (k < 0) {
            if (errno == EINTR) continue;
            s->in.resize(old);
            src_seterr(s, std::string("read failed: ") + strerror(errno));
            return -1;
        }
        s->in.resize(old + (size_t)k);
        if (k == 0) { s->raw_eof = 1; return 0; }
        return 1;
    }
}

static ssize_t src_produce(seqsrc_t *s, uint8_t *buf, size_t n)
{
    if (!s->err.empty()) return -1;
    if (s->kind == SEQSRC_BGZF && s->pool) {
        BgzfPool *p = s->pool;
        size_t got = 0;
        while (got < n) {
            if (!p->cur || p->cur_pos >= p->cur->out.size()) {
                delete p->cur;
                p->cur = pool_next_block(p);
                p->cur_pos = 0;
                if (!p->cur) {
                    if (p->failed) { src_seterr(s, p->err); return -1; }
                    break;
                }
                continue;
            }
            size_t k = p->cur->out.size() - p->cur_pos;
            if (k > n - got) k = n - got;
            memcpy(buf + got, p->cur->out.data() + p->cur_pos, k);
            p->cur_pos += k;
            got += k;
        }
        return (ssize_t)got;
    }

    if (s->kind == SEQSRC_RAW) {
        size_t got = 0;
        while (got < n) {
            if (s->in_pos == s->in.size()) {
                int r = src_raw_refill(s);
                if (r < 0) return -1;
                if (r == 0) break;
            }
            size_t k = s->in.size() - s->in_pos;
            if (k > n - got) k = n - got;
            memcpy(buf + got, s->in.data() + s->in_pos, k);
            s->in_pos += k;
            got += k;
        }
        return (ssize_t)got;
    }

    // gzip (and single-threaded BGZF) stream
    size_t got = 0;
    while (got < n) {
        if (s->in_pos == s->in.size()) {
            int r = src_raw_refill(s);
            if (r < 0) return -1;
            if (r == 0) {
                if (s->mid_member) {
                    src_seterr(s, "input ends inside a compressed block (file is truncated)");
                    return -1;
                }
                break;
            }
        }
        s->zs.next_in = s->in.data() + s->in_pos;
        s->zs.avail_in = (uInt)(s->in.size() - s->in_pos);
        s->zs.next_out = buf + got;
        s->zs.avail_out = (uInt)(n - got);
        uInt ai = s->zs.avail_in, ao = s->zs.avail_out;
        int ret = inflate(&s->zs, Z_NO_FLUSH);
        size_t consumed = ai - s->zs.avail_in, produced = ao - s->zs.avail_out;
        s->in_pos += consumed;
        got += produced;
        s->cur_member += produced;
        if (consumed || produced) s->mid_member = true;
        if (ret == Z_STREAM_END) {
            s->last_member = s->cur_member;
            s->cur_member = 0;
            s->mid_member = false;
            if (inflateReset(&s->zs) != Z_OK) { src_seterr(s, "inflateReset failed"); return -1; }
            continue;
        }
        if (ret == Z_BUF_ERROR) continue;
        if (ret != Z_OK) {
            src_seterr(s, std::string("gzip inflate failed: ") + (s->zs.msg ? s->zs.msg : "corrupt data"));
            return -1;
        }
    }
    return (ssize_t)got;
}

seqsrc_t *seqsrc_open(const char *path, int inflate_threads, char *err, size_t err_n)
{
    if (err && err_n) err[0] = 0;
    seqsrc_t *s = new seqsrc_t();
    s->path = (!path || strcmp(path, "-") == 0) ? "stdin" : path;
    int url_kind = path ? seqsrc_url_kind(path) : 0;
    if (!path || strcmp(path, "-") == 0) {
        s->fd = STDIN_FILENO;
    } else if (url_kind) {
        s->fd = spawn_fetcher(path, url_kind, &s->child, &s->fetcher, err, err_n);
        if (s->fd < 0) { delete s; return nullptr; }
        s->close_fd = 1;
    } else {
        s->fd = open(path, O_RDONLY);
        if (s->fd < 0) {
            if (err) snprintf(err, err_n, "cannot open %s: %s", path, strerror(errno));
            delete s;
            return nullptr;
        }
        s->close_fd = 1;
    }

    // sniff
    while (s->in.size() - s->in_pos < 64 && !s->raw_eof) {
        if (src_raw_refill(s) < 0) {
            if (err) snprintf(err, err_n, "%s: %s", s->path.c_str(), s->err.c_str());
            if (s->close_fd) close(s->fd);
            delete s;
            return nullptr;
        }
    }
    const uint8_t *m = s->in.data() + s->in_pos;
    size_t got = s->in.size() - s->in_pos;

    if (got == 0 && s->child) {
        close(s->fd); s->fd = -1;
        int st = reap_fetcher(s->child); s->child = 0;
        if (err) snprintf(err, err_n, "%s: no data from %s (exit status %d; see its message above)",
                          s->path.c_str(), s->fetcher.c_str(), st);
        delete s;
        return nullptr;
    }

    if (got >= 2 && m[0] == 0x1f && m[1] == 0x8b) {
        s->kind = SEQSRC_GZIP;
        if (got >= 12 && (m[3] & 4)) {
            uint16_t xlen = le16(m + 10);
            if (got >= 12u + xlen) {
                for (size_t i = 12; i + 4 <= 12u + xlen;) {
                    uint16_t slen = le16(m + i + 2);
                    if (m[i] == 'B' && m[i + 1] == 'C' && slen == 2) { s->kind = SEQSRC_BGZF; break; }
                    i += 4 + slen;
                }
            }
        }
        if (s->kind == SEQSRC_BGZF && inflate_threads > 1) {
            BgzfPool *p = new BgzfPool();
            p->fd = s->fd;
            p->n_threads = inflate_threads;
            p->carry.assign(s->in.begin() + (long)s->in_pos, s->in.end());
            s->in.clear(); s->in_pos = 0;
            s->pool = p;
            pool_start(p);
        } else {
            memset(&s->zs, 0, sizeof s->zs);
            if (inflateInit2(&s->zs, 15 + 16) != Z_OK) {
                if (err) snprintf(err, err_n, "inflateInit2 failed");
                if (s->close_fd) close(s->fd);
                delete s;
                return nullptr;
            }
            s->zs_init = true;
        }
    } else {
        s->kind = SEQSRC_RAW;
    }
    return s;
}

int seqsrc_kind(const seqsrc_t *s) { return s->kind; }

const char *seqsrc_error(const seqsrc_t *s) { return s->err.c_str(); }

ssize_t seqsrc_read(seqsrc_t *s, void *vbuf, size_t n)
{
    uint8_t *buf = (uint8_t *)vbuf;
    size_t got = 0;
    while (got < n && s->ahead_pos < s->ahead.size()) buf[got++] = s->ahead[s->ahead_pos++];
    if (s->ahead_pos == s->ahead.size()) { s->ahead.clear(); s->ahead_pos = 0; }
    if (got < n && !s->eof) {
        ssize_t k = src_produce(s, buf + got, n - got);
        if (k < 0) return got ? (ssize_t)got : -1;
        if ((size_t)k < n - got) s->eof = true;
        got += (size_t)k;
    }
    return (ssize_t)got;
}

ssize_t seqsrc_peek(seqsrc_t *s, void *vbuf, size_t n)
{
    size_t have = s->ahead.size() - s->ahead_pos;
    if (have < n && !s->eof) {
        std::vector<uint8_t> tmp(n - have);
        ssize_t k = src_produce(s, tmp.data(), tmp.size());
        if (k < 0) return -1;
        if ((size_t)k < tmp.size()) s->eof = true;
        s->ahead.insert(s->ahead.end(), tmp.begin(), tmp.begin() + k);
        have = s->ahead.size() - s->ahead_pos;
    }
    size_t k = have < n ? have : n;
    memcpy(vbuf, s->ahead.data() + s->ahead_pos, k);
    return (ssize_t)k;
}

int seqsrc_close(seqsrc_t *s)
{
    if (!s) return 0;
    int rc = s->err.empty() ? 0 : 1;
    if (s->pool) {
        if (rc == 0 && s->eof) {
            std::lock_guard<std::mutex> lk(s->pool->mu);
            if (s->pool->reader_done && !s->pool->failed && s->pool->last_isize != 0) {
                s->err = "BGZF input ends without its end-of-file block (file is truncated)";
                rc = 1;
            }
        }
        pool_stop(s->pool);
        delete s->pool;
        s->pool = nullptr;
    } else if (s->kind == SEQSRC_BGZF && rc == 0 && s->eof && s->last_member != 0) {
        s->err = "BGZF input ends without its end-of-file block (file is truncated)";
        rc = 1;
    }
    if (s->zs_init) inflateEnd(&s->zs);
    if (s->close_fd && s->fd >= 0) close(s->fd);
    if (s->child) {
        int st = reap_fetcher(s->child); s->child = 0;
        if (rc == 0 && s->eof && st != 0) {
            s->err = s->fetcher + " exited with status " + std::to_string(st) + " (download failed or incomplete)";
            rc = 1;
        }
    }
    if (rc) fprintf(stderr, "[ERROR] %s: %s\n", s->path.c_str(), s->err.c_str());
    delete s;
    return rc;
}

static int seqsrc_kread(seqsrc_t *s, void *buf, int len)
{
    ssize_t k = seqsrc_read(s, buf, (size_t)len);
    return k < 0 ? 0 : (int)k;
}

KSEQ_INIT(seqsrc_t *, seqsrc_kread)

struct seqrec_s {
    seqsrc_t *src;
    kseq_t *ks;
};

seqrec_t *seqrec_open(seqsrc_t *src)
{
    seqrec_t *r = new seqrec_t();
    r->src = src;
    r->ks = kseq_init(src);
    return r;
}

int seqrec_next(seqrec_t *r, const char **name, const char **seq, int32_t *len)
{
    ssize_t l = kseq_read(r->ks);
    if (l < 0) {
        if (*seqsrc_error(r->src)) return -1;
        if (l == -2) { fprintf(stderr, "[ERROR] malformed FASTQ record (truncated quality string)\n"); return -1; }
        return 0;
    }
    char *nm = r->ks->name.s;
    if (nm) { char *t = strchr(nm, '\t'); if (t) *t = 0; }
    *name = nm ? nm : "";
    *seq = r->ks->seq.s;
    *len = (int32_t)l;
    return 1;
}

void seqrec_close(seqrec_t *r)
{
    if (!r) return;
    kseq_destroy(r->ks);
    delete r;
}

enum { SR_SAM = 0, SR_SAMGZ, SR_BAM, SR_CRAM, SR_URL };

#define SR_BUFSZ 65536

struct sr_t {
    int   fd;
    int   close_fd;
    pid_t child;
    seqsrc_t *src;
    int   fmt;

    uint8_t *in;
    size_t   in_cap, in_len, in_pos;
    int      raw_eof;
    int      phys_eof;

    int      bgzf;
    z_stream bz;
    int      bz_init;
    off_t    blk_next;         /* file offset of the next block */
    int      last_isize_zero;  /* last block decompressed to nothing = EOF marker */
    uint8_t *cbuf;
    size_t   cbuf_cap;

    int      indexed;
    int32_t  idx_tid;
    int32_t  last_refid;
    int      no_records;

    uint8_t *out;
    size_t   out_cap, out_len, out_pos;

    int      gz;
    z_stream zs;
    int      zs_init;
    int      mid_member;
    size_t   cur_member, last_member;

    int  eof;
    int  err;
    char errmsg[256];

    uint8_t *rec;    size_t rec_cap;
    char    *line;   size_t line_cap;
    char    *seqbuf; size_t seq_cap;

    /* BAM reference list */
    char  **refname;
    int32_t n_ref;

    char *region;
    int   filter_rname;
};

static char g_open_err[256];

static void sr_seterr(sr_t *r, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(r->errmsg, sizeof(r->errmsg), fmt, ap);
    va_end(ap);
    r->err = 1;
}

static inline int32_t sr_le32(const uint8_t *p)
{
    return (int32_t)((uint32_t)p[0] | (uint32_t)p[1] << 8 |
                     (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24);
}

static inline uint16_t sr_le16(const uint8_t *p)
{
    return (uint16_t)((uint32_t)p[0] | (uint32_t)p[1] << 8);
}

static int sr_grow(void *pbuf, size_t *cap, size_t need)
{
    uint8_t **buf = (uint8_t **)pbuf;
    if (*cap >= need) return 0;
    size_t nc = *cap ? *cap : 1024;
    while (nc < need) nc <<= 1;
    uint8_t *p = (uint8_t *)realloc(*buf, nc);
    if (!p) return -1;
    *buf = p;
    *cap = nc;
    return 0;
}

static int sr_raw_refill(sr_t *r)
{
    if (r->raw_eof) return 0;
    if (r->in_pos) {
        memmove(r->in, r->in + r->in_pos, r->in_len - r->in_pos);
        r->in_len -= r->in_pos;
        r->in_pos = 0;
    }
    if (r->in_len == r->in_cap) return 0;
    if (r->src) {
        ssize_t n = seqsrc_read(r->src, r->in + r->in_len, r->in_cap - r->in_len);
        if (n < 0) { sr_seterr(r, "%s", seqsrc_error(r->src)); return -1; }
        if (n == 0) { r->raw_eof = 1; return 0; }
        r->in_len += (size_t)n;
        return 1;
    }
    for (;;) {
        ssize_t n = read(r->fd, r->in + r->in_len, r->in_cap - r->in_len);
        if (n < 0) {
            if (errno == EINTR) continue;
            sr_seterr(r, "read failed: %s", strerror(errno));
            return -1;
        }
        if (n == 0) { r->raw_eof = 1; return 0; }
        r->in_len += (size_t)n;
        return 1;
    }
}

static ssize_t sr_readfd(sr_t *r, void *buf, size_t n)
{
    uint8_t *p = (uint8_t *)buf;
    size_t got = 0;
    while (got < n) {
        ssize_t k = read(r->fd, p + got, n - got);
        if (k < 0) {
            if (errno == EINTR) continue;
            sr_seterr(r, "read failed: %s", strerror(errno));
            return -1;
        }
        if (k == 0) break;
        got += (size_t)k;
    }
    if (got == 0) return 0;
    if (got < n) { sr_seterr(r, "input ends mid-block (file is truncated)"); return -1; }
    return (ssize_t)n;
}

static int sr_bgzf_block(sr_t *r)
{
    uint8_t hdr[12];
    ssize_t got = sr_readfd(r, hdr, 12);
    if (got <= 0) return got == 0 ? 0 : -1;

    if (hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 8 || !(hdr[3] & 4)) {
        sr_seterr(r, "not a BGZF block (BAM is expected to be BGZF-compressed)");
        return -1;
    }
    uint16_t xlen = sr_le16(hdr + 10);
    if (sr_grow(&r->cbuf, &r->cbuf_cap, xlen ? xlen : 1) < 0) { sr_seterr(r, "out of memory"); return -1; }
    if (xlen && sr_readfd(r, r->cbuf, xlen) <= 0) {
        if (!r->err) sr_seterr(r, "truncated BGZF header");
        return -1;
    }

    int have_bsize = 0;
    uint16_t bsize = 0;
    for (uint32_t i = 0; i + 4 <= xlen; ) {
        uint8_t si1 = r->cbuf[i], si2 = r->cbuf[i + 1];
        uint16_t slen = sr_le16(r->cbuf + i + 2);
        if (i + 4 + slen > xlen) break;
        if (si1 == 66 && si2 == 67 && slen == 2) { bsize = sr_le16(r->cbuf + i + 4); have_bsize = 1; break; }
        i += 4 + slen;
    }
    if (!have_bsize) { sr_seterr(r, "BGZF block has no BC size field"); return -1; }

    if ((size_t)bsize + 1 < (size_t)12 + xlen + 8) { sr_seterr(r, "bad BGZF block size"); return -1; }
    size_t cdata_len = (size_t)bsize + 1 - 12 - xlen - 8;

    if (sr_grow(&r->cbuf, &r->cbuf_cap, cdata_len + 8) < 0) { sr_seterr(r, "out of memory"); return -1; }
    if (cdata_len && sr_readfd(r, r->cbuf, cdata_len) <= 0) {
        if (!r->err) sr_seterr(r, "truncated BGZF payload");
        return -1;
    }
    uint8_t trailer[8];
    if (sr_readfd(r, trailer, 8) <= 0) { if (!r->err) sr_seterr(r, "truncated BGZF trailer"); return -1; }
    uint32_t isize = (uint32_t)sr_le32(trailer + 4);

    r->blk_next += (off_t)bsize + 1;
    r->last_isize_zero = (isize == 0);
    if (isize == 0) return 1;

    if (isize > r->out_cap - r->out_len) {
        if (sr_grow(&r->out, &r->out_cap, r->out_len + isize) < 0) { sr_seterr(r, "out of memory"); return -1; }
    }
    inflateReset(&r->bz);
    r->bz.next_in   = r->cbuf;
    r->bz.avail_in  = (uInt)cdata_len;
    r->bz.next_out  = r->out + r->out_len;
    r->bz.avail_out = (uInt)isize;
    int ret = inflate(&r->bz, Z_FINISH);
    if (ret != Z_STREAM_END || r->bz.avail_out != 0) {
        sr_seterr(r, "corrupt BGZF block: %s", r->bz.msg ? r->bz.msg : zError(ret));
        return -1;
    }
    r->out_len += isize;
    return 1;
}

static int sr_fill(sr_t *r)
{
    if (r->err) return -1;
    if (r->out_pos) {
        memmove(r->out, r->out + r->out_pos, r->out_len - r->out_pos);
        r->out_len -= r->out_pos;
        r->out_pos = 0;
    }
    if (r->out_len == r->out_cap) return 1;

    if (r->bgzf) {
        for (;;) {
            int g = sr_bgzf_block(r);
            if (g < 0) return -1;
            if (g == 0) { r->phys_eof = 1; return 0; }
            if (r->out_len > r->out_pos) return 1;   /* skip empty blocks */
        }
    }

    if (!r->gz) {
        if (r->in_pos == r->in_len) {
            if (sr_raw_refill(r) < 0) return -1;
            if (r->in_pos == r->in_len) { r->phys_eof = 1; return 0; }
        }
        size_t k = r->in_len - r->in_pos;
        if (k > r->out_cap - r->out_len) k = r->out_cap - r->out_len;
        memcpy(r->out + r->out_len, r->in + r->in_pos, k);
        r->in_pos += k;
        r->out_len += k;
        return 1;
    }

    for (;;) {
        if (r->in_pos == r->in_len) {
            if (sr_raw_refill(r) < 0) return -1;
            if (r->in_pos == r->in_len) {
                if (r->mid_member) {
                    sr_seterr(r, "input ends inside a compressed block (file is truncated)");
                    return -1;
                }
                r->phys_eof = 1;
                return 0;
            }
        }
        r->zs.next_in   = r->in + r->in_pos;
        r->zs.avail_in  = (uInt)(r->in_len - r->in_pos);
        r->zs.next_out  = r->out + r->out_len;
        r->zs.avail_out = (uInt)(r->out_cap - r->out_len);

        uInt ai = r->zs.avail_in, ao = r->zs.avail_out;
        int ret = inflate(&r->zs, Z_NO_FLUSH);
        size_t consumed = ai - r->zs.avail_in;
        size_t produced = ao - r->zs.avail_out;

        r->in_pos    += consumed;
        r->out_len   += produced;
        r->cur_member += produced;
        if (consumed || produced) r->mid_member = 1;

        if (ret == Z_STREAM_END) {
            r->last_member = r->cur_member;
            r->cur_member  = 0;
            r->mid_member  = 0;
            inflateReset(&r->zs);
            if (produced) return 1;
            continue;
        }
        if (ret != Z_OK) {
            sr_seterr(r, "zlib: %s", r->zs.msg ? r->zs.msg : zError(ret));
            return -1;
        }
        if (produced) return 1;
    }
}

static size_t sr_avail(sr_t *r, size_t n)
{
    while (r->out_len - r->out_pos < n) {
        int g = sr_fill(r);
        if (g < 0) return (size_t)-1;
        if (g == 0) break;
    }
    return r->out_len - r->out_pos;
}

static int sr_readn(sr_t *r, void *dst, size_t n)
{
    uint8_t *d = (uint8_t *)dst;
    size_t got = 0;
    while (got < n) {
        if (r->out_pos == r->out_len) {
            int g = sr_fill(r);
            if (g < 0) return -1;
            if (g == 0) break;
        }
        size_t k = r->out_len - r->out_pos;
        if (k > n - got) k = n - got;
        memcpy(d + got, r->out + r->out_pos, k);
        r->out_pos += k;
        got += k;
    }
    if (got == 0) return 0;
    if (got < n) {
        sr_seterr(r, "input ends mid-record (file is truncated)");
        return -1;
    }
    return 1;
}

static int sr_skipn(sr_t *r, size_t n)
{
    uint8_t scratch[4096];
    while (n) {
        size_t k = n > sizeof(scratch) ? sizeof(scratch) : n;
        int g = sr_readn(r, scratch, k);
        if (g <= 0) {
            if (g == 0) sr_seterr(r, "input ends inside the header (file is truncated)");
            return -1;
        }
        n -= k;
    }
    return 1;
}

static int sr_getline(sr_t *r, size_t *len_out)
{
    size_t len = 0;
    int any = 0;
    for (;;) {
        if (r->out_pos == r->out_len) {
            int g = sr_fill(r);
            if (g < 0) return -1;
            if (g == 0) { if (!any) return 0; break; }
        }
        any = 1;
        uint8_t *p = r->out + r->out_pos;
        size_t avail = r->out_len - r->out_pos;
        uint8_t *nl = (uint8_t *)memchr(p, '\n', avail);
        size_t take = nl ? (size_t)(nl - p) : avail;
        if (sr_grow(&r->line, &r->line_cap, len + take + 1) < 0) {
            sr_seterr(r, "out of memory");
            return -1;
        }
        memcpy(r->line + len, p, take);
        len += take;
        r->out_pos += take + (nl ? 1 : 0);
        if (nl) break;
    }
    if (len && r->line[len - 1] == '\r') --len;
    r->line[len] = 0;
    if (len_out) *len_out = len;
    return 1;
}

static void sr_bam_aux(const uint8_t *p, const uint8_t *end, char *hp, char *ps)
{
    memcpy(hp, "NA", 3);
    memcpy(ps, "NA", 3);
    while (p + 3 <= end) {
        char t0 = (char)p[0], t1 = (char)p[1];
        uint8_t type = p[2];
        const uint8_t *v = p + 3;
        long long ival = 0;
        int is_int = 0;
        size_t vlen;

        switch (type) {
        case 'A': vlen = 1; break;
        case 'c': if (v + 1 > end) return; ival = (int8_t)v[0];            is_int = 1; vlen = 1; break;
        case 'C': if (v + 1 > end) return; ival = v[0];                    is_int = 1; vlen = 1; break;
        case 's': if (v + 2 > end) return; ival = (int16_t)sr_le16(v);     is_int = 1; vlen = 2; break;
        case 'S': if (v + 2 > end) return; ival = sr_le16(v);              is_int = 1; vlen = 2; break;
        case 'i': if (v + 4 > end) return; ival = sr_le32(v);              is_int = 1; vlen = 4; break;
        case 'I': if (v + 4 > end) return; ival = (uint32_t)sr_le32(v);    is_int = 1; vlen = 4; break;
        case 'f': vlen = 4; break;
        case 'd': vlen = 8; break;
        case 'Z': case 'H': {
            const uint8_t *q = v;
            while (q < end && *q) ++q;
            if (q >= end) return;
            vlen = (size_t)(q - v) + 1;
            break;
        }
        case 'B': {
            if (v + 5 > end) return;
            uint8_t st = v[0];
            int32_t cnt = sr_le32(v + 1);
            size_t sz = (st == 'c' || st == 'C') ? 1 : (st == 's' || st == 'S') ? 2 : 4;
            if (cnt < 0) return;
            vlen = 5 + (size_t)cnt * sz;
            break;
        }
        default: return;
        }
        if (v + vlen > end) return;

        if (is_int) {
            if      (t0 == 'H' && t1 == 'P') snprintf(hp, 32, "%lld", ival);
            else if (t0 == 'P' && t1 == 'S') snprintf(ps, 32, "%lld", ival);
        }
        p = v + vlen;
    }
}

static int sr_bam_header(sr_t *r)
{
    uint8_t buf[4];

    if (sr_readn(r, buf, 4) != 1) {
        if (!r->err) sr_seterr(r, "truncated BAM magic");
        return -1;
    }
    if (sr_readn(r, buf, 4) != 1) { if (!r->err) sr_seterr(r, "truncated BAM header"); return -1; }
    int32_t l_text = sr_le32(buf);
    if (l_text < 0) { sr_seterr(r, "bad BAM header text length"); return -1; }
    if (l_text && sr_skipn(r, (size_t)l_text) < 0) return -1;

    if (sr_readn(r, buf, 4) != 1) { if (!r->err) sr_seterr(r, "truncated BAM header"); return -1; }
    int32_t n_ref = sr_le32(buf);
    if (n_ref < 0) { sr_seterr(r, "bad BAM reference count"); return -1; }

    r->refname = (char **)calloc((size_t)n_ref ? (size_t)n_ref : 1, sizeof(char *));
    if (!r->refname) { sr_seterr(r, "out of memory"); return -1; }
    r->n_ref = n_ref;

    for (int32_t i = 0; i < n_ref; ++i) {
        if (sr_readn(r, buf, 4) != 1) { if (!r->err) sr_seterr(r, "truncated BAM reference list"); return -1; }
        int32_t l_name = sr_le32(buf);
        if (l_name <= 0) { sr_seterr(r, "bad BAM reference name length"); return -1; }
        char *nm = (char *)malloc((size_t)l_name);
        if (!nm) { sr_seterr(r, "out of memory"); return -1; }
        if (sr_readn(r, nm, (size_t)l_name) != 1) {
            free(nm);
            if (!r->err) sr_seterr(r, "truncated BAM reference name");
            return -1;
        }
        nm[l_name - 1] = 0;
        r->refname[i] = nm;
        if (sr_readn(r, buf, 4) != 1) { if (!r->err) sr_seterr(r, "truncated BAM reference length"); return -1; }
    }
    return 0;
}

static const char sr_nt16[] = "=ACMGRSVTWYHKDBN";

static int sr_read_bam(sr_t *r, sr_rec_t *rec)
{
    uint8_t bs[4];
    int g = sr_readn(r, bs, 4);
    if (g <= 0) return g;

    int32_t block_size = sr_le32(bs);
    if (block_size < 32) { sr_seterr(r, "bad BAM record size %d", block_size); return -1; }
    if (sr_grow(&r->rec, &r->rec_cap, (size_t)block_size) < 0) { sr_seterr(r, "out of memory"); return -1; }
    if (sr_readn(r, r->rec, (size_t)block_size) != 1) {
        if (!r->err) sr_seterr(r, "input ends mid-record (file is truncated)");
        return -1;
    }

    const uint8_t *p = r->rec, *end = r->rec + block_size;
    int32_t  refID       = sr_le32(p);
    uint8_t  l_read_name = p[8];
    uint8_t  mapq        = p[9];
    uint16_t n_cigar     = sr_le16(p + 12);
    uint16_t flag        = sr_le16(p + 14);
    int32_t  l_seq       = sr_le32(p + 16);

    if (l_read_name == 0 || l_seq < 0) { sr_seterr(r, "malformed BAM record"); return -1; }

    const uint8_t *q = p + 32;
    if (q + l_read_name > end) { sr_seterr(r, "malformed BAM record (name)"); return -1; }
    const char *qname = (const char *)q;
    q += l_read_name;

    if ((size_t)(end - q) < 4 * (size_t)n_cigar) { sr_seterr(r, "malformed BAM record (cigar)"); return -1; }
    q += 4 * (size_t)n_cigar;

    size_t nseq = ((size_t)l_seq + 1) / 2;
    if ((size_t)(end - q) < nseq + (size_t)l_seq) { sr_seterr(r, "malformed BAM record (seq)"); return -1; }
    const uint8_t *seqp = q;
    q += nseq + (size_t)l_seq;

    if (sr_grow(&r->seqbuf, &r->seq_cap, (size_t)l_seq + 1) < 0) { sr_seterr(r, "out of memory"); return -1; }
    for (int32_t i = 0; i < l_seq; ++i)
        r->seqbuf[i] = sr_nt16[(seqp[i >> 1] >> ((~i & 1) << 2)) & 0xf];
    r->seqbuf[l_seq] = 0;

    sr_bam_aux(q, end, rec->hp, rec->ps);

    r->last_refid = refID;
    rec->qname = qname;
    rec->rname = (refID >= 0 && refID < r->n_ref) ? r->refname[refID] : "*";
    rec->flag  = flag;
    rec->mapq  = mapq;
    rec->l_seq = l_seq;
    rec->seq   = r->seqbuf;
    return 1;
}

static int sr_read_sam(sr_t *r, sr_rec_t *rec)
{
    for (;;) {
        size_t len;
        int g = sr_getline(r, &len);
        if (g <= 0) return g;
        if (len == 0 || r->line[0] == '@') continue;

        char *fields[11];
        char *p = r->line, *aux = NULL;
        int nf;
        for (nf = 0; nf < 11; ++nf) {
            fields[nf] = p;
            char *t = strchr(p, '\t');
            if (t) { *t = 0; p = t + 1; }
            else   { p = NULL; break; }
        }
        if (nf < 10) { sr_seterr(r, "malformed SAM line (%d fields)", nf + 1); return -1; }
        aux = p;

        rec->qname = fields[0];
        rec->flag  = (uint16_t)strtoul(fields[1], NULL, 10);
        rec->rname = fields[2];
        rec->mapq  = (uint8_t)strtoul(fields[4], NULL, 10);

        char *seq = fields[9];
        if (seq[0] == '*' && seq[1] == 0) {
            seq[0] = 0;
            rec->l_seq = 0;
        } else {
            size_t n = strlen(seq);
            for (size_t i = 0; i < n; ++i) seq[i] = (char)toupper((unsigned char)seq[i]);
            rec->l_seq = (int32_t)n;
        }
        rec->seq = seq;

        memcpy(rec->hp, "NA", 3);
        memcpy(rec->ps, "NA", 3);
        for (char *a = aux; a && *a; ) {
            char *t = strchr(a, '\t');
            if (t) *t = 0;
            if (a[0] && a[1] && a[2] == ':' && a[3] && a[4] == ':') {
                if      (a[0] == 'H' && a[1] == 'P' && a[3] == 'i') snprintf(rec->hp, 32, "%s", a + 5);
                else if (a[0] == 'P' && a[1] == 'S' && a[3] == 'i') snprintf(rec->ps, 32, "%s", a + 5);
            }
            a = t ? t + 1 : NULL;
        }
        return 1;
    }
}

static int sr_bgzf_seek(sr_t *r, uint64_t voff)
{
    off_t    blk    = (off_t)(voff >> 16);
    uint32_t within = (uint32_t)(voff & 0xffff);

    if (lseek(r->fd, blk, SEEK_SET) < 0) {
        sr_seterr(r, "seek failed: %s", strerror(errno));
        return -1;
    }
    r->blk_next  = blk;
    r->out_len   = r->out_pos = 0;
    r->raw_eof   = 0;
    r->phys_eof  = 0;

    int g = sr_fill(r);
    if (g < 0) return -1;
    if (g == 0 || within > r->out_len) {
        sr_seterr(r, "index points past the end of the file (is the .bai stale?)");
        return -1;
    }
    r->out_pos = within;
    return 0;
}

static int sr_bai_start(const char *path, int32_t tid, uint64_t *voff)
{
    char buf1[4096], buf2[4096];
    const char *cand[2];
    int nc = 0;

    snprintf(buf1, sizeof(buf1), "%s.bai", path);
    cand[nc++] = buf1;
    size_t pl = strlen(path);
    if (pl > 4 && strcmp(path + pl - 4, ".bam") == 0 && pl < sizeof(buf2)) {
        memcpy(buf2, path, pl - 4);
        memcpy(buf2 + pl - 4, ".bai", 5);
        cand[nc++] = buf2;
    }

    uint8_t *d = NULL;
    size_t dl = 0;
    for (int i = 0; i < nc && !d; ++i) {
        FILE *f = fopen(cand[i], "rb");
        if (!f) continue;
        if (fseek(f, 0, SEEK_END) == 0) {
            long sz = ftell(f);
            if (sz > 8 && fseek(f, 0, SEEK_SET) == 0) {
                uint8_t *p = (uint8_t *)malloc((size_t)sz);
                if (p && fread(p, 1, (size_t)sz, f) == (size_t)sz) { d = p; dl = (size_t)sz; }
                else free(p);
            }
        }
        fclose(f);
    }
    if (!d) return -1;

    size_t p = 0;
    if (dl < 8 || memcmp(d, "BAI\1", 4) != 0) goto bad;
    p = 4;
    int32_t n_ref;
    n_ref = sr_le32(d + p); p += 4;
    if (tid >= n_ref) goto bad;

    for (int32_t i = 0; i < n_ref; ++i) {
        if (p + 4 > dl) goto bad;
        int32_t n_bin = sr_le32(d + p); p += 4;
        if (n_bin < 0) goto bad;
        for (int32_t b = 0; b < n_bin; ++b) {
            if (p + 8 > dl) goto bad;
            p += 4;
            int32_t n_chunk = sr_le32(d + p); p += 4;
            if (n_chunk < 0 || (size_t)n_chunk * 16 > dl - p) goto bad;
            p += (size_t)n_chunk * 16;
        }
        if (p + 4 > dl) goto bad;
        int32_t n_intv = sr_le32(d + p); p += 4;
        if (n_intv < 0 || (size_t)n_intv * 8 > dl - p) goto bad;

        if (i == tid) {
            uint64_t best = (uint64_t)-1;
            for (int32_t k = 0; k < n_intv; ++k) {
                const uint8_t *e = d + p + (size_t)k * 8;
                uint64_t v = (uint64_t)(uint32_t)sr_le32(e) |
                             ((uint64_t)(uint32_t)sr_le32(e + 4) << 32);
                if (v && v < best) best = v;
            }
            free(d);
            *voff = best;
            return 0;
        }
        p += (size_t)n_intv * 8;
    }
bad:
    free(d);
    return -1;
}

static int sr_spawn_samtools(sr_t *r, const char *path, const char *ref,
                             const char *region, int n_threads)
{
    int pfd[2];
    if (pipe(pfd) != 0) {
        snprintf(g_open_err, sizeof(g_open_err), "pipe failed: %s", strerror(errno));
        return -1;
    }

    char thr[16];
    snprintf(thr, sizeof(thr), "%d", n_threads > 1 ? n_threads : 1);

    const char *argv[12];
    int n = 0;
    argv[n++] = "samtools";
    argv[n++] = "view";
    argv[n++] = "-u";     /* uncompressed BAM */
    if (n_threads > 1)   { argv[n++] = "-@"; argv[n++] = thr; }
    if (ref && *ref)     { argv[n++] = "-T"; argv[n++] = ref; }
    argv[n++] = path;
    if (region && *region) argv[n++] = region;
    argv[n] = NULL;

    signal(SIGPIPE, SIG_IGN);

    pid_t pid = fork();
    if (pid < 0) {
        close(pfd[0]); close(pfd[1]);
        snprintf(g_open_err, sizeof(g_open_err), "fork failed: %s", strerror(errno));
        return -1;
    }
    if (pid == 0) {
        close(pfd[0]);
        if (dup2(pfd[1], STDOUT_FILENO) < 0) _exit(127);
        close(pfd[1]);
        execvp("samtools", (char *const *)argv);
        _exit(127);   /* samtools not on PATH */
    }
    close(pfd[1]);
    r->fd       = pfd[0];
    r->close_fd = 1;
    r->child    = pid;
    return 0;
}

const char *sr_open_error(void) { return g_open_err; }

int sr_path_is_url(const char *path)
{
    if (!path) return 0;
    static const char *schemes[] = {
        "http://", "https://", "ftp://",
        "s3://", "s3+http://", "s3+https://",
        "gs://", "gs+http://", "gs+https://",
    };
    for (size_t i = 0; i < sizeof schemes / sizeof *schemes; ++i)
        if (strncmp(path, schemes[i], strlen(schemes[i])) == 0) return 1;
    return 0;
}

static int sr_child_probe(sr_t *r, const char *path, int is_url)
{
    while (r->in_len == r->in_pos && !r->raw_eof) {
        if (sr_raw_refill(r) < 0) {
            snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg);
            return -1;
        }
    }
    if (r->in_len > r->in_pos) return 0;

    int st = 0;
    while (waitpid(r->child, &st, 0) < 0 && errno == EINTR) { }
    r->child = -1;
    int code = WIFEXITED(st) ? WEXITSTATUS(st) : -1;
    if (code == 0) {
        r->no_records = 1;
        return 0;
    }
    if (code == 127) {
        snprintf(g_open_err, sizeof(g_open_err),
                 is_url ? "could not run samtools, which is required for URL input; "
                          "install it or download the file and pass a local path"
                        : "could not run samtools, which is required for CRAM input; "
                          "install it or convert to BAM (samtools view -b)");
    } else if (is_url) {
        snprintf(g_open_err, sizeof(g_open_err),
                 "samtools could not read %s (exit %d); check that the URL is reachable "
                 "and, if a chromosome was requested, that its index (.bai/.crai) "
                 "sits beside it on the server", path, code);
    } else {
        snprintf(g_open_err, sizeof(g_open_err),
                 "samtools could not read %s (exit %d); check the reference given "
                 "to -r and that a .crai index exists", path, code);
    }
    return -1;
}

static int sr_child_bam_setup(sr_t *r, const char *path)
{
    r->zs.zalloc = Z_NULL; r->zs.zfree = Z_NULL; r->zs.opaque = Z_NULL;
    if (inflateInit2(&r->zs, 15 + 16) != Z_OK) {
        snprintf(g_open_err, sizeof(g_open_err), "inflateInit2 failed");
        return -1;
    }
    r->zs_init = 1;
    r->gz = 1;

    size_t avail = sr_avail(r, 4);
    if (avail == (size_t)-1) { snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg); return -1; }
    if (avail < 4 || memcmp(r->out + r->out_pos, "BAM\1", 4) != 0) {
        snprintf(g_open_err, sizeof(g_open_err),
                 "samtools view -u did not produce BAM for %s; the stream is unreadable", path);
        return -1;
    }
    if (sr_bam_header(r) < 0) { snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg); return -1; }
    return 0;
}

static int sr_is_sam_header(const uint8_t *m, size_t n)
{
    static const char *rt[5] = { "HD", "SQ", "RG", "PG", "CO" };
    if (n < 4 || m[0] != '@' || m[3] != '\t') return 0;
    for (int i = 0; i < 5; ++i)
        if (memcmp(m + 1, rt[i], 2) == 0) return 1;
    return 0;
}

/* classify decompressed leading bytes. */
static int sr_classify(const uint8_t *m, size_t n)
{
    if (n >= 4 && memcmp(m, "BAM\1", 4) == 0) return SR_FMT_BAM;
    if (n >= 4 && memcmp(m, "CRAM", 4) == 0) return SR_FMT_CRAM;
    if (n >= 1 && m[0] == '>')  return SR_FMT_FASTA;
    if (n >= 1 && m[0] == '@')  return sr_is_sam_header(m, n) ? SR_FMT_SAM : SR_FMT_FASTQ;
    return SR_FMT_SAM;
}

int sr_classify_bytes(const void *m, size_t n) { return sr_classify((const uint8_t *)m, n); }

int sr_probe_format(const char *path)
{
    if (!path || strcmp(path, "-") == 0) return SR_FMT_UNKNOWN;

    gzFile fp = gzopen(path, "rb");
    if (!fp) return SR_FMT_UNKNOWN;
    uint8_t buf[16];
    int n = gzread(fp, buf, (unsigned)sizeof buf);
    gzclose(fp);
    if (n <= 0) return SR_FMT_UNKNOWN;
    return sr_classify(buf, (size_t)n);
}

const char *sr_fmt_name(int fmt)
{
    switch (fmt) {
    case SR_FMT_SAM:   return "sam";
    case SR_FMT_BAM:   return "bam";
    case SR_FMT_CRAM:  return "cram";
    case SR_FMT_FASTA: return "fasta";
    case SR_FMT_FASTQ: return "fastq";
    case SR_FMT_URL:   return "url";
    default:           return "unknown";
    }
}

int sr_fmt_parse(const char *name)
{
    if (!name) return SR_FMT_UNKNOWN;
    if (strcmp(name, "sam")   == 0) return SR_FMT_SAM;
    if (strcmp(name, "bam")   == 0) return SR_FMT_BAM;
    if (strcmp(name, "cram")  == 0) return SR_FMT_CRAM;
    if (strcmp(name, "fasta") == 0 || strcmp(name, "fa") == 0) return SR_FMT_FASTA;
    if (strcmp(name, "fastq") == 0 || strcmp(name, "fq") == 0) return SR_FMT_FASTQ;
    return SR_FMT_UNKNOWN;
}

int sr_fmt_is_alignment(int fmt)
{
    return fmt == SR_FMT_SAM || fmt == SR_FMT_BAM || fmt == SR_FMT_CRAM ||
           fmt == SR_FMT_URL;
}

static int sr_reject_sequence(const uint8_t *m, size_t n, const char *path)
{
    int fmt = sr_classify(m, n);
    if (fmt != SR_FMT_FASTA && fmt != SR_FMT_FASTQ) return 0;
    snprintf(g_open_err, sizeof(g_open_err),
             "%s is %s, not alignments; pass --input-fmt %s, or omit "
             "--input-fmt to detect the format from the file's contents",
             strcmp(path, "-") == 0 ? "stdin" : path,
             sr_fmt_name(fmt), sr_fmt_name(fmt));
    return 1;
}

static sr_t *sr_alloc(void)
{
    sr_t *r = (sr_t *)calloc(1, sizeof(sr_t));
    if (!r) { snprintf(g_open_err, sizeof(g_open_err), "out of memory"); return NULL; }
    r->fd = -1;
    r->child = -1;
    r->in_cap = r->out_cap = SR_BUFSZ;
    r->in  = (uint8_t *)malloc(r->in_cap);
    r->out = (uint8_t *)malloc(r->out_cap);
    if (!r->in || !r->out) {
        snprintf(g_open_err, sizeof(g_open_err), "out of memory");
        free(r->in); free(r->out); free(r);
        return NULL;
    }
    return r;
}

static void sr_free(sr_t *r);

sr_t *sr_open_source(seqsrc_t *src, const char *label)
{
    g_open_err[0] = 0;
    sr_t *r = sr_alloc();
    if (!r) return NULL;
    r->src = src;

    while (r->in_len - r->in_pos < 4 && !r->raw_eof) {
        if (sr_raw_refill(r) < 0) {
            snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg);
            sr_free(r);
            return NULL;
        }
    }
    size_t got = r->in_len - r->in_pos;
    const uint8_t *m = r->in + r->in_pos;
    if (got >= 4 && memcmp(m, "BAM\1", 4) == 0) {
        r->fmt = SR_BAM;
        if (sr_bam_header(r) < 0) {
            snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg);
            sr_free(r);
            return NULL;
        }
    } else if (got >= 4 && memcmp(m, "CRAM", 4) == 0) {
        snprintf(g_open_err, sizeof(g_open_err), "%s is CRAM, which the reads engine cannot read; "
                 "convert it with samtools view -b (or samtools fastq)", label);
        sr_free(r);
        return NULL;
    } else {
        if (sr_reject_sequence(m, got, label)) { sr_free(r); return NULL; }
        r->fmt = SR_SAM;
    }
    return r;
}

sr_t *sr_open(const char *path, const char *ref, const char *region, int n_threads)
{
    g_open_err[0] = 0;

    sr_t *r = sr_alloc();
    if (!r) return NULL;
    if (region && *region) {
        r->region = strdup(region);
        if (!r->region) { snprintf(g_open_err, sizeof(g_open_err), "out of memory"); goto fail; }
    }

    int is_stdin;
    is_stdin = (!path || strcmp(path, "-") == 0);

    if (!is_stdin && sr_path_is_url(path)) {
        if (sr_spawn_samtools(r, path, ref, region, n_threads) < 0) goto fail;
        r->fmt = SR_URL;
        r->filter_rname = 0;   /* samtools already applied the region */
        if (sr_child_probe(r, path, 1) < 0) goto fail;
        if (!r->no_records && sr_child_bam_setup(r, path) < 0) goto fail;
        return r;
    }

    if (is_stdin) {
        r->fd = STDIN_FILENO;
    } else {
        r->fd = open(path, O_RDONLY);
        if (r->fd < 0) {
            snprintf(g_open_err, sizeof(g_open_err), "cannot open %s: %s", path, strerror(errno));
            goto fail;
        }
        r->close_fd = 1;
    }

    /* sniff the format from the leading bytes. */
    while (r->in_len - r->in_pos < 4 && !r->raw_eof) {
        if (sr_raw_refill(r) < 0) {
            snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg);
            goto fail;
        }
    }
    {
        size_t got = r->in_len - r->in_pos;
        const uint8_t *m = r->in + r->in_pos;

        if (got >= 4 && memcmp(m, "CRAM", 4) == 0) {
            if (is_stdin) {
                snprintf(g_open_err, sizeof(g_open_err),
                         "CRAM cannot be read from stdin; pass the .cram path so samtools can index it");
                goto fail;
            }
            close(r->fd);
            r->fd = -1;
            r->close_fd = 0;
            r->in_len = r->in_pos = 0;
            r->raw_eof = 0;
            if (sr_spawn_samtools(r, path, ref, region, n_threads) < 0) goto fail;
            r->fmt = SR_CRAM;
            r->filter_rname = 0;
            if (sr_child_probe(r, path, 0) < 0) goto fail;
            if (!r->no_records && sr_child_bam_setup(r, path) < 0) goto fail;
            return r;
        }

        if (got >= 2 && m[0] == 0x1f && m[1] == 0x8b) {
            r->zs.zalloc = Z_NULL; r->zs.zfree = Z_NULL; r->zs.opaque = Z_NULL;
            if (inflateInit2(&r->zs, 15 + 16) != Z_OK) {
                snprintf(g_open_err, sizeof(g_open_err), "inflateInit2 failed");
                goto fail;
            }
            r->zs_init = 1;
            r->gz = 1;

            size_t avail = sr_avail(r, 4);
            if (avail == (size_t)-1) { snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg); goto fail; }
            if (avail >= 4 && memcmp(r->out + r->out_pos, "BAM\1", 4) == 0) {
                r->fmt = SR_BAM;

                if (!is_stdin && lseek(r->fd, 0, SEEK_SET) == 0) {
                    inflateEnd(&r->zs);
                    r->zs_init = 0;
                    r->gz = 0;
                    memset(&r->bz, 0, sizeof(r->bz));
                    if (inflateInit2(&r->bz, -15) != Z_OK) {
                        snprintf(g_open_err, sizeof(g_open_err), "inflateInit2 (raw) failed");
                        goto fail;
                    }
                    r->bz_init = 1;
                    r->bgzf = 1;
                    r->in_len = r->in_pos = r->out_len = r->out_pos = 0;
                    r->raw_eof = r->mid_member = 0;
                    r->cur_member = r->last_member = 0;
                    r->blk_next = 0;
                    if (sr_grow(&r->out, &r->out_cap, 2 * SR_BUFSZ) < 0) {
                        snprintf(g_open_err, sizeof(g_open_err), "out of memory");
                        goto fail;
                    }
                }

                if (sr_bam_header(r) < 0) { snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg); goto fail; }

                if (r->bgzf && r->region) {
                    int32_t tid = -1;
                    for (int32_t i = 0; i < r->n_ref; ++i)
                        if (strcmp(r->refname[i], r->region) == 0) { tid = i; break; }
                    if (tid >= 0) {
                        uint64_t voff;
                        if (sr_bai_start(path, tid, &voff) == 0) {
                            if (voff == (uint64_t)-1) {
                                r->no_records = 1;
                            } else if (sr_bgzf_seek(r, voff) < 0) {
                                snprintf(g_open_err, sizeof(g_open_err), "%s", r->errmsg);
                                goto fail;
                            }
                            r->indexed = 1;
                            r->idx_tid = tid;
                        } else {
                            fprintf(stderr, "[INFO] no .bai index for %s; scanning the whole file for %s\n", path, r->region);
                        }
                    }
                }
            } else {
                if (sr_reject_sequence((const uint8_t *)(r->out + r->out_pos), avail, path))
                    goto fail;
                r->fmt = SR_SAMGZ;
            }
        } else {
            if (sr_reject_sequence(m, got, path)) goto fail;
            r->fmt = SR_SAM;
        }
    }

    r->filter_rname = (r->region != NULL && !r->indexed);
    return r;

fail:
    sr_free(r);
    return NULL;
}

static void sr_free(sr_t *r)
{
    if (!r) return;
    if (r->zs_init) inflateEnd(&r->zs);
    if (r->bz_init) inflateEnd(&r->bz);
    if (r->close_fd && r->fd >= 0) close(r->fd);
    for (int32_t i = 0; i < r->n_ref; ++i) free(r->refname[i]);
    free(r->refname);
    free(r->region);
    free(r->rec);
    free(r->line);
    free(r->seqbuf);
    free(r->cbuf);
    free(r->in);
    free(r->out);
    free(r);
}

int sr_read(sr_t *r, sr_rec_t *rec)
{
    if (!r || r->err) return -1;
    if (r->no_records) { r->eof = 1; return 0; }
    for (;;) {
        int is_bam = (r->fmt == SR_BAM || r->fmt == SR_CRAM || r->fmt == SR_URL);
        int g = is_bam ? sr_read_bam(r, rec) : sr_read_sam(r, rec);
        if (g < 0) return -1;
        if (g == 0) { r->eof = 1; return 0; }

        if (r->indexed) {
            if (r->last_refid < 0 || r->last_refid > r->idx_tid) { r->eof = 1; return 0; }
            if (r->last_refid < r->idx_tid) continue;
        } else if (r->filter_rname && strcmp(rec->rname, r->region) != 0) {
            continue;
        }
        return 1;
    }
}

int sr_has_ref(const sr_t *r, const char *name)
{
    if (!r || !r->refname) return -1;
    for (int32_t i = 0; i < r->n_ref; ++i)
        if (strcmp(r->refname[i], name) == 0) return 1;
    return 0;
}

const char *sr_error(const sr_t *r) { return r ? r->errmsg : ""; }

const char *sr_format(const sr_t *r)
{
    if (!r) return "?";
    switch (r->fmt) {
    case SR_BAM:   return "BAM";
    case SR_CRAM:  return "CRAM";
    case SR_URL:   return "URL (via samtools)";
    case SR_SAMGZ: return "SAM/gz";
    default:       return "SAM";
    }
}

int sr_close(sr_t *r)
{
    if (!r) return 0;
    int had_err = r->err;
    int rc = had_err ? 1 : 0;

    if (!r->err && r->phys_eof &&
        (r->fmt == SR_BAM || ((r->fmt == SR_CRAM || r->fmt == SR_URL) && r->gz))) {
        int missing = r->bgzf ? !r->last_isize_zero : (r->mid_member || r->last_member != 0);
        if (missing) {
            sr_seterr(r, "BAM ends without a BGZF end-of-file marker (file is truncated)");
            rc = 1;
        }
    }

    if (r->child > 0) {
        if (r->fd >= 0) { close(r->fd); r->fd = -1; r->close_fd = 0; }
        int st = 0;
        while (waitpid(r->child, &st, 0) < 0 && errno == EINTR) { }
        if (WIFEXITED(st)) {
            int code = WEXITSTATUS(st);
            if (code == 127) {
                sr_seterr(r, "could not run samtools (needed for CRAM and URL input); is it on your PATH?");
                rc = 1;
            } else if (code != 0) {
                if (r->fmt == SR_URL)
                    sr_seterr(r, "samtools exited with status %d while streaming the URL; "
                                 "the output is incomplete (network error, or a bad URL/region)", code);
                else
                    sr_seterr(r, "samtools exited with status %d; CRAM decoding failed "
                                 "(check the reference passed to -r and that a .crai index exists)", code);
                rc = 1;
            }
        } else if (WIFSIGNALED(st)) {
            if (!(WTERMSIG(st) == SIGPIPE && !r->eof)) {
                sr_seterr(r, "samtools was killed by signal %d", WTERMSIG(st));
                rc = 1;
            }
        }
    }

    if (r->zs_init) inflateEnd(&r->zs);
    if (r->bz_init) inflateEnd(&r->bz);
    if (r->close_fd && r->fd >= 0) close(r->fd);
    for (int32_t i = 0; i < r->n_ref; ++i) free(r->refname[i]);
    free(r->refname);
    free(r->region);
    free(r->rec);
    free(r->line);
    free(r->seqbuf);
    free(r->cbuf);
    free(r->in);
    free(r->out);
    if (rc && !had_err) fprintf(stderr, "[ERROR] %s\n", r->errmsg);
    free(r);
    return rc;
}
