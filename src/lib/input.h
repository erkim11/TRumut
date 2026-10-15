/*
input.h: trumut input readers.
*/
#ifndef INPUT_H
#define INPUT_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { SEQSRC_RAW = 0, SEQSRC_GZIP = 1, SEQSRC_BGZF = 2 };

typedef struct seqsrc_s seqsrc_t;

seqsrc_t *seqsrc_open(const char *path, int inflate_threads, char *err, size_t err_n);

int seqsrc_url_kind(const char *path);

int seqsrc_kind(const seqsrc_t *s);

ssize_t seqsrc_read(seqsrc_t *s, void *buf, size_t n);

ssize_t seqsrc_peek(seqsrc_t *s, void *buf, size_t n);

const char *seqsrc_error(const seqsrc_t *s);

int seqsrc_close(seqsrc_t *s);

typedef struct seqrec_s seqrec_t;

seqrec_t *seqrec_open(seqsrc_t *src);

int seqrec_next(seqrec_t *r, const char **name, const char **seq, int32_t *len);

void seqrec_close(seqrec_t *r);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SR_FSECONDARY     0x100
#define SR_FSUPPLEMENTARY 0x800

typedef struct {
    const char *qname;
    const char *rname;
    uint16_t    flag;
    uint8_t     mapq;
    int32_t     l_seq;
    const char *seq;
    char        hp[32];  /* HP aux tag, or "NA" */
    char        ps[32];  /* PS aux tag, or "NA" */
} sr_rec_t;

typedef struct sr_t sr_t;

sr_t *sr_open(const char *path, const char *ref, const char *region, int n_threads);

struct seqsrc_s;
sr_t *sr_open_source(struct seqsrc_s *src, const char *label);

int sr_read(sr_t *r, sr_rec_t *rec);

int sr_close(sr_t *r);

const char *sr_error(const sr_t *r);

const char *sr_open_error(void);

int sr_has_ref(const sr_t *r, const char *name);

const char *sr_format(const sr_t *r);

// format probing
enum {
    SR_FMT_UNKNOWN = -1,
    SR_FMT_SAM = 0,
    SR_FMT_BAM,
    SR_FMT_CRAM,
    SR_FMT_FASTA,
    SR_FMT_FASTQ,
    SR_FMT_URL
};

int sr_path_is_url(const char *path);

/* detect a file's format from its leading bytes. */
int sr_probe_format(const char *path);

/* "sam", "bam", "cram", "fasta", "fastq", or "unknown". */
const char *sr_fmt_name(int fmt);

int sr_fmt_parse(const char *name);

int sr_classify_bytes(const void *m, size_t n);

int sr_fmt_is_alignment(int fmt);

#ifdef __cplusplus
}
#endif

#endif
