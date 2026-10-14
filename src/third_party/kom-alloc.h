/* PROVENANCE (tr-calling)
 *   Upstream:  lh3/kmer-map -- https://github.com/lh3/kmer-map
 *              file kmer-map.c, lines 14-24 (kom_malloc / kom_calloc /
 *              kom_realloc / kom_grow) and line 173 (ba_strdup).
 *   Author:    Heng Li (lh3)
 *   License:   lh3/kmer-map carries NO LICENSE file as of 2026-07-31.
 *              Reused here with the author's knowledge -- lh3 oversees the
 *              tr-calling project. If tr-calling is ever released publicly,
 *              ask lh3 to add an explicit LICENSE to kmer-map (every other
 *              lh3 repository is MIT), or replace this file: the macros
 *              below are a few lines of malloc wrappers and are trivially
 *              reimplementable.
 *   Local modifications:
 *     - kom_grow / kom_malloc / kom_calloc / kom_realloc: verbatim.
 *     - kom_strdup: adapted from kmer-map's ba_strdup(), which allocates
 *       from a block arena. This version uses plain malloc instead, since
 *       tr-calling frees these strings individually.
 *   Exact upstream commit was not recorded when this code was vendored.
 *
 *   Used by: src/detection.cpp, src/tr-detection.cpp,
 *            src/tr-detection-fa.cpp
 */

#ifndef TR_THIRD_PARTY_KOM_ALLOC_H
#define TR_THIRD_PARTY_KOM_ALLOC_H

#include <stdlib.h>
#include <string.h>

/* memory management macros */
#define kom_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define kom_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define kom_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

/* dynamic array growth macro */
#define kom_grow(type, ptr, __i, __m) do { \
        if ((__i) >= (__m)) { \
            (__m) = (__i) + 1; \
            (__m) += ((__m)>>1) + 16; \
            (ptr) = kom_realloc(type, (ptr), (__m)); \
        } \
    } while (0)

/* string duplication (strdup() doesn't conform to C99) */
static inline char *kom_strdup(const char *src)
{
    size_t len = strlen(src);
    char *dst = kom_malloc(char, len + 1);
    memcpy(dst, src, len + 1);
    return dst;
}

#endif /* TR_THIRD_PARTY_KOM_ALLOC_H */
