/* PROVENANCE (tr-calling)
 *   Upstream:  lh3/kmer-map -- https://github.com/lh3/kmer-map
 *              file kmer-map.c, line 122 (kom_nt4_table).
 *   Author:    Heng Li (lh3)
 *   License:   lh3/kmer-map carries NO LICENSE file as of 2026-07-31.
 *              Reused here with the author's knowledge -- lh3 oversees the
 *              tr-calling project. See third_party/kom-alloc.h for the same
 *              note; if tr-calling is released publicly, ask lh3 to add an
 *              explicit LICENSE to kmer-map.
 *   Local modifications: none. Byte-identical to upstream.
 *   Exact upstream commit was not recorded when this code was vendored.
 *
 *   NOTE: this table is NOT the same as minimap2's seq_nt4_table, which
 *   additionally maps 'U'/'u' (85/117) to 3 for uracil. This one maps them
 *   to 4 (unknown). Do not swap one for the other.
 *
 *   NOTE: indices 0-3 map to 0-3 rather than to 4, so the table is
 *   idempotent on already-encoded input.
 *
 *   Used by: include/detection.h (re-export), src/detection.cpp
 */

#ifndef TR_THIRD_PARTY_NT4_TABLE_H
#define TR_THIRD_PARTY_NT4_TABLE_H

#include <stdint.h>

// nt4 encoding table: A/C/G/T -> 0/1/2/3; others -> 4
inline uint8_t kom_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

#endif /* TR_THIRD_PARTY_NT4_TABLE_H */
