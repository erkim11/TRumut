/* The MIT License

   Copyright (c) 2008-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* PROVENANCE (tr-calling)
 *   Upstream:  lh3/seqtk -- https://github.com/lh3/seqtk
 *              file seqtk.c: comp_tab[] at line 226, and the in-place
 *              reverse-complement loop at lines 1511-1518.
 *   Author:    Heng Li (lh3)
 *   License:   MIT (notice above is seqtk's LICENSE, reproduced in full).
 *   Exact upstream commit was not recorded when this code was vendored.
 *
 *   Local modifications to comp_tab -> kom_comp_table:
 *     1. EXTENDED from 128 to 256 entries. Entries 128-255 are identity.
 *        seqtk's table is only 128 long, so indexing it with a high-bit
 *        byte reads out of bounds; the 256-entry version cannot.
 *     2. BUG FIX at index 96 (backtick). seqtk has 64 there, which looks
 *        like a copy-paste slip -- its row for 96-111 was duplicated from
 *        the row for 64-79, whose leading entry legitimately is 64. This
 *        version maps 96 -> 96 (identity), consistent with every other
 *        non-nucleotide byte in the table. DO NOT "fix" this back to 64.
 *     All other entries in 0-127 are identical to seqtk's.
 *
 *   Local modifications to the revcomp loop -> kom_revcomp():
 *     - Extracted from seqtk's inline loop into a standalone function
 *       taking (uint64_t len, char *seq) instead of operating on a kseq_t.
 *       Algorithm is unchanged: swap-and-complement inward from both ends,
 *       then complement the middle base when the length is odd.
 *
 *   Used by: include/detection.h (re-export), src/detection.cpp
 */

#ifndef TR_THIRD_PARTY_SEQTK_REVCOMP_H
#define TR_THIRD_PARTY_SEQTK_REVCOMP_H

#include <stdint.h>

// complement table for reverse complement
inline uint8_t kom_comp_table[256] = {
      0,   1,    2,   3,    4,   5,    6,   7,    8,   9,  10,  11,   12,  13,  14,  15,
     16,  17,   18,  19,   20,  21,   22,  23,   24,  25,  26,  27,   28,  29,  30,  31,
     32,  33,   34,  35,   36,  37,   38,  39,   40,  41,  42,  43,   44,  45,  46,  47,
     48,  49,   50,  51,   52,  53,   54,  55,   56,  57,  58,  59,   60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

// reverse complement a sequence in place
inline void kom_revcomp(uint64_t len, char *seq)
{
    uint64_t i;
    for (i = 0; i < len>>1; ++i) {
        uint8_t t = (uint8_t)seq[len - i - 1];
        seq[len - i - 1] = (char)kom_comp_table[(uint8_t)seq[i]];
        seq[i] = (char)kom_comp_table[t];
    }
    if (len&1) seq[len>>1] = (char)kom_comp_table[(uint8_t)seq[len>>1]];
}

#endif /* TR_THIRD_PARTY_SEQTK_REVCOMP_H */
