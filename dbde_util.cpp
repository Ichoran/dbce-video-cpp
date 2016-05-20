/** DBDE Util: converts Dynamic Bit Depth Encoding videos *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences   *
  * This file distributed under the Apache License 2.0    **/

#include <x86intrin.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "dbde_util.h"

//#define LOG_ME

void p(__m128i x) {
    printf("%016lx\n%016lx\n\n", *((int64_t*)&x), *((int64_t*)(((char*)&x)+sizeof(int64_t))));
}

int dbde_pack_8x8(uint8_t *image, int stride, uint8_t *target) {
    // Pack all the data into SSE data structures
    __m128i ab = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i cd = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i ef = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i gh = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image));

#ifdef LOG_ME
    printf("STARTING\n");
    p(ab);
    p(cd);
    p(ef);
    p(gh);
#endif

    // Compute min/max across the whole array (16 wide)
    __m128i hi = _mm_max_epu8(ab, cd);
    __m128i lo = _mm_min_epu8(ab, cd);
    hi = _mm_max_epu8(hi, ef);
    lo = _mm_min_epu8(lo, ef);
    hi = _mm_max_epu8(hi, gh);
    lo = _mm_min_epu8(lo, gh);

#ifdef LOG_ME
    printf("HI/LO INITIAL\n");
    p(hi);
    p(lo);
#endif

    // For min: reduce to 8 wide, then use minpos instruction (just works)
    lo = _mm_min_epu8(lo, _mm_srli_epi16(lo, 8)); // stuvwxyz min 0s0u0x0z
    lo = _mm_minpos_epu16(lo);                    // Answer is in low byte of low int16_t (high bytes are all 0 from shift)
    int I0 = _mm_cvtsi128_si32(lo) & 0xFF;        // Probably faster than _mm_extract_epi8(lo, 0) & 0xFF

    // For max: convert to min by subtracting from 255, then use min technique
    hi = _mm_sub_epi8(_mm_set1_epi8(-1), hi);
    hi = _mm_min_epu8(hi, _mm_srli_epi16(hi, 8));
    hi = _mm_minpos_epu16(hi);
    int I1 = 255 - (_mm_extract_epi16(hi, 0) & 0xFF);

    if (I0 == I1) return I0;

    I1 -= I0;
    lo = _mm_set1_epi8(I0);
    ab = _mm_sub_epi8(ab, lo);                 // Packed 16 across
    cd = _mm_sub_epi8(cd, lo);                 // Packed 16 across
    ef = _mm_sub_epi8(ef, lo);                 // Packed 16 across
    gh = _mm_sub_epi8(gh, lo);                 // Packed 16 across

#ifdef LOG_ME
    printf("\nI0 = %08x I1 = %08x\nSUBTRACTED\n", I0, I1);
    p(lo);
    p(ab);
    p(cd);
    p(ef);
    p(gh);
#endif

    if ((I1 & 0x80) != 0) {
        // All bits required
        *((__m128i*)target) = ab; target += sizeof(__m128i);
        *((__m128i*)target) = cd; target += sizeof(__m128i);
        *((__m128i*)target) = ef; target += sizeof(__m128i);
        *((__m128i*)target) = gh;
        return 0x800 | I0;
    }
    else {
        int k = (I1 >= 0x10) ?
                    (((I1&0x40) == 0) ? (((I1&0x20) == 0) ? 5 : 6) : 7) :
                    (((I1&0x0C)==0) ? (((I1&0x02)==0) ? 1 : 2) : (((I1&0x08)==0) ? 3 : 4));
        __m128i x;
        x = _mm_set1_epi16((1 << (k+8)) + 1);
        ab = _mm_maddubs_epi16(x, ab);             // Now 8 shorts across
        cd = _mm_maddubs_epi16(x, cd);             // Now 8 shorts across
        ef = _mm_maddubs_epi16(x, ef);             // Now 8 shorts across
        gh = _mm_maddubs_epi16(x, gh);             // Now 8 shorts across

#ifdef LOG_ME
        printf("\n8SHORTS\n"); p(ab); p(cd); p(ef); p(gh);
#endif

        x = _mm_set1_epi32((1 << (2*k+16)) + 1);
        ab = _mm_madd_epi16(x, ab);                // Now 4 ints across
        cd = _mm_madd_epi16(x, cd);                // Now 4 ints across
        ef = _mm_madd_epi16(x, ef);                // Now 4 ints across
        gh = _mm_madd_epi16(x, gh);                // Now 4 ints across

#ifdef LOG_ME
        printf("\n4INTS\n"); p(ab); p(cd); p(ef); p(gh);
#endif

        int *i;
        int kk = 0;
        uint64_t p = 0;
        k = k << 2;
        __m128i *packs[] = { &ab, &cd, &ef, &gh };
        for (int m = 0; m < 4; m++) {
            i = (int*)packs[m];
            for (int mm = 0; mm < 4; mm++) {
                p = p | (((uint64_t)(*i)) << kk);
                kk += k;
#ifdef LOG_ME
                printf("%016lx %08x\n", p, kk);
#endif
                if (kk >= 64) {
                    *((uint64_t*)target) = p;
                    target += sizeof(uint64_t);
                    p = (k > 64) ? ((*i) >> (96-kk)) : 0;
                    kk -= 64;
                }
                i++;
            }
        }
#ifdef LOG_ME
        printf("%x\n",I1);
#endif
        return (k << 14) | I0;  // Note we already multiplied k by 4, so we shift two less!
    }
}

/** MUST have rightmargin >=  1 && downmargin >= 1, and must also have at least one of rightmargin < 8 && downmargin < 8! */
int dbde_pack_8x8_partial(uint8_t *image, int stride, int rightmargin, int downmargin, uint8_t *target) {
    uint8_t full[64];
    if (rightmargin >= 8) {
        int j = 0;
        for (; j < downmargin; j++) *(((uint64_t*)full) + j) = *((uint64_t*)(image + j*stride));
        for (; j < 8; j++) *(((uint64_t*)full) + j) = *((uint64_t*)(image + (downmargin-1)*stride));
    }
    else {
#ifdef LOG_ME
        printf("DOWN %d RIGHT %d\n", downmargin, rightmargin);
#endif
        uint8_t *fu = full;
        int j = 0;
        for (; j < downmargin; j++) {
            int k = 0;
            for (; k < rightmargin; k++) {
                *fu = *(image + k);
#ifdef LOG_ME
                printf("%02x", *fu);
#endif
                fu += 1;
            }
            uint8_t c = *(fu-1);
            for (; k < 8; k++) {
                *fu = c;
                fu += 1;
            }
            image += stride;
#ifdef LOG_ME
            printf("\n");
#endif
        }
        if (j < 8) {
            uint64_t c8 = *(((uint64_t*)full) + (j-1));
            for (; j < 8; j++) *(((uint64_t*)full) + j) = c8;
        }
    }
#ifdef LOG_ME
    uint64_t *fx = (uint64_t*)full;
    for (int fi = 0; fi < 8; fi++) printf("%016lx\n", *(fx++));
#endif
    return dbde_pack_8x8(full, 8, target);
}

video_header dbde_unpack_video_header(uint8_t *encoded) {
    video_header vh = *((video_header*)encoded);
    return vh;
}


int dbde_pack_video_header(video_header vh, uint8_t *encoded) {
    video_header* target = (video_header*)encoded;
    *target = vh;
    return sizeof(video_header);
}


frame_header dbde_unpack_frame_header(uint8_t *encoded, int length) {
    frame_header fh;
    if (length >= sizeof(frame_header)) fh = *((frame_header*)encoded);
    else { fh.u64s = 0; }
    return fh;
}

int dbde_pack_frame_header(frame_header fh, uint8_t *encoded) {
    frame_header* target = (frame_header*)encoded;
    *target = fh;
    return sizeof(frame_header);
}

dbde_data dbde_unpack_data(uint8_t *encoded) {
    dbde_data dd;
    dd.bits_len = *((int32_t*)encoded);
    encoded += sizeof(int32_t);
    dd.bits = encoded;
    encoded += dd.bits_len;
    dd.mins_len = *((int32_t*)encoded);
    encoded += sizeof(int32_t);
    dd.mins = encoded;
    encoded += dd.mins_len;
    dd.data_len = *((int32_t*)encoded);
    encoded += dd.data_len * sizeof(uint64_t);
    return dd;
}

int dbde_pack_data(dbde_data dd, uint8_t *encoded) {
    *((int32_t*)encoded) = dd.bits_len;
    encoded += sizeof(int32_t);
    memcpy(encoded, dd.bits, dd.bits_len);
    encoded += dd.bits_len;
    *((int32_t*)encoded) = dd.mins_len;
    encoded += sizeof(int32_t);
    memcpy(encoded, dd.mins, dd.mins_len);
    encoded += dd.mins_len;
    *((int32_t*)encoded) = dd.data_len;
    encoded += sizeof(int32_t);
    memcpy(encoded, dd.data, sizeof(uint64_t)*dd.data_len);
    return 3*sizeof(uint32_t) + dd.bits_len + dd.mins_len * sizeof(uint64_t)*dd.data_len;
}

/*
void dbde_decode_image(dbde_data dd, int W, int H, uint8_t *image);

uint8_t* dbde_decode_image_alloc(dbde_data dd, int W, int H);

void dbde_encode_image(uint8_t *image, int W, int H, dbde_data& dd);

dbde_data* dbde_encode_image_alloc(uint8_t *image, int W, int H);
*/

int main(int argc, char **argv) {
    uint8_t data[] = { 4, 2, 3, 4, 5, 6, 7, 8,
                       9,10,11,15,13,14, 2,12,
                       4, 4, 3, 7, 9,15, 2, 6,
                      11,14,33, 3, 4,15, 4, 15,
                      15,14,13,12,11,10, 9, 8,
                       7, 6, 5, 4, 3, 3, 3, 3,
                      22,18, 4,16, 8, 5,12,11,
                      14,16,14,16,14,16,14,16
                     };

    uint8_t example[] = 
        {   25, 27, 23, 29, 22, 24, 29, 23, 25, 24,
            22, 24, 21, 25, 22, 27, 28, 21, 27, 26,
            25, 26, 22, 29, 25, 20, 28, 23, 26, 25,
            19, 23, 25, 21, 28, 19, 22, 25, 25, 27,
            27, 25, 30, 28, 25, 23, 27, 26, 24, 24,
            31, 30, 31, 28, 29, 26, 24, 25, 27, 26,
            30, 28, 32, 25, 28, 27, 28, 27, 26, 26,
            29, 31, 31, 32, 29, 29, 25, 22, 24, 25,
            31, 34, 33, 31, 30, 29, 28, 28, 26, 26,
            34, 34, 35, 35, 33, 28, 29, 28, 26, 26
        };
    uint8_t out[64];
    memset(out, 0, 64);

    int64_t ta0 = __rdtsc();
    int i = dbde_pack_8x8(example, 10, out);
    int64_t tz0 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf("%016lx\n", *((int64_t*)(out + 8*j)));
    }

    printf("\n\n\n");

    int64_t ta1 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8, 10, 2, 8, out);
    int64_t tz1 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf("%016lx\n", *((int64_t*)(out + 8*j)));
    }

    printf("\n\n\n");

    int64_t ta2 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8*10, 10, 8, 2, out);
    int64_t tz2 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf("%016lx\n", *((int64_t*)(out + 8*j)));
    }

    printf("\n\n\n");

    int64_t ta3 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8*10 + 8, 10, 2, 2, out);
    int64_t tz3 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf("%016lx\n", *((int64_t*)(out + 8*j)));
    }

    uint32_t *big = (uint32_t*)calloc(1024*1024 + 2048, 4);
    for (int bi = 0; bi < 1024*1024 + 2048; bi++) big[bi] = rand();
    uint32_t *ix = ((uint32_t*)&(out[0]));
    uint32_t si = 0;
    uint32_t n = 0;
    int64_t ta4 = __rdtsc();
    for (int ni = 0; ni < 1024*1024/16; ni++) {
        dbde_pack_8x8(((uint8_t*)&(big[ni*16])), 2048, out);
        si += *ix + *(ix+1) + *(ix+2) + *(ix+3) + *(ix+4) + *(ix+5) + *(ix+6) + *(ix+7);
        n += 1;
    }
    int64_t tz4 = __rdtsc();
    printf("%d %d\n", si, n);

    printf("\n\n%d %d %d %d\n%f\n", (int)(tz0-ta0), (int)(tz1-ta1), (int)(tz2-ta2), (int)(tz3-ta3), (tz4-ta4)/((double)n));
}
