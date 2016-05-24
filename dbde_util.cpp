/** DBDE Util: converts Dynamic Bit Depth Encoding videos *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences   *
  * This file distributed under the Apache License 2.0    **/

#include <x86intrin.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "dbde_util.h"

/************************************************
 * Routines to pack image data into DBDE format *
 ************************************************/

uint32_t dbde_pack_8x8(uint8_t *image, int stride, uint8_t *target) {
    // Pack all the data into SSE data structures
    __m128i ab = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i cd = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i ef = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)); image += 2*stride;
    __m128i gh = _mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image));

    // Compute min/max across the whole array (16 wide)
    __m128i hi = _mm_max_epu8(ab, cd);
    __m128i lo = _mm_min_epu8(ab, cd);
    hi = _mm_max_epu8(hi, ef);
    lo = _mm_min_epu8(lo, ef);
    hi = _mm_max_epu8(hi, gh);
    lo = _mm_min_epu8(lo, gh);

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

    if ((I1 & 0x80) != 0) {
        // All bits required
        _mm_storeu_si128((__m128i*)target, ab); target += sizeof(__m128i);
        _mm_storeu_si128((__m128i*)target, cd); target += sizeof(__m128i);
        _mm_storeu_si128((__m128i*)target, ef); target += sizeof(__m128i);
        _mm_storeu_si128((__m128i*)target, gh);
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

        x = _mm_set1_epi32((1 << (2*k+16)) + 1);
        ab = _mm_madd_epi16(x, ab);                // Now 4 ints across
        cd = _mm_madd_epi16(x, cd);                // Now 4 ints across
        ef = _mm_madd_epi16(x, ef);                // Now 4 ints across
        gh = _mm_madd_epi16(x, gh);                // Now 4 ints across


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
                if (kk >= 64) {
                    *((uint64_t*)target) = p;
                    target += sizeof(uint64_t);
                    p = (k > 64) ? ((*i) >> (96-kk)) : 0;
                    kk -= 64;
                }
                i++;
            }
        }
        return (k << 6) | I0;  // Note we already multiplied k by 4, so we shift two less!
    }
}

/** MUST have rightmargin >=  1 && downmargin >= 1, and must also have at least one of rightmargin < 8 && downmargin < 8! */
uint32_t dbde_pack_8x8_partial(uint8_t *image, int stride, int rightmargin, int downmargin, uint8_t *target) {
    uint8_t full[64];
    if (rightmargin >= 8) {
        int j = 0;
        for (; j < downmargin; j++) *(((uint64_t*)full) + j) = *((uint64_t*)(image + j*stride));
        for (; j < 8; j++) *(((uint64_t*)full) + j) = *((uint64_t*)(image + (downmargin-1)*stride));
    }
    else {
        uint8_t *fu = full;
        int j = 0;
        for (; j < downmargin; j++) {
            int k = 0;
            for (; k < rightmargin; k++) {
                *fu = *(image + k);
                fu += 1;
            }
            uint8_t c = *(fu-1);
            for (; k < 8; k++) {
                *fu = c;
                fu += 1;
            }
            image += stride;
        }
        if (j < 8) {
            uint64_t c8 = *(((uint64_t*)full) + (j-1));
            for (; j < 8; j++) *(((uint64_t*)full) + j) = c8;
        }
    }
    return dbde_pack_8x8(full, 8, target);
}

size_t dbde_pack_image(uint8_t *image, int W, int H, uint8_t *target) {
    int w = (W+7)/8;
    int h = (H+7)/8;
    int sz = 2*(4 + w*h)+4;
    *((int*)target) = w*h; target += 4;
    uint8_t *bd = target; target += w*h;
    *((int*)target) = w*h; target += 4;
    uint8_t *mi = target; target += w*h;
    int *n64 = (int*)target; target += 4;
    *n64 = 0;
    int ww = W/8;
    int hh = H/8;
    int y = 0;
    for (; y < hh; y++) {
        int x = 0;
        for (; x < ww; x++) {
            uint32_t bdmi = dbde_pack_8x8(image + 8*(x + W*y), W, target);
            uint8_t bi = bdmi >> 8;
            target = target + ((int)bi)*8;
            *(bd++) = bi;
            *(mi++) = (uint8_t)(bdmi & 0xFF);
            *n64 += bi;
        }
        if (x < w) {
            uint32_t bdmi = dbde_pack_8x8_partial(image + 8*(x + W*y), W, W&0x7, 8, target);
            uint8_t bi = bdmi >> 8;
            target = target + ((int)bi)*8;
            *(bd++) = bi;
            *(mi++) = (uint8_t)(bdmi & 0xFF);
            *n64 += bi;
        }
    }
    if (y < h) {
        for (int x = 0; x < w; x++) {
            uint32_t bdmi = dbde_pack_8x8_partial(image + 8*(x + W*y), W, (x==ww) ? (W&0x7) : 8, H&0x7, target);
            uint8_t bi = bdmi >> 8;
            target = target + ((int)bi)*8;
            *(bd++) = bi;
            *(mi++) = (uint8_t)(bdmi & 0xFF);
            *n64 += bi;
        }
    }
    return sz + 8 * (*n64);
}

size_t dbde_pack_frame_header(frame_header fh, uint8_t *target) {
    *((frame_header*)target) = fh;
    return sizeof(frame_header);
}

size_t dbde_pack_frame(uint64_t index, uint8_t *image, int W, int H, uint8_t *target) {
    frame_header fh = { 2, 0ll, 0ll };
    fh.index = index;
    size_t sz = dbde_pack_frame_header(fh, target);
    sz += dbde_pack_image(image, W, H, target + sz);
    return sz;
}

size_t dbde_pack_video_header(video_header vh, uint8_t *target) {
    *((video_header*)target) = vh;
    return sizeof(video_header);
}


/**********************************************
 * Routines to unpack DBDE data into an image *
 **********************************************/

void dbde_unpack_8x8(uint8_t depth, uint8_t minval, uint8_t* packed, size_t stride, uint8_t *image) {
    __m128i lo = _mm_set1_epi8(minval);
    if (depth == 0) {
        printf("Depth is 0\n");
        int64_t l = _mm_cvtsi128_si64(lo);
        for (int j = 0; j < 8; j++) { *((int64_t*)image) = l; image += stride; }        
    }
    else {
        uint8_t temp[64];
        if (depth < 8) {
            uint32_t *data = (uint32_t*)packed;
            uint64_t v = *(data++);
            int m = (1 << depth) - 1;
            int nb = 32;
            for (int i = 0; i < 64; i++) {
                if (nb < depth) {
                    v = v | ((*(data++)) << nb);
                    nb += 32;
                }
                temp[i] = (uint8_t)(v & m);
                v = v >> depth;
                nb -= depth;
            }
            packed = temp;
        }
        __m128i x;
        x = _mm_add_epi8(lo, _mm_loadu_si128((__m128i*)packed)); packed += sizeof(__m128i);
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
        x = _mm_add_epi8(lo, _mm_loadu_si128((__m128i*)packed)); packed += sizeof(__m128i);
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
        x = _mm_add_epi8(lo, _mm_loadu_si128((__m128i*)packed)); packed += sizeof(__m128i);
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
        x = _mm_add_epi8(lo, _mm_loadu_si128((__m128i*)packed));
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1);      
    }
}

void dbde_unpack_8x8_partial(uint8_t depth, uint8_t minval, uint8_t* packed, size_t stride, int rightmargin, int downmargin, uint8_t* image) {
    uint8_t img[64];
    dbde_unpack_8x8(depth, minval, packed, 8, img);
    for (int y = 0; y < downmargin; y++) {
        for (int x = 0; x < rightmargin; x++) {
            image[stride*y + x] = img[8*y + x];
        }
    }
}

size_t dbde_unpack_image(uint8_t *packed, int W, int H, uint8_t *image) {
    uint8_t *pack = packed;
    int w = (W+7)/8;
    int h = (H+7)/8;
    int32_t nb = *((int32_t*)pack); pack += sizeof(int32_t);
    if (nb != w*h) return 0;
    uint8_t *b = pack; pack += nb;
    int32_t nm = *((int32_t*)pack); pack += sizeof(int32_t);
    if (nm != w*h) return 0;
    uint8_t *m = pack; pack += nm;
    int32_t n64 = *((int32_t*)pack); pack += sizeof(int32_t);
    for (int i=0; i < w*h; i++) n64 -= b[i];
    if (n64 != 0) return 0;
    int ww = W/8;
    int hh = H/8;
    int y = 0;
    for (; y < hh; y++) {
        int x = 0;
        for (; x < ww; x++) {
            uint8_t bi = *(b++);
            dbde_unpack_8x8(bi, *(m++), pack, W, image + 8*(y*W + x));
            pack += 8*((int)bi);
        }
        if (x < w) {
            uint8_t bi = *(b++);
            dbde_unpack_8x8_partial(bi, *(m++), pack, W, W & 0x7, 8, image + 8*(y*W + x));
            pack += 8*((int)bi);
        }
    }
    if (y < h) {
        for (int x = 0; x < w; x++) {
            uint8_t bi = *(b++);
            dbde_unpack_8x8_partial(bi, *(m++), pack, W, (x==ww) ? (W & 0x7) : 8, H & 0x7, image + 8*(y*W + x));
            pack += 8*((int)bi);
        }
    }
    return pack - packed;
}

frame_header dbde_unpack_frame(uint8_t **packed, int W, int H, uint8_t *image) {
    frame_header fh = *((frame_header*)*packed);
    if (fh.u64s != 2) { fh.u64s = -1; return fh; }
    *packed += sizeof(frame_header);
    size_t n = dbde_unpack_image(*packed, W, H, image);
    if (n == 0) { fh.u64s = -1; }
    else { *packed += n; }
    return fh;
}

video_header dbde_unpack_frame(uint8_t **packed) {
    video_header vh = *((video_header*)*packed);
    if (vh.u64s != 3) { vh.u64s = -1; }
    else { *packed += sizeof(video_header); }
    return vh;
}

#ifdef DBDE_UTIL_MAIN

int main(int argc, char **argv) {
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
    uint8_t out[64], img[64];
    memset(out, 0, 64);

    int64_t ta0 = __rdtsc();
    int i = dbde_pack_8x8(example, 10, out);
    dbde_unpack_8x8((i >> 8), (i & 0xFF), out, 8, img);
    int64_t tz0 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf(
            "%016lx          %016lx          %016lx\n",
            *((int64_t*)(out + 8*j)),
            *((int64_t*)(img + 8*j)),
            *((int64_t*)(example + 10*j))
        );
    }

    printf("\n\n\n");

    int64_t ta1 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8, 10, 2, 8, out);
    dbde_unpack_8x8((i >> 8), (i & 0xFF), out, 8, img);
    int64_t tz1 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        int64_t ex = 0;
        for (int k = 0; k < 8; k++) ex = ex | (((uint64_t)example[8 + 10*j + (k > 1 ? 1 : k)]) << (8*k));
        printf(
            "%016lx          %016lx          %016lx\n",
            *((int64_t*)(out + 8*j)),
            *((int64_t*)(img + 8*j)),
            ex
        );
    }

    printf("\n\n\n");

    int64_t ta2 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8*10, 10, 8, 2, out);
    dbde_unpack_8x8((i >> 8), (i & 0xFF), out, 8, img);
    int64_t tz2 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        printf(
            "%016lx          %016lx          %016lx\n",
            *((int64_t*)(out + 8*j)),
            *((int64_t*)(img + 8*j)),
            *((int64_t*)(example + 80 + 10*(j > 1 ? 1 : j)))
        );
    }

    printf("\n\n\n");

    int64_t ta3 = __rdtsc();
    i = dbde_pack_8x8_partial(example + 8*10 + 8, 10, 2, 2, out);
    dbde_unpack_8x8((i >> 8), (i & 0xFF), out, 8, img);
    int64_t tz3 = __rdtsc();
    printf("%x\n", i);
    for (int j = 0; j < 8; j++) { 
        int64_t ex = 0;
        for (int k = 0; k < 8; k++) ex = ex | (((uint64_t)example[88 + 10*(j > 1 ? 1 : j) + (k > 1 ? 1 : k)]) << (8*k));
        printf(
            "%016lx          %016lx          %016lx\n",
            *((int64_t*)(out + 8*j)),
            *((int64_t*)(img + 8*j)),
            ex
        );
    }

    printf("\n\n\n");

#define XRES 2536
#define YRES 2048

#ifdef PIECEWISE
    uint8_t *big = (uint8_t*)malloc(XRES*YRES + (XRES+YRES));
    uint8_t *targ = (uint8_t*)malloc(8*(XRES*YRES + (XRES*YRES)/4 + XRES+YRES));
    for (int bi = 0; bi < XRES*YRES + (XRES+YRES); bi++) big[bi] = rand();
    uint32_t *ix = ((uint32_t*)&(out[0]));
    uint32_t si = 0;
    uint32_t n = 0;
    int64_t ta4 = __rdtsc();
    for (int nj = 0; nj < YRES/8; nj++) {
        for (int ni = 0; ni < XRES/8; ni++) {
            dbde_pack_8x8(&(big[nj*XRES*8 + ni*8]), XRES, out);
            si += *ix + *(ix+1) + *(ix+2) + *(ix+3) + *(ix+4) + *(ix+5) + *(ix+6) + *(ix+7);
            n += 1;
        }
    }
    int64_t tz4 = __rdtsc();
    int64_t ta5 = tz4;
    int64_t tz5 = ta5;
    printf("%d %d\n", si, n);
#else
    uint8_t *big = (uint8_t*)malloc(XRES*YRES + (XRES+YRES));
    uint8_t *targ = (uint8_t*)malloc(XRES*YRES + (XRES*YRES)/16 + 12 + XRES+YRES + 20);
    uint8_t *re = (uint8_t*)malloc(XRES*YRES + (XRES+YRES));
    for (int bi = 0; bi < XRES*YRES + (XRES+YRES); bi++) big[bi] = rand() & 0x7F;
    uint32_t si = 0;
    uint32_t n = 0;
    int64_t ta4 = __rdtsc();
    int nX = XRES;
    int nY = YRES;
    si = dbde_pack_frame(1, big, nX, nY, targ);
    n = *((int*)(targ + si/2));
    int64_t tz4 = __rdtsc();
    int64_t ta5 = tz4;
    uint8_t *targv = targ;
    dbde_unpack_frame(&targv, nX, nY, re);
    printf("%ld %ld\n",(long)si ,targv - targ);
    int64_t tz5 = __rdtsc();
    int wrong = 0;
    for (int y = 0; y < nY; y++)
        for (int x = 0; x < nX; x++)
            if (big[y*nX+x] != re[y*nX+x]) wrong++;
    printf("%x %d   !! %d of %d\n", si, n, wrong, nX*nY);
    n = ((nX+7)/8) * ((nY+7)/8);
#endif    

#undef YRES
#undef XRES

    printf(
        "\n\n%d %d %d %d\n%f %f\n%f %f\n",
        (int)(tz0-ta0),
        (int)(tz1-ta1),
        (int)(tz2-ta2),
        (int)(tz3-ta3),
        (tz4-ta4)/((double)n),
        1.0/((tz4-ta4)*3e-10),
        (tz5-ta5)/((double)n),
        1.0/((tz5-ta5)*3e-10)
    );
}

#endif
