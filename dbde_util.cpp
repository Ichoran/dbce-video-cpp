/** DBDE Util: converts Dynamic Bit Depth Encoding videos *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences   *
  * This file distributed under the Apache License 2.0    **/

#include <smmintrin.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "dbde_util.h"

/************************************************
 * Routines to pack image data into DBDE format *
 ************************************************/

#ifdef DBDE_INVERT_ENDIAN
#define ENDIAN(x) _mm_shuffle_epi8(x, _mm_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8))
#else
#define ENDIAN(x) x
#endif


uint32_t dbde_pack_8x8(uint8_t *image, int stride, uint8_t *target) {
    // Pack all the data into SSE data structures
    __m128i ab = ENDIAN(_mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image))); image += 2*stride;
    __m128i cd = ENDIAN(_mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image))); image += 2*stride;
    __m128i ef = ENDIAN(_mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image))); image += 2*stride;
    __m128i gh = ENDIAN(_mm_set_epi64x(*((int64_t*)(image + stride)), *((int64_t*)image)));

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

        uint32_t *i;
        int kk = 0;
        uint64_t p = 0;
        k = k << 2;
        __m128i *packs[] = { &ab, &cd, &ef, &gh };
        for (int m = 0; m < 4; m++) {
            i = (uint32_t*)packs[m];
            for (int mm = 0; mm < 4; mm++) {
                p = p | (((uint64_t)(*i)) << kk);
                kk += k;
                if (kk >= 32) {
                    *((uint32_t*)target) = (uint32_t)p;
                    target += sizeof(uint32_t);
                    p = (kk > 32) ? (p >> 32) : 0;
                    kk -= 32;
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
    size_t offset = 0;
    *((int32_t*)target) = fh.u64s; offset += 4;
    *((uint64_t*)(target + offset)) = fh.index; offset += 8;
    *((uint64_t*)(target + offset)) = fh.reserved0; offset += 8;
    return offset;
}

size_t dbde_pack_frame(uint64_t index, uint8_t *image, int W, int H, uint8_t *target) {
    frame_header fh = { 2, 0ll, 0ll };
    fh.index = index;
    size_t sz = dbde_pack_frame_header(fh, target);
    sz += dbde_pack_image(image, W, H, target + sz);
    return sz;
}

size_t dbde_pack_video_header(video_header vh, uint8_t *target) {
    size_t offset = 0;
    *((int32_t*)target) = vh.u64s; offset += 4;
    *((uint64_t*)(target + offset)) = vh.height; offset += 8;
    *((uint64_t*)(target + offset)) = vh.width; offset += 8;
    *((uint64_t*)(target + offset)) = vh.frame_hz; offset += 8;
    return offset;
}


/**********************************************
 * Routines to unpack DBDE data into an image *
 **********************************************/

void dbde_unpack_8x8(uint8_t depth, uint8_t minval, uint8_t* packed, size_t stride, uint8_t *image) {
    __m128i lo = _mm_set1_epi8(minval);
    if (depth == 0) {
#ifdef __x86_64__
        int64_t l = _mm_cvtsi128_si64(lo);
        for (int j = 0; j < 8; j++) { *((int64_t*)image) = l; image += stride; }        
#else
        int32_t l = _mm_cvtsi128_si32(lo);
        for (int j = 0; j < 8; j++) { *((int32_t*)image) = l; *(((int32_t*)image)+1) = l; image += stride; }
#endif
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
                    v = v | (((uint64_t)(*(data++))) << nb);
                    nb += 32;
                }
                temp[i] = (uint8_t)(v & m);
                v = v >> depth;
                nb -= depth;
            }
            packed = temp;
        }
        __m128i x;
        x = _mm_add_epi8(lo, ENDIAN(_mm_loadu_si128((__m128i*)packed))); packed += sizeof(__m128i);
#ifdef __x86_64__ 
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
#else
        *((int32_t*)image) = _mm_cvtsi128_si32(x); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 1); image += stride;
        *((int32_t*)image) = _mm_extract_epi32(x, 2); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 3); image += stride;
#endif
        x = _mm_add_epi8(lo, ENDIAN(_mm_loadu_si128((__m128i*)packed))); packed += sizeof(__m128i);
#ifdef __x86_64__ 
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
#else
        *((int32_t*)image) = _mm_cvtsi128_si32(x); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 1); image += stride;
        *((int32_t*)image) = _mm_extract_epi32(x, 2); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 3); image += stride;
#endif
        x = _mm_add_epi8(lo, ENDIAN(_mm_loadu_si128((__m128i*)packed))); packed += sizeof(__m128i);
#ifdef __x86_64__ 
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1); image += stride;
#else
        *((int32_t*)image) = _mm_cvtsi128_si32(x); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 1); image += stride;
        *((int32_t*)image) = _mm_extract_epi32(x, 2); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 3); image += stride;
#endif
        x = _mm_add_epi8(lo, ENDIAN(_mm_loadu_si128((__m128i*)packed)));
#ifdef __x86_64__ 
        *((int64_t*)image) = _mm_cvtsi128_si64(x); image += stride;
        *((int64_t*)image) = _mm_extract_epi64(x, 1);
#else
        *((int32_t*)image) = _mm_cvtsi128_si32(x); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 1); image += stride;
        *((int32_t*)image) = _mm_extract_epi32(x, 2); *(((int32_t*)image)+1) = _mm_extract_epi32(x, 3);
#endif
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

frame_header dbde_unpack_frame_header(uint8_t **packed) {
    frame_header fh;
    fh.u64s = *((int*)*packed); *packed += 4;
    fh.index = *((uint64_t*)*packed); *packed += 8;
    fh.reserved0 = *((uint64_t*)*packed); *packed += 8;
    if (fh.u64s != 2) fh.u64s = -1;
    return fh;
}

frame_header dbde_unpack_frame(uint8_t **packed, int W, int H, uint8_t *image) {
    frame_header fh = dbde_unpack_frame_header(packed);
    size_t n = dbde_unpack_image(*packed, W, H, image);
    if (n == 0) { fh.u64s = -1; }
    else { *packed += n; }
    return fh;
}

video_header dbde_unpack_video_header(uint8_t **packed) {
    video_header vh;
    vh.u64s = *((int*)*packed); *packed += 4;
    vh.height = *((uint64_t*)*packed); *packed += 8;
    vh.width = *((uint64_t*)*packed); *packed += 8;
    vh.frame_hz = *((uint64_t*)*packed); *packed += 8;
    if (vh.u64s != 3) { vh.u64s = -1; }
    return vh;
}


dbde_file_walker dbde_start_file_walk(const char* name, int frames_buffered, video_header *vh) {
    if (frames_buffered < 1) frames_buffered = 2;
    uint8_t tinybuf[28];   // Just the header!
    FILE *f = fopen(name, "rb");
    if (!f) return (dbde_file_walker){NULL, 0, 0, 0, 0, 1, 1, NULL};
    dbde_file_walker w = (dbde_file_walker){f, 0, 0, 0, 0, 1, 1, NULL};
    int n = fread(tinybuf, 1, 28, w.fptr);
    uint8_t *buf = tinybuf;
    *vh = dbde_unpack_video_header(&buf);
    if (vh->u64s != 3 || (buf - tinybuf) != 28) { fclose(w.fptr); w.fptr = NULL; return w; }
    uint64_t npix = ((vh->height * vh->width) * frames_buffered);
    w.N = npix + (npix/8) + ((uint64_t)frames_buffered)*32;
    if (vh->height == 0 || vh->width == 0 ||
        vh->height > 0x37FFFFFF || vh->width > 0x37FFFFFF ||
        (vh->height * vh->width) > 0x37FFFFFF ||
        w.N >= 0x7FFFFFFF
    ) {
        fclose(w.fptr);
        w.fptr = NULL;
        return w;
    }
    w.buffer = (uint8_t*)malloc(w.N);
    w.n = fread(w.buffer, 1, w.N, w.fptr);
    w.i = 0;
    if (ferror(w.fptr)) { fclose(w.fptr); w.fptr = NULL; }
    w.width = (int32_t)vh->width;
    w.height = (int32_t)vh->height;
    return w;
}

// Make sure there is always plenty of buffer past the current index.
// Ideally we achieve this by reading the file, but whatever.
bool dbde_advance_file_buffer(dbde_file_walker &w) {
    int required = w.height * w.width;
    required += required/8 + 32;
    if (w.n - w.i < required) {
        if (w.i < w.n) memmove(w.buffer, w.buffer + w.i, w.n - w.i);
        w.n -= w.i;
        w.i = 0;
    }
    size_t n = feof(w.fptr) ? 0 : fread(w.buffer + w.n, 1, w.N - w.n, w.fptr);
    if (ferror(w.fptr)) return false;
    w.n += n;
    return true;
}

bool dbde_walk_a_file(dbde_file_walker *walker, frame_header *fh, uint8_t *image) {
    size_t npix = (size_t)(walker->width*walker->height);
    if (walker->n - walker->i < (npix + npix/8 + 32)) {
        if (!dbde_advance_file_buffer(*walker)) { fclose(walker->fptr); walker->fptr = NULL; return false; }
        if (walker->n - walker->i < 20) return false;
    }
    uint8_t *buf = walker->buffer + walker->i;
    *fh = dbde_unpack_frame(&buf, walker->width, walker->height, image);
    if (fh->u64s != 2) { fclose(walker->fptr); walker->fptr = NULL; return false; }
    if (buf - walker->buffer > walker->n) { printf("BUFFER OVERRUN!!!\n"); fclose(walker->fptr); walker->fptr = NULL; return false; }
    if (buf - walker->buffer < walker->i) { printf("BUFFER UNDERRUN!!\n"); fclose(walker->fptr); walker->fptr = NULL; return false; }
    walker->i = buf - walker->buffer;
    return true;
}

void dbde_end_file_walk(dbde_file_walker *walker) {
    if (walker->fptr) fclose(walker->fptr);
    walker->fptr = NULL;
}

