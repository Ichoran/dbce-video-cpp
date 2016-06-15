/** DBDE Util Test: tests of proper DBDE Util function  *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences *
  * This file distributed under the Apache License 2.0  **/

#include <smmintrin.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "dbde_util.h"

void dbde_print_ascii(int W, int H, int ox, int oy, int nx, int ny, uint8_t *image) {
    int *im = (int*)malloc(32*32*sizeof(int));
    int y = 0;
    for (int j = 0; j < ny; j++) {
        while ((j+1)*32 > (y+1)*ny) y++;
        int x = 0;
        for (int i = 0; i < nx; i++) {
            while ((i+1)*32 > (x+1)*nx) x++;
            im[y*32 + x] += image[(j+oy)*W + (i+ox)];
        }
    }
    int I = 0;
    int Z = im[0]; 
    for (int j = 0; j < 32; j++) for (int i = 0; i < 32 ; i++) {
        if (I < im[j*32+i]) I = im[j*32+i];
        if (Z > im[j*32+i]) Z = im[j*32+i];
    }
    I = (I-Z)/10;
    if (I < 1) I = 1;
    for (int j = 0; j < 32 ; j++) {
        for (int i = 0; i < 32 ; i++) {
            switch((im[32*j+i]-Z)/I) {
                case 0: printf("  "); break;
                case 1: printf("__"); break;
                case 2: printf(".."); break;
                case 3: printf(",,"); break;
                case 4: printf("++"); break;
                case 5: printf("oo"); break;
                case 6: printf("OO"); break;
                case 7: printf("XX"); break;
                case 8: printf("&&"); break;
                case 9: printf("##"); break;
                default: printf("@@");
            }
        }
        printf("\n");
    }
}

void dbde_dump_pgm(int W, int H, uint8_t *image, const char* filename) {
    FILE *f = fopen(filename, "wb");
    fprintf(f, "P2\n");
    fprintf(f, "%d %d\n",W, H);
    fprintf(f, "%d\n", 255);
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            if (i==0) fprintf(f, "%d", image[j*W+i]);
            else fprintf(f, " %d", image[j*W+i]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

bool dbde_util_unit_test() {
#define DUUT_FAIL printf("%d %d %d %d %d %d\n", nX, nY, bI, mI, w, h); free(src); free(pack); free(dst); return false
    int nX = 8 + (rand()%(16384-8));
    int nY = 8 + (rand()%(16384-8));
    nX = 8;
    nY = 8;
    if (nY*nX > 32*1024*1024) nY = (32*1024*1024)/nX;
    int bI = rand()%9;
    int mI = (bI < 8) ? rand()%(256 - (1 << bI)) : 0;
    int w = (nX+7)/8;
    int h = (nY+7)/8;
    uint8_t *src = (uint8_t*)malloc(nX * nY);
    uint8_t *pack = (uint8_t*)malloc(w*h*64 + w*h*2 + sizeof(frame_header) + 12 + 1024);
    uint8_t *dst = (uint8_t*)malloc(nX * nY);
    for (int i = 0; i < nX*nY; i++) src[i] = (uint8_t)(mI + (rand() % (1 << bI)));
    uint64_t ix = rand() | (((uint64_t) rand()) << 32);
    int si = dbde_pack_frame(ix, src, nX, nY, pack);
    uint8_t *packv = pack;
    frame_header fh = dbde_unpack_frame(&packv, nX, nY, dst);
    if (fh.u64s != 2) { printf("Error.\n"); DUUT_FAIL; }
    if (fh.index != ix) { printf("Bad index %llx (needed %llx)\n", (long long)fh.index, (long long)ix); DUUT_FAIL; }
    int wrong = 0;
    for (int i = 0; i < nX * nY; i++) if (src[i] != dst[i]) wrong += 1;
    if (wrong > 0) { printf("%d of %d pixels mismatched, %d bits %d min\n", wrong, nX*nY, bI, mI); DUUT_FAIL; }
    if (packv - pack > w*h*64 + w*h*2 + sizeof(frame_header) + 12) {
        printf("Overshot by %ld bytes\n", (packv - pack) - (w*h*64 + w*h*2 + sizeof(frame_header) + 12));
        DUUT_FAIL;
    }
    return true;
#undef DUUT_FAIL
}

bool mycmp_vh(video_header u, video_header v) {
    if (u.u64s != v.u64s) { printf("u64s: %d != %d", u.u64s, v.u64s); return false; }
    if (u.width != v.width) { printf("width: %lld != %lld", (long long)u.width, (long long)v.width); return false; }
    if (u.height != v.height) { printf("height: %lld != %lld", (long long)u.height, (long long)v.height); return false; }
    if (u.frame_hz != v.frame_hz) { printf("frameHz: %lld != %lld", (long long)u.frame_hz, (long long)v.frame_hz); return false; }
    return true;
}

bool mycmp_fh(frame_header f, frame_header g) {
    if (f.u64s != g.u64s) { printf("u64s: %d != %d", f.u64s, g.u64s); return false; }
    if (f.index != g.index) { printf("index: %lld != %lld", (long long)f.index, (long long)g.index); return false; }
    if (f.reserved0 != g.reserved0) { printf("reserved0: %lld != %lld", (long long)f.reserved0, (long long)g.reserved0); return false; }
    return true;
}

bool mycmp(const void* u, const void* v, size_t n) {
    bool same = true;
    for (int i = 0; i < n; i++) {
        if (*(((char*)u)+i) != *(((char*)v)+i)) {
            int j = i & 0xFFFFFFF0;
            printf("%06d: ",j);
            for (int k = j; k < j+16; k++) printf("%4d", *(((uint8_t*)u)+k));
            printf("\n");
            printf("%06d: ",j);
            for (int k = j; k < j+16; k++) printf("%4d", *(((uint8_t*)v)+k));
            printf("\n");
            printf("        ");
            for (int k = j; k < i; k++) printf("    ");
            printf(" ^^^\n");
            i = j+15;
            same = false;
        }
    }
    return same;
}

bool dbde_util_unit_minimal_8x16() {
    uint8_t image[128] = {
        0, 1, 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        2, 3, 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
        4, 5, 6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
        6, 7, 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
        5, 6, 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21,
        3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 18, 20,
        1, 2, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 15, 17, 19
    };
    uint8_t packed[128] = {
        3, 0, 0, 0,                // Video header
        8, 0, 0, 0, 0, 0, 0, 0,    // Height
        16, 0, 0, 0, 0, 0, 0, 0,   // Width
        1, 0, 0, 0, 0, 0, 0, 0,    // Frame rate
        2, 0, 0, 0,                // Frame header
        1, 0, 0, 0, 0, 0, 0, 0,    // Frame number
        0, 0, 0, 0, 0, 0, 0, 0,    // Unused
        2, 0, 0, 0,                // Number of uint64s of data
        4, 4,                      // Bits per block
        2, 0, 0, 0,                // Number of minimum values
        0, 8,                      // Minimum values
        8, 0, 0, 0,                // Number of packed U64s
        0x10, 0x32, 0x54, 0x76,    // Row 1 block 1
        0x32, 0x54, 0x76, 0x98,    // Row 2 block 1
        0x54, 0x76, 0x98, 0xBA,    // Row 3 block 1
        0x76, 0x98, 0xBA, 0xDC,    // Row 4 block 1
        0x87, 0xA9, 0xCB, 0xED,    // Row 5 block 1
        0x65, 0x87, 0xA9, 0xCB,    // Row 6 block 1
        0x43, 0x65, 0x87, 0xA9,    // Row 7 block 1
        0x21, 0x43, 0x65, 0x87,    // Row 8 block 1
        0x10, 0x32, 0x54, 0x76,    // Row 1 block 2
        0x32, 0x54, 0x76, 0x98,    // Row 2 block 2
        0x54, 0x76, 0x98, 0xBA,    // Row 3 block 2
        0x76, 0x98, 0xBA, 0xDC,    // Row 4 block 2
        0x87, 0xA9, 0xCB, 0xED,    // Row 5 block 2
        0x65, 0x87, 0xA9, 0xDB,    // Row 6 block 2
        0x43, 0x65, 0x87, 0xCA,    // Row 7 block 2
        0x21, 0x43, 0x75, 0xB9,    // Row 8 block 2
    };
    video_header vh = (video_header){3, 8, 16, 1};
    frame_header fh = (frame_header){2, 1, 0};
    uint8_t new_im[128];
    uint8_t z0[1024];
    uint8_t new_pk[128];
    uint8_t z1[1024];
    video_header new_vh;
    frame_header new_fh;
    uint8_t *b = packed;
    new_vh = dbde_unpack_video_header(&b);
    if (b - packed != 28) { printf("Read %d bytes instead of 28 to unpack video header\n", (int)(b - packed)); exit(1); }
    if (!mycmp_vh(vh, new_vh)) { printf("Unpacked video header wrong\n"); exit(1); }
    new_fh = dbde_unpack_frame_header(&b);
    if (b - packed != 48) { printf("Read %d bytes instead of 20 to unpack frame header\n", (int)(b - packed) - 28); exit(1); }
    if (!mycmp_fh(fh, new_fh)) { printf("Unpacked frame header wrong\n"); exit(1); }
    b = packed + 28;  // Re-read frame header so we can get whole frame at once
    new_fh = dbde_unpack_frame(&b, (int)vh.width, (int)vh.height, new_im);
    if (!mycmp_fh(fh, new_fh)) { printf("Unpacked frame header wrong on second pass\n"); exit(1); }
    if (b - packed != 128) { printf("Read %d bytes instead of 128 during unpacking\n", (int)(b-packed)); exit(1); }
    if (!mycmp(image, new_im, 128)) { printf("Image data does not match\n"); exit(1); }
    b = new_pk;
    int n;
    n = dbde_pack_video_header(vh, b); b += n;
    n = dbde_pack_frame(fh.index, image, (int)vh.width, (int)vh.height, b);
    if (n != 100) { printf("Wrong number of bytes written: %d\n", n); exit(1); }
    if (!mycmp(packed, new_pk, 128)) { printf("Packed data does not match\n"); exit(1); }
}

// Include down here to avoid letting the guys above know about operations we may have available on a dev machine
#include <x86intrin.h>

int main(int argc, char **argv) {
    uint8_t example[] = 
        {   25, 27, 23, 29, 22, 24, 29, 23, 25, 24,
            22, 24, 21, 25, 22, 27, 28, 21, 27, 26,
            25, 26, 22, 29, 25, 20, 28, 23, 26, 25,
            19, 23, 25, 41, 28, 19, 22, 25, 25, 27,
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

    dbde_util_unit_minimal_8x16();

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
    for (int bi = 0; bi < XRES*YRES + (XRES+YRES); bi++) big[bi] = rand() & 0xFF;
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

    for (int i = 0; i < 1024; i++) if (!dbde_util_unit_test()) { printf("Failed iteration %d\n", i); return 1; }

#ifdef DBDE_READ_FILE_TEST
    printf("Reading %s\n",DBDE_READ_FILE_TEST);
    int64_t ta6 = __rdtsc();
    video_header vh;
    frame_header fh;
    uint8_t *image = (uint8_t*)malloc(4096*4096);
    dbde_file_walker dfw = dbde_start_file_walk(DBDE_READ_FILE_TEST, 16, &vh);
    if (!dfw.fptr) printf("FAILED OPEN\n");
    else {
        int i = 0;
        while (dbde_walk_a_file(&dfw, &fh, image)) {
            i++;
            if (!(i%100)) {
                if (i == 100) dbde_dump_pgm(vh.width, vh.height, image, "frame_100.pgm");
                dbde_print_ascii(
                    vh.width, vh.height,
                    (vh.width > 1024) ? 768 : 0, (vh.height > 1024) ? 768 : 0,
                    (vh.width > 1024) ? 128 : vh.width, (vh.height > 1024) ? 128 : vh.height,
                    image
                );
                printf("%d\n",i);
            }
        }
        dbde_end_file_walk(&dfw);
    }
    int64_t tz6 = __rdtsc();
    printf("Elapsed: %e\n", (double)(tz6-ta6));
#endif
}
