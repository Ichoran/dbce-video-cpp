/** DBDE Util: converts Dynamic Bit Depth Encoding videos *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences   *
  * This file distributed under the Apache License 2.0    **/

#ifndef DBDE_UTIL
#define DBDE_UTIL

struct video_header {
    uint32_t u64s;
    uint64_t height;
    uint64_t width;
    uint64_t frame_hz;
};

struct frame_header {
    uint32_t u64s;
    uint64_t index;
    uint64_t reserved0;
};

struct dbde_data {
    uint32_t bits_len;
    uint8_t *bits;
    uint32_t mins_len;
    uint8_t *mins;
    uint32_t data_len;
    uint64_t *data;
};

video_header dbde_unpack_video_header(uint8_t *encoded);

int dbde_pack_video_header(video_header vh, uint8_t *encoded);

frame_header dbde_unpack_frame_header(uint8_t *encoded, int length);

int dbde_pack_frame_header(frame_header fh, uint8_t *encoded);

dbde_data dbde_unpack_data(uint8_t *encoded, int length);

int dbde_pack_data(dbde_data dd, uint8_t *encoded);

void dbde_decode_image(dbde_data dd, int W, int H, uint8_t *image);

uint8_t* dbde_decode_image_alloc(dbde_data dd, int W, int H);

void dbde_encode_image(uint8_t *image, int W, int H, dbde_data& dd);

dbde_data* dbde_encode_image_alloc(uint8_t *image, int W, int H);

#endif
