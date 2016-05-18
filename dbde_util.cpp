/** DBDE Util: converts Dynamic Bit Depth Encoding videos *
  * Copyright 2016 by Rex Kerr and Calico Life Sciences   *
  * This file distributed under the Apache License 2.0    **/

#include <smmintrin.h>
#include <stdint.h>
#include <string.h>
#include "dbde_util.h"

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
