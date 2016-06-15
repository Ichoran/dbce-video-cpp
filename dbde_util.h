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

uint32_t dbde_pack_8x8(uint8_t *image, int stride, uint8_t *target);
uint32_t dbde_pack_8x8_partial(uint8_t *image, int stride, int rightmargin, int downmargin, uint8_t *target);

size_t dbde_pack_image(uint8_t *image, int W, int H, uint8_t *target);
size_t dbde_pack_frame_header(frame_header fh, uint8_t *target);
size_t dbde_pack_frame(uint64_t index, uint8_t *image, int W, int H, uint8_t *target);

size_t dbde_pack_video_header(video_header vh, uint8_t *target);

void dbde_unpack_8x8(uint8_t depth, uint8_t minval, uint8_t* packed, size_t stride, uint8_t *image);
void dbde_unpack_8x8_partial(uint8_t depth, uint8_t minval, uint8_t* packed, size_t stride, int rightmargin, int downmargin, uint8_t* image);

size_t dbde_unpack_image(uint8_t* packed, int W, int H, uint8_t *image);
frame_header dbde_unpack_frame(uint8_t **packed, int W, int H, uint8_t *image);

video_header dbde_unpack_video_header(uint8_t **packed);

struct dbde_file_walker {
    FILE *fptr;      // File we're reading from
    int32_t frames;  // Number of frames read
    size_t i;        // Index of unread data in buffer
    size_t n;        // Index of good data in buffer
    size_t N;        // Size of buffer in memory
    int32_t width;   // Width of images in file
    int32_t height;  // Height of images in file
    uint8_t *buffer; // The buffer
};

dbde_file_walker dbde_start_file_walk(const char* name, int frames_buffered, video_header *vh);
bool dbde_walk_a_file(dbde_file_walker *walker, frame_header* fh, uint8_t *image);
void dbde_end_file_walk(dbde_file_walker *walker);

#endif
