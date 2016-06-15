# dbde-video-cpp
A library for Dynamic Bit Depth Encoded videos in C++ (useful for certain kinds of scientific imaging).

The format was created at HHMI Janelia by Lakshmi Ramasamy as a quick FPGA-friendly compression scheme.

This is a C++ library to convert frames encoded in DBDE format to and from arrays of pixels, using only the operations found on a modern CPU (at least SSE4.1).

## The format

### Entire file

A DBDE video is encoded as a single video header followed by zero or more blocks consisting of a frame header followed by data:

|DBDE file|
|---------|
|video header
|frame 1 header
|frame 1 data
| ...
|frame N header
|frame N data

The header does _not_ specify how many frames there are.  You know that you are done reading the file when you run out of bytes.

### Video header

The video header is 28 bytes long, and is packed as follows.  All values are stored in little-endian format.

| Data type |  Size (bytes) |  Meaning |
|-----------|---------|----------|
| I32       | 4       | Number of additional U64s in video header (always 3 for now)
| U64       | 8       | Image height
| U64       | 8       | Image width
| U64       | 8       | Frame rate in Hz (not needed for parsing)

### Frame header

The frame header is 20 bytes and contains no information essential for recovery of the video stream (aside from declaring how much space the header takes, in case a larger header is useful in the future)

| Data type |  Size (bytes) |  Meaning |
|-----------|---------|----------|
| I32       | 4       | Number of additional U64s in frame header (always 2 for now)
| U64       | 8       | Frame number (frames are in order but some may have been dropped)
| U64       | 8       | Reserved

### Frame data

#### How images are packed

Images are tiled in 8x8 blocks in row major order and are "constant padded" in that if number of rows or columns is not divisible by 8, the part of the 8x8 block that is missing is filled in left-to-right with the last valid value to get full rows, and then any missing rows are filled in with the last full row.

Within each tile, the minimum and maximum value is computed.  The minimum number of bits needed to store that range is used, and the bits are packed into U64s (least significant first).  Note that because there are integer numbers of bits and 64 pixels, there will always be an integer number of U64s to store the entire data.

#### Data layout

Let us use `W` to denote the width in pixels of the entire image, and `H` to denote the height.  Then there are `w = ceil(W/8)` blocks across and `h = ceil(H/8)` blocks down, for a total of `h*w` blocks.

| Data type |  Size (bytes) |  Meaning |
|-----------|---------|----------|
| I32       | 4       | Number of bytes in bit depth arrays (always `h*w`)
|  U8       | `h*w`   | Bit depth in each block needed to cover min-max range (row-major order)
| I32       | 4       | Number of bytes in per-block minimum intensity
| U8        | `h*w`   | Minimum intensities (note: could expand size to handle higher bit depth images!)
| I32       | 4       | Number of U64s in data block (equal to the sum of the number of bits)
| U64       | varies  | Bits for pixels in each block (row-major order within blocks, then sequentially across blocks in row-major order)

#### An example

Consider the following 10x10 8-bit image data (rows and columns are numbered):

|     |0  |  1|  2|  3|  4|  5|  6|  7|  8|  9|
|-----|---|---|---|---|---|---|---|---|---|---|
|**0**|25|27|23|29|22|24|29|23|25|24
|**1**|22|24|21|25|22|27|28|21|27|26
|**2**|25|26|22|29|25|20|28|23|26|25
|**3**|19|23|25|21|28|19|22|25|25|27
|**4**|27|25|30|28|25|23|27|26|24|24
|**5**|31|30|31|28|29|26|24|25|27|26
|**6**|30|28|32|25|28|27|28|27|26|26
|**7**|29|31|31|32|29|29|25|22|24|25
|**8**|31|34|33|31|30|29|28|28|26|26
|**9**|34|34|35|35|33|28|29|28|26|26

First, we consider the first 8x8 block and find the maximum and minimum value therein:

|     |0  |  1|  2|  3|  4|  5|  6|  7|  8|  9|
|-----|---|---|---|---|---|---|---|---|---|---|
|**0**|`25`|`27`|`23`|`29`|`22`|`24`|`29`|`23`|25|24
|**1**|`22`|`24`|`21`|`25`|`22`|`27`|`28`|`21`|27|26
|**2**|`25`|`26`|`22`|`29`|`25`|`20`|`28`|`23`|26|25
|**3**|**_19_**|`23`|`25`|`21`|`28`|`19`|`22`|`25`|25|27
|**4**|`27`|`25`|`30`|`28`|`25`|`23`|`27`|`26`|24|24
|**5**|`31`|`30`|`31`|`28`|`29`|`26`|`24`|`25`|27|26
|**6**|`30`|`28`|**_32_**|`25`|`28`|`27`|`28`|`27`|26|26
|**7**|`29`|`31`|`31`|`32`|`29`|`29`|`25`|`22`|24|25
|**8**|31|34|33|31|30|29|28|28|26|26
|**9**|34|34|35|35|33|28|29|28|26|26

We find that the maximum and minimum pixels cover a range of 14, so 4 bits are enough to represent the data.  We remember these numbers: **4** bits, **19** is the minimum value, and subtract off the minimum from that part of the image and convert to bits (here written as hexidecimal, 0-F):

|     |0 | 1| 2| 3| 4| 5| 6| 7
|-----|---|---|---|---|---|---|---|---
|**0**|6|8|4|A|3|5|A|4
|**1**|3|5|2|6|3|8|9|2
|**2**|6|7|3|A|6|1|9|4
|**3**|0|4|6|2|9|0|3|6
|**4**|8|6|B|9|6|4|8|7
|**5**|C|B|C|9|A|7|5|6
|**6**|B|9|D|6|9|8|9|8
|**7**|A|C|C|D|A|A|6|3

We pack these into U64s in least-significant order, which means (as we write hexidecimal in most-significant order) we have the four U64s 0x298362534A53A486, 0x630926404916A376, 0x657A9CBC78469B68, and 0x36AADCCA89896D9B for our data.

We then move right and consider the next (partial) block:

|     |0  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|-----|---|---|---|---|---|---|---|---|---|---|
|**0**|25|27|23|29|22|24|29|23|`25`|**_24_**
|**1**|22|24|21|25|22|27|28|21|**_27_**|`26`
|**2**|25|26|22|29|25|20|28|23|`26`|`25`
|**3**|19|23|25|21|28|19|22|25|`25`|`27`
|**4**|27|25|30|28|25|23|27|26|`24`|`24`
|**5**|31|30|31|28|29|26|24|25|`27`|`26`
|**6**|30|28|32|25|28|27|28|27|`26`|`26`
|**7**|29|31|31|32|29|29|25|22|`24`|`25`
|**8**|31|34|33|31|30|29|28|28|26|26
|**9**|34|34|35|35|33|28|29|28|26|26

We find that we need **2** bits per pixel and our minimum value is **24**.  We copy the last value rightwards, subtract off the minimum, and pad with zeros:

|     | 8 | 9 |
|-----|---|---
|**0**|1|0|0|0|0|0|0|0|0|0
|**1**|3|2|2|2|2|2|2|2|2|2
|**2**|2|1|1|1|1|1|1|1|1|1
|**3**|1|3|3|3|3|3|3|3|3|3
|**4**|0|0|0|0|0|0|0|0|0|0
|**5**|3|2|2|2|2|2|2|2|2|2
|**6**|2|2|2|2|2|2|2|2|2|2
|**7**|0|1|1|1|1|1|1|1|1|1

We then pack these into two U64s like so

0xFFFD5556AAAB0001, 0x5554AAAAAAAB0000.

We then move down to the next row:

|     | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|-----|---|---|---|---|---|---|---|---|---|---|
|**0**|25|27|23|29|22|24|29|23|25|24
|**1**|22|24|21|25|22|27|28|21|27|26
|**2**|25|26|22|29|25|20|28|23|26|25
|**3**|19|23|25|21|28|19|22|25|25|27
|**4**|27|25|30|28|25|23|27|26|24|24
|**5**|31|30|31|28|29|26|24|25|27|26
|**6**|30|28|32|25|28|27|28|27|26|26
|**7**|29|31|31|32|29|29|25|22|24|25
|**8**|`31`|`34`|`33`|`31`|`30`|`29`|**_28_**|`28`|26|26
|**9**|`34`|`34`|**_35_**|`35`|`33`|`28`|`29`|`28`|26|26

We need **3** bits and our minimum value is **28**.  Our data is:

|     | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|-----|---|---|---|---|---|---|---|---|
|**8**|3|6|5|3|2|1|0|0
|**9**|6|6|7|6|5|0|1|0

followed by a bunch of repeats of row **9**.  We pack these into three U64s, which turn out to be 0x5DF6045DF600A773, 0xF6045DF6045DF604, 0x045DF6045DF6045D.

Finally, the last block at the bottom-right corner is all 26, so we need **0** bits and our minimum value is **26**.

Thus, we encode the image data as:

```
0x00000004            // Number of blocks
0x04 0x02 0x03 0x00   // Number of bits per block
0x00000004            // Number of min values
0x13 0x18 0x1C 0x1A   // Minimum values (1/byte)
0x298362534A53A486    // Image data
0x630926404916A376
0x657A9CBC78469B68
0x36AADCCA89896D9B
0xFFFD5556AAAB0001
0x5554AAAAAAAB0000
0x5DF6045DF600A773
0xF6045DF6045DF604
0x045DF6045DF6045D
```
