// bitmap file drawer library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_BMP_HPP
#define INCLUDE_BMP_HPP

#include "Common/Type.hpp"

using byte_t = char;
using word_t = uint16_t;
using dword_t = uint32_t;

typedef struct{
    word_t bfType;
    dword_t bfSize;
    word_t bfReserved1;
    word_t bfReserved2;
    dword_t bfOffBits;
}BitMapHeader;

typedef struct{
    dword_t biSize;
    int32_t width;
    int32_t height;
    word_t biPlanes;
    word_t biBitCount;
    dword_t biCompression;
    dword_t biSizeImage;
    int32_t biXPelsPerMeter;
    int32_t biYPelsPerMeter;
    dword_t biClrUsed;
    dword_t biClrImportant;
}BitMapInfoHeader;

typedef struct{
    byte_t r;
    byte_t g;
    byte_t b;
    byte_t reserved;
}Pixel;

class Bitmap{
public:
    Bitmap(const char* fileName, dataSize_t width, dataSize_t height);
    ~Bitmap();
    Pixel& operator[](dataSize_t idx);
    void Flush();
private:
    Pixel* pixels;
    BitMapHeader header;
    BitMapInfoHeader infoHeader;
    dataSize_t width;
    dataSize_t height;
    const char* fileName;
};

#endif