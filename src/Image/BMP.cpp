#include "Image/BMP.hpp"
#include "Common/DebugAssert.hpp"
#include <fstream>

Bitmap::Bitmap(const char* fileName, dataSize_t width, dataSize_t height){
    constexpr dataSize_t headerSize = 14;
    constexpr dataSize_t infoHeaderSize = 40;

    this->fileName = fileName;
    this->pixels = new Pixel[width * height];
    this->width = width;
    this->height = height;
    dataSize_t uWidth = width * sizeof(Pixel);
    uWidth = ((uWidth + 3)/4) * 4;
    header.bfType = 0x4d42; //BM
    header.bfSize = headerSize + infoHeaderSize + uWidth * height + 2;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits = headerSize + infoHeaderSize;

    infoHeader.biSize = 40;
    infoHeader.width = width;
    infoHeader.height = height;
    infoHeader.biPlanes = 1;
    infoHeader.biBitCount = 8 * 3; //RGB
    infoHeader.biCompression = 0;
    infoHeader.biSizeImage =  uWidth * height;
    infoHeader.biXPelsPerMeter = 3779; //Approximately 96 DPI
    infoHeader.biYPelsPerMeter = 3779; //Approximately 96 DPI
    infoHeader.biClrUsed = 0;
    infoHeader.biClrImportant = 0;
}

Bitmap::~Bitmap(){
    Flush();
    delete[] pixels;
}

Pixel& Bitmap::operator[](dataSize_t idx){
    return pixels[idx];
}

void Bitmap::Flush(){
    std::ofstream file;
    file.open(fileName);
    DebugAssert(__FILE__, __LINE__,"File can not be opend", file.is_open());
    
    file.write((byte_t*)&header.bfType, sizeof(word_t));
    file.write((byte_t*)&header.bfSize, sizeof(dword_t));
    file.write((byte_t*)&header.bfReserved1, sizeof(word_t));
    file.write((byte_t*)&header.bfReserved2, sizeof(word_t));
    file.write((byte_t*)&header.bfOffBits, sizeof(dword_t));

    file.write((byte_t*)&infoHeader.biSize, sizeof(dword_t));
    file.write((byte_t*)&infoHeader.width, sizeof(int32_t));
    file.write((byte_t*)&infoHeader.height, sizeof(int32_t));
    file.write((byte_t*)&infoHeader.biPlanes, sizeof(word_t));
    file.write((byte_t*)&infoHeader.biBitCount, sizeof(word_t));
    file.write((byte_t*)&infoHeader.biCompression, sizeof(dword_t));
    file.write((byte_t*)&infoHeader.biSizeImage, sizeof(dword_t));
    file.write((byte_t*)&infoHeader.biXPelsPerMeter, sizeof(int32_t));
    file.write((byte_t*)&infoHeader.biYPelsPerMeter, sizeof(int32_t));
    file.write((byte_t*)&infoHeader.biClrUsed, sizeof(dword_t));
    file.write((byte_t*)&infoHeader.biClrImportant, sizeof(dword_t));
    dataSize_t padding = width * sizeof(byte_t) * 3;
    padding = ((padding + 3)/4) * 4 - width * sizeof(byte_t) * 3;
    for(int i = height - 1; i >= 0; --i){
        for(int j = 0; j < width; ++j){
            file.write(reinterpret_cast<byte_t*>(&pixels[i * width + j].b), sizeof(byte_t));
            file.write(reinterpret_cast<byte_t*>(&pixels[i * width + j].g), sizeof(byte_t));
            file.write(reinterpret_cast<byte_t*>(&pixels[i * width + j].r), sizeof(byte_t));
        }
        for(int j = 0; j < padding; ++j)
            file.write(0, padding);
    }
    file.close();
}