#pragma once

#include <png.h>

class Palette {
public:
  int width;
  int height;
  png_bytepp data;

  Palette(int width, int height) {
    this->width = width;
    this->height = height;
    this->data = new png_bytep [height];
    for(int i=0;i<height;i++) {
      this->data[i] = new png_byte [3*width];
      for(int j=0;j<width*3;j++) {
        this->data[i][j] = 0xFF;      // init with white color
      }
    }
  }
  ~Palette() {
    for(int i=0;i<height;i++) delete this->data[i];
    delete this->data;
  }

  // x,y starts with 0
  void setColor(int x, int y, png_byte r, png_byte g, png_byte b) {
    this->data[y][3*x] = r;
    this->data[y][3*x+1] = g;
    this->data[y][3*x+2] = b;
  }
};

void writeImage(const char *filename, Palette &pal) {
  int width = 100, height = 100;
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;
	
	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, pal.width, pal.height,
			8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	// Set title
	// if (title != NULL) {
	// 	png_text title_text;
	// 	title_text.compression = PNG_TEXT_COMPRESSION_NONE;
	// 	title_text.key = "Title";
	// 	title_text.text = title;
	// 	png_set_text(png_ptr, info_ptr, &title_text, 1);
	// }

	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep) malloc(3 * pal.width * sizeof(png_byte));

	// Write image data
  png_write_image(png_ptr, pal.data);
	// int x, y;
	// for (y=0 ; y<pal.height ; y++) {
	// 	for (x=0 ; x<pal.width ; x++) {
	// 		setRGB(&(row[x*3]), pal.data[y][x]);
	// 	}
	// 	png_write_row(png_ptr, row);
	// }

	// End write
	png_write_end(png_ptr, NULL);

	
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);

}
