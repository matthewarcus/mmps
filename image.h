// $Revision: 1.4 $
// image.h
// (C) 2004 by Matthew Arcus

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#if !defined IMAGE_H
#define IMAGE_H

#include <string.h>
#include "utils.h"

// Define this if its OK to cast an unaligned char pointer to RGB *
//#define UNALIGNED
typedef unsigned char Byte;
typedef Byte Color;

// rgb colors
class Rgb
{
 public:
  inline Rgb(Color r_ = 0, Color g_ = 0, Color b_ = 0) : r(r_), g(g_), b(b_) {};
  Rgb(const char *s);
  inline void Dim(double d) {
    r = Color(r*d);
    g = Color(g*d);
    b = Color(b*d);
  };
  inline void Dim() { r = r/2; g = g/2; b = b/2; };
  bool Match(const Rgb &rgb, double fuzz);
 public:
  Color r;
  Color g;
  Color b;
};

static const Rgb black(0,0,0);
static const Rgb white(255,255,255);
static const Rgb red(255,0,0);
static const Rgb green(0,255,0);
static const Rgb blue(0,0,255);
static const Rgb cyan(0,255,255);
static const Rgb magenta(255,0,255);
static const Rgb yellow(255,255,0);
static const Rgb orange(255,150,0);

class RgbSpec : public OptSpec
{
 public:
  RgbSpec (const char *optstring_, Rgb &value_)
    : OptSpec(optstring_), value(value_) {};
  int ReadOpt(int &argc, char** &argv);
 private:
  Rgb &value;
};

class LinePlotter
{
 public:
  virtual bool GetXY(double t, double&x, double&y) const = 0;
  virtual ~LinePlotter(){};
};

class MapFunction
{
 public:
  virtual double Scale(long width, long height) = 0;
  virtual void InitY(double y) = 0;
  virtual Rgb operator()(double x, double y) const = 0;
  virtual ~MapFunction(){};
};

// An image, typically the data is mmapped to a file
class Image {
 public:
  Image(long width = 0, long height = 0);
  ~Image();
  void Copy (const Image& image);
  void Read (const char *filename, bool writeable = false);
  void Create(long width, long height);
  void Write (const char *filename) const;
  void Close ();
  void Map (MapFunction& f, double xoff, double yoff);
  // Add a latitude, longitude grid to the image
  void AddGrid(double gridx, double gridy, const Rgb& color);
  void DrawLine(long x0, long y0, long x1, long y2, const Rgb& color);
  void PlotLine(double t0, double t1, const LinePlotter& lineplotter,
		const Rgb& color, int numinitialpoints = 29); // Why 29?

  void PlotPoint(double x, double y, int r, const Rgb &rgb);

  inline long Height() const { return height; };
  inline long Width() const { return width; };
  
  void Fill (const Rgb& value);
  inline void CheckSetRgb(long i, long j, const Rgb& rgb) {
    if (i >= 0 && i < width && j >= 0 && j < height) {
      SetRgb(j * width + i, rgb);
    }
  }
  inline void CheckGetRgb(long i, long j, Rgb& rgb) {
    if (i >= 0 && i < width && j >= 0 && j < height) {
      GetRgb(j * width + i, rgb);
    }
    else {
      rgb = Rgb(0,0,0);
    }
  }
  inline void SetRgb(long i, long j, const Rgb& rgb) {
    SetRgb(j * width + i, rgb);
  }
  inline void GetRgb(long i, long j, Rgb& rgb) const {
    GetRgb(j * width + i, rgb);
  }
  inline void SetRgb(long i, const Rgb& rgb) {
#if defined UNALIGNED
    Rgb* p = reinterpret_cast<Rgb *>(data+i+i+i);
    *p = rgb;
#else
    Byte *p = reinterpret_cast<Byte*>(data+i+i+i);
    *p++ = rgb.r;
    *p++ = rgb.g;
    *p = rgb.b;
#endif
  }
  inline void GetRgb(long i, Rgb& rgb) const {
#if defined UNALIGNED
    Rgb* p = reinterpret_cast<Rgb*>(data+i+i+i);
    rgb = *p;
#else
    Byte *p = reinterpret_cast<Byte*>(data+i+i+i);
    rgb.r = *p++;
    rgb.g = *p++;
    rgb.b = *p;
#endif
  }
  inline void Fill (Byte value) {
    memset(data, value, 3 * height * width);
  }
 private:
  void PlotLineAux(double t0, double t1, const LinePlotter& lineplotter,
		   const Rgb& color, int depth);

 private:
  Image(const Image&);
  Image &operator=(const Image&);
 private:
  long width;
  long height;
  char *data;
  void *mmapdata;
  unsigned long mmapsize;
  int mmapfd;
};
#endif
