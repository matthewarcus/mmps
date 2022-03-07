// $Revision: 1.2 $
// combine.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "constants.h"
#include "image.h"

const char* usage = 
"Usage: %s file1 file1";

template <class T>
inline T max (T i, T j) {
  if (i < j) return j;
  else return i;
}

int main(int argc, char* argv[])
{
  Image image1;
  Image image2;
  const char* file1;
  const char* file2;
  const char* outfilename = NULL;
  Rgb background = black;
  int width;
  int height;
  double fuzz = 2.0;

  char *progname = argv[0];
  argc--;argv++;
  OptSpec *optspecs [] = {
    // General options
    new RgbSpec ("-bg", background),
    new DoubleSpec ("-fuzz", fuzz),
    new StringSpec ("-out", outfilename),
    NULL // Finish up
  };

  while (argc > 0) {
    int result = OptSpec::ReadOpts(argc, argv, optspecs);
    if (result != 1) {
      break;
    }
    argc--; argv++;
  }
  
  if (argc < 2) {
    error(usage, progname);
  }
  file1 = argv[0];
  file2 = argv[1];
  image1.Read(file1);
  image2.Read(file2);

  width = image1.Width();
  height = image1.Height();
  if (image2.Width() != width || image2.Height() != height) {
    error("Images must have same dimensions");
  }
  Image outImage(width, height);
  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {
      Rgb rgb1;
      Rgb rgb2;
      Rgb out;
      image1.GetRgb(col,row,rgb1);
      image2.GetRgb(col,row,rgb2);
      if (rgb1.Match(background,fuzz)) {
	out = rgb2;
      } else if (rgb2.Match(background, fuzz)) {
	out = rgb1;
      } else {
#if 0
	// Average the rgb values, this produces alarming discontinuities
	out.r = (rgb1.r + rgb2.r)/2;
	out.g = (rgb1.g + rgb2.g)/2;
	out.b = (rgb1.b + rgb2.b)/2;
#endif
	// But taking the maximum of the rgb values is smooth
	out.r = max(rgb1.r,rgb2.r);
	out.g = max(rgb1.g,rgb2.g);
	out.b = max(rgb1.b,rgb2.b);

      }
      outImage.SetRgb(col,row,out);
    }
 }
  outImage.Write(outfilename);
  return 0;
}
