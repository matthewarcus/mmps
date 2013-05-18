// $Revision: 1.2 $
// circularize.cpp
// (C) 2004 by Matthew Arcus

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "constants.h"
#include "image.h"

const char* usage = 
"Usage: %s [-bg r:g:b][-fuzz x%%][-f filename][-out filename]\n";

static inline double square(double x) {
  return x * x;
}

int main(int argc, char* argv[])
{
  Image inImage;
  const char* infilename = "images/earth.ppm";
  const char* outfilename = NULL;
  Rgb background = black;
  Rgb ibackground = black;
  int width;
  int height;
  double fuzz = 2.0;

  // array of pointers. Dynamic allocation is less typing than
  // stack allocating everything.
  OptSpec *optspecs [] = {
    // General options
    new RgbSpec ("-bg", background),
    new RgbSpec ("-ibg", ibackground),
    new DoubleSpec ("-fuzz", fuzz),
    new StringSpec ("-f", infilename),
    new StringSpec ("-out", outfilename),
    NULL // Finish up
  };
  
  char *progname = argv[0];
  argc--;argv++;
  while (argc > 0) {
    int result = OptSpec::ReadOpts(argc, argv, optspecs);
    if (result != 1) {
      error(usage, progname);
    }
    argc--; argv++;
  }

  inImage.Read(infilename);

  height = inImage.Height();
  // Make sure image is big enough to fit a full circular image
  width = inImage.Height();
  double radius = height/2.0;

  Image outImage(width, height);
  //fprintf(stderr, "Size %d\n", height);

  for (int row = 0; row < inImage.Height(); row++) {
    int first = -1;
    int last = -1;
    int column;
    // First scan the row in the input image to find the bounds.
    for (column = 0; column < inImage.Width(); column++) {
      Rgb rgb;
      inImage.GetRgb(column, row, rgb);
      if (!rgb.Match(ibackground, fuzz)) {
	first = column;
	break;
      }
    }
    if (first < 0) {
      // No non-background pixel found, what to do?
      first = last = inImage.Width()/2;
    } else {
      // Now look from the right
      for (column = inImage.Width()-1; column >= 0; column--) {
	Rgb rgb;
	inImage.GetRgb(column, row, rgb);
	if (!rgb.Match(ibackground, fuzz)) {
	  last = column+1;
	  break;
	}
      }
    }
    //fprintf(stderr, "%d %d\n", first, last);
    // Work out the half width of the image on this scan line
    double hw = sqrt(square(radius) - square(radius - (row + 0.5)));
    double a = (last-first)/(2.0*hw);
    double b = first - a * (radius - hw);
    for (column = 0; column < width; column++) {
      // Find the corresponding pixel in the input image
      int p = int(floor(a * column + b));
      if (p >= first && p < last) {
	Rgb rgb;
	inImage.GetRgb(p,row,rgb);
	outImage.SetRgb(column,row,rgb);
      } else {
	outImage.SetRgb(column,row,background);
      }
    }
  }
  outImage.Write(outfilename);
  return 0;
}
