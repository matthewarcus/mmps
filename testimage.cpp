// $Revision: 1.2 $
// testimage.cpp
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
#include <math.h>

#include "utils.h"
#include "image.h"

const long double pi = 3.14159265358979323846;

Rgb foreground(0,0,0);
Rgb background(255, 255, 255);

const char* usage = "Usage: %s [-w n] [-h n]\n";

class SpiralPlotter : public LinePlotter
{
  bool GetXY(double t, double&x, double&y) const
  {
    x = t * cos (1/t);
    y = t * sin (1/t);

    x *= 200;
    y *= 200;

    x += 200;
    y += 200;
    return true;
  }
};

class ThingPlotter : public LinePlotter
{
public:
  ThingPlotter(double a_, double b_) : a(a_), b(b_) {};
  bool GetXY(double t, double&x, double&y) const
  {
    x = cos (a * t);
    y = sin (b * t);

    x *= 200;
    y *= 200;

    x += 200;
    y += 200;
    return true;
  }
private:
  double a;
  double b;
};

class CirclePlotter : public LinePlotter
{
  bool GetXY(double t, double&x, double&y) const
  {
    x = cos(t);
    y = sin(t);

    x *= 200;
    y *= 200;

    x += 250;
    y += 200;
    return true;
  }
};

class TanPlotter : public LinePlotter
{
  bool GetXY(double t, double&x, double&y) const
  {
    x = t/(pi*2);
    y = tan(t)/(pi*2);

    x *= 200;
    y *= 200;

    x += 200;
    y += 200;
    return true;
  }
};

int main(int argc, char* argv[])
{
  int width = 360;
  int height = 180;
  int equalarea = false;

  char *progname = argv[0];
  argc--;argv++;
  while (argc > 0) {
    if (strcmp(argv[0],"-w") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      width = atoi(argv[0]);
    } else if (strcmp(argv[0], "-h") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      height = atoi(argv[0]);
    } else if (strcmp(argv[0], "-bg") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      int r,g,b;
      sscanf(argv[0],"%d:%d:%d",&r,&g,&b);
      background.r = r;
      background.g = g;
      background.b = b;
    } else if (strcmp(argv[0], "-fg") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      int r,g,b;
      sscanf(argv[0],"%d:%d:%d",&r,&g,&b);
      foreground.r = r;
      foreground.g = g;
      foreground.b = b;
    } else if (strcmp(argv[0], "-equalarea") == 0) {
      equalarea = true;
    } else {
      error(usage, progname);
    }
    argc--; argv++;
  }

  Image outImage;
  outImage.Create(width, height);
  outImage.Fill(0); // Set to all black

#if 0
  for (int i = 0; i < 62; i+=1) {
    double x0 = cos(i * 2 * pi/64);
    double y0 = sin(i * 2 * pi/64);
    double x1 = cos((i + 1) * 2 * pi/64);
    double y1 = sin((i + 1) * 2 * pi/64);
    x0 *= 200;
    y0 *= 200;
    x1 *= 200;
    y1 *= 200;

    if (i%2 == 1) {
      x0 += 200;
      y0 += 200;
      x1 += 200;
      y1 += 200;
    } else {
      x0 += 196;
      y0 += 196;
      x1 += 196;
      y1 += 196;
    }

    outImage.DrawLine(x0,y0,x1,y1,Rgb(255,0,0));
  }
#endif
#if 0
  SpiralPlotter plotter;
  CirclePlotter plotter2;
  TanPlotter plotter3;
  ThingPlotter plotter4(10,11);
  ThingPlotter plotter5(9,17);
  outImage.PlotLine(0, pi/2, plotter, Rgb(0,0,255));
  outImage.PlotLine(0, 2 * pi, plotter2, Rgb(0,255,0));
  outImage.PlotLine(-2*pi, 2*pi, plotter3, Rgb(255,0,0));
  outImage.PlotLine(0, 2 * pi, plotter4, Rgb(0,255,255), 197);
  outImage.PlotLine(0, 2 * pi, plotter5, Rgb(255,0,255), 197);
#endif
#if 0
  for (int i = 0; i < 20; i++) {
    outImage.DrawLine(0, 0, 200, i * 10, Rgb(255,0,0));
    outImage.DrawLine(i * 10, 0, 0, 200, Rgb(0,255,0));
  }
#endif

#if 1
  if (equalarea) {
    double xsize = width / 36.0;
    int ysize = height / 18;
    int index, x, y;
    for (index = 0, y = 0; y < height; y++) {
      double k = cos (pi * (0.5 * height - y)/height);
      for (x = 0; x < width; x++, index++) {
	double realx = k * (0.5 * width - x);
	bool neg = (realx < 0.0);
	realx = fabs(realx);
	int i = int(floor(realx/xsize))%2;
	if (neg) i = 1 - i;
	if (i == (y/ysize)%2){
	  outImage.SetRgb(index,foreground);
	} else {
	  outImage.SetRgb(index,background);
	}
    }
    }
  } else {
    int xsize = width / 36;
    int ysize = height / 18;
    int index = 0;
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++, index++) {
	if ((x/xsize)%2 == (y/ysize)%2) {
	  outImage.SetRgb(index,foreground);
	} else {
	  outImage.SetRgb(index,background);
	}
      }
    }
  }
#endif
  outImage.Write(NULL);
  return 0;
}
