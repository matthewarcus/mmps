// $Revision: 1.3 $
// image.cpp
// (C) 2004 by Matthew Arcus

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
// Need to include this eg. for cygwin to avoid CR/LF problems
// O_BINARY hopefully will be defined in fcntl.h if we need it.
#if defined O_BINARY
#include <io.h>
#endif
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "utils.h"
#include "image.h"

//#define ROUND floor
#define ROUND roundint

const unsigned int LINELENGTH = 128;

const long double pi = 3.14159265358979323846;
inline double radians(double deg) { return deg * pi / 180.0; };

double greenwichlat = radians(51.0 + (28.0 + 38.0/60.0)/60.0);

static inline double square(double x) {
  return x * x;
}

Rgb::Rgb(const char *s)
{
  int r_,g_,b_;
  if (sscanf(s,"%d:%d:%d",&r_,&g_,&b_) == 3){
    r = r_;
    g = g_;
    b = b_;
    // Should use a lookup table here really.
  } else {
    const char *p = s;
    const char *d = "dark";
    bool dark = false;
    if (memcmp(p,d,strlen(d)) == 0) {
      p += strlen(d);
      dark = true;
    }
    if (strcmp(p,"black") == 0) *this = black;
    else if (strcmp(p,"white") == 0) *this = white;
    else if (strcmp(p,"red") == 0) *this = red;
    else if (strcmp(p,"green") == 0) *this = green;
    else if (strcmp(p,"blue") == 0) *this = blue;
    else if (strcmp(p,"cyan") == 0) *this = cyan;
    else if (strcmp(p,"magenta") == 0) *this = magenta;
    else if (strcmp(p,"yellow") == 0) *this = yellow;
    else if (strcmp(p,"orange") == 0) *this = orange;
    else fprintf(stderr,"Invalid color specification: %s\n",s);
    if (dark) {
      double f = 0.6;
      r = (unsigned char)(r * f);
      g = (unsigned char)(g * f);
      b = (unsigned char)(b * f);
    }
  }
}

bool Rgb::Match(const Rgb &rgb, double fuzz)
{
  double total = (square(r - rgb.r) + square(g - rgb.g) + square(b - rgb.b)) / 
                 (3 * 255.0 * 255.0);
  return (total <= fuzz/100.0);
}
  
// Read an Rgb spec
int RgbSpec::ReadOpt(int &argc, char** &argv) {
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    argc--; argv++;
    if (argc == 0) {
      fprintf(stderr,"Need a parameter for %s\n", optstring);
      return -1;
    } else {
      value = Rgb(argv[0]);
      return 1;
    }
  }
}

Image::Image(long w, long h)
  : width(0), height(0), data(0), mmapdata(0), mmapsize(0)
{
  if (w != 0) {
    Create(w, h);
  }
}

Image::~Image() 
{
  Close();
}

void Image::Create(long w, long h) 
{
  width = w;
  height = h;
  data = new char[(long) width * height * 3];
}

void Image::Copy (const Image& image) {
  if (data == NULL) {
    width = image.width;
    height = image.height;
    data = new char[(long) width * height * 3];
    memcpy(data, image.data,(long) width * height * 3);
  } else {
    for (long y = 0; y < height && y < image.height; y++) {
      for (long x = 0; x < width && x < image.width; x++) {
	Rgb rgb;
	image.GetRgb(x,y,rgb);
	SetRgb(x,y,rgb);
      }
    }
  }
}

void Image::Close()
{
  if (mmapdata != NULL) {
    if (munmap(mmapdata, mmapsize)) {
      syserror("Error unmapping input file");
    }
    close(mmapfd);
    mmapdata = NULL;
    mmapsize = 0;
    mmapfd = 0;
    data = NULL;
  } else if (data != NULL) {
    delete (data);
    data = NULL;
  }
}

void Image::Fill (const Rgb& value)
{
  for (long i = 0; i < (long) width * height; i++) {
    SetRgb(i,value);
  }
}

// Read a line from filebuffer, ignoring blank & comment lines
static void GetLine (char* &dataptr, char* dataend, 
	      char* linebuffer, int length)
{
  do {
    int i = 0;
    if (dataptr == dataend) {
      error("Unexpected end of file\n");
    }
    while (dataptr < dataend && i <length-1) {
      char c = *dataptr;
      dataptr++;
      if (c == '\n') {
	break;
      }
      linebuffer[i] = c;
      i++;
    }
    linebuffer[i] = 0;
  } while (linebuffer[0] == '#' || linebuffer[0] == 0);
}

// Read an image from a file
// Expects format:
// P6
// <width> <height>
// <ncolors>
// <raw data in bytes triplets>
//
// In header, blank lines & lines starting with '#' are ignored
void Image::Read (const char* filename, bool writeable)
{
  int fd = -1;
  struct stat filestat;
  size_t filesize = 0;
  char *filebuffer = NULL;
  char *fileend = NULL;
  char *fileptr = NULL;
  char linebuffer[LINELENGTH];
  long w; long h; int ncolors; long datasize;

  fd = open(filename, writeable ? O_RDWR : O_RDONLY);
  if (fd < 0) {
    syserror("Couldn't open input file");
  }
  if (fstat(fd, &filestat) != 0) {
    syserror("Couldn't stat file");
  }
  filesize = filestat.st_size;
  // Now try to mmap the file
  int prot = PROT_READ;
  if (writeable) prot |= PROT_WRITE;
  filebuffer = static_cast<char *>(mmap(0, filesize, prot, MAP_SHARED, fd, 0));
  if (filebuffer == NULL) {
    syserror("Couldn't mmap file");
  }
  mmapdata = filebuffer;
  mmapsize = filesize;
  mmapfd = fd;

  if (strncmp(filebuffer,"P6",2) != 0) {
    error("Invalid file\n");
  }
  fileptr = filebuffer + 2;
  fileend = filebuffer + filesize;
  GetLine(fileptr, fileend, linebuffer, LINELENGTH);
  if (sscanf(linebuffer, "%ld %ld", &w, &h) != 2){
    error("Invalid file\n");
  }
  GetLine(fileptr, fileend, linebuffer, LINELENGTH);
  if (sscanf(linebuffer, "%d", &ncolors) != 1){
    error("Invalid file\n");
  }
  if (ncolors != 255){
    error("Unexpected number of colors\n");
  }
  datasize = (long) w*h*3;
  // Now we hope that fileptr is pointing at the rgb data
  if (fileend - fileptr != datasize) {
    error("Not enough data in file\n");
  }
  width = w;
  height = h;
  data = fileptr;
}

void Image::Write (const char* filename) const
{
  if (mmapdata != NULL) {
    fprintf(stderr, "mmapped image - skipping WriteImage\n");
  } else {
    int fd = 1; // stdout
    char header[128];
    size_t datasize;
    if (filename != NULL) {
      fd = creat(filename, 0666);
      if (fd < 0) {
	syserror("Can't open output file");
      }
    }
#if defined O_BINARY
    // Ghastly hack to ensure CR/LF substitution not done, eg. on cygwin
    // Set the output type to binary, not text.
    setmode(fd,O_BINARY);
#endif
    sprintf(header, "P6\n%ld %ld\n255\n", width, height);
    if (write(fd, header, strlen(header)) < 0) {
      syserror("write failed");
    }
    datasize = (long) width * height * 3;
    char *dataptr = data;
    while (datasize > 0) {
      size_t MAXDATA = 64UL*1024*1024;
      ssize_t written = write(fd, dataptr, std::min(MAXDATA,datasize));
      fprintf(stderr,"%zd bytes written\n", written);
      if (written < 0) {
	syserror("write failed");
      }
      datasize -= written;
      dataptr += written;
    }
    if (filename != NULL) {
      close(fd);
    }
  }
}

// Add a latitude, longitude grid to the image
void Image::AddGrid(double gridx, double gridy, const Rgb& color)
{
  double ystep = gridy * height/180.0;
  double xstep = gridx * width/360.0;
  // Quicker this way around with big images
  for (long y = 0; y < height; y++) {
    for (double xgrad = 0.0; ROUND(xgrad) < width; xgrad += xstep) {
      long x = long(ROUND(xgrad));
      SetRgb((long) y * width + x, color);
    }
  }
  for (double ygrad = 0.0; ROUND(ygrad) < height; ygrad += ystep) {
    long y = long(ROUND(ygrad));
    for (long x = 0; x < width; x++) {
      SetRgb((long) y * width + x, color);
    }
  }
#if 0
  long y = long(ROUND(height/2.0 - greenwichlat * height / pi));
  for (long x = 0; x < width; x++) {
    SetRgb((long) y * width + x, red);
  }
#endif
}

void Image::Map (MapFunction& f, double xoff, double yoff)
{
  double hh = height/2;
  double hw = width/2;
  double scale = f.Scale(width,height);
  //fprintf(stderr, "%f %f %f\n", xoff, yoff, scale);
  for (long y = 0; y < height; y++) {
    double y0 = (hh-y)*scale+yoff;
    f.InitY(y0);
    for (long x = 0; x < width; x++) {
      double x0 = (x-hw)*scale+xoff;
      //fprintf (stderr, "%d %d %f %f\n", x,y,x0,y0);
      SetRgb((long) y*width+x, f(x0,y0));
    }
  }
}

void Image::DrawLine(long x0, long y0, long x1, long y1, const Rgb& rgb)
{
  long nx = abs(x1 - x0) + 1;
  long ny = abs(y1 - y0) + 1;
  if (ny > nx) {
    int xinc = (x0 > x1) ? -1 : 1;
    int yinc = (y0 > y1) ? -1 : 1;
    long x = x0;
    long c = 0;
    for (long y = y0; y != y1+yinc; y+=yinc) {
      CheckSetRgb(x,y,rgb);
      if ((xinc == 1 && x > x1) || (xinc == -1 && x < x1) ) {
	fprintf(stderr, "Oops 1: %ld %ld %ld %ld %ld\n", x0, y0, x1, y1, x);
	error("Foo");
      }
      c += nx;
      if (c >= ny) {
	c -= ny;
	x += xinc;
      }
    }
  } else {
    int xinc = (x0 > x1) ? -1 : 1;
    int yinc = (y0 > y1) ? -1 : 1;
    long y = y0;
    long c = 0;
    for (long x = x0; x != x1+xinc; x+=xinc) {
      CheckSetRgb(x,y,rgb);
      if ((yinc == 1 && y > y1) || (yinc == -1 && y < y1) ) {
	fprintf(stderr, "Oops 2: %ld %ld %ld %ld %ld\n", x0, y0, x1, y1, y);
	error("Foo");
      }
      c += ny;
      if (c >= nx) {
	c -= nx;
	y += yinc;
      }
    }
  }
}

void Image::PlotPoint(double x, double y, int r, const Rgb &rgb)
{
  //fprintf(stderr,"%g %g %g\n",x,y,r);
  for (int j = -r; j <= r; j ++) {
    for (int i = -r; i <= r; i++) {
      CheckSetRgb(long(ROUND(x+i)),long(ROUND(y+j)),rgb);
    }
  }
}

// Plot a parametrically defined line betwen t0 and t1
// Break line up in to fixed number of chunks before calling
// the recursive version to get the detail.
void Image::PlotLine(double t0, double t1, 
		     const LinePlotter& lineplotter,
		     const Rgb& color,
		     int numinitialpoints)
{
  double inc = (t1-t0)/numinitialpoints;
  for (int i = 0; i < numinitialpoints; i++) {
    PlotLineAux(t0 + i * inc, t0 + (i+1) * inc,
		lineplotter, color, 10);
  }
}

void Image::PlotLineAux(double t0, double t1, 
			const LinePlotter& lineplotter,
			const Rgb& color,
			int depth)
{
  if (depth == 0) {
    return;
  }
  double x0 = 0.0, y0 = 0.0;
  double x1 = 0.0, y1 = 0.0;
  // The bools tell us if the inverse mapping is even defined
  bool defined0 = lineplotter.GetXY(t0, x0, y0);
  bool defined1 = lineplotter.GetXY(t1, x1, y1);
  // This really ought to take note of the curvature.
  if (defined0 && defined1 && fabs (x0-x1) + fabs(y0-y1) < 10) {
    long ix0 = long(ROUND(x0));
    long iy0 = long(ROUND(y0));
    long ix1 = long(ROUND(x1));
    long iy1 = long(ROUND(y1));
    DrawLine(ix0,iy0,ix1,iy1,color);
    //CheckSetRgb(ix0,iy0,Rgb(255,255,255));
    //CheckSetRgb(ix1,iy1,Rgb(255,255,255));
  } else {
    int margin = 50;
    // Consider whether to recurse
    // If both points are off to one side or above or below the image,
    // don't recurse, we assume that the line doesn't pass through.
    long sidex0 = (x0 > width + margin) ? 1 : (x0 < -margin) ? -1 : 0;
    long sidey0 = (y0 > height + margin) ? 1 : (y0 < -margin) ? -1 : 0;
    long sidex1 = (x1 > width + margin) ? 1 : (x1 < -margin) ? -1 : 0;
    long sidey1 = (y1 > height + margin) ? 1 : (y1 < -margin) ? 1 : 0;
    if ((defined0 || defined1) && sidex0 * sidex1 + sidey0 * sidey1 <= 0) {
#if 0
      double t2 = (t0 + t0 + t1) / 3;
      PlotLineAux(t0, t2, lineplotter, color, depth-1);
      double t3 = (t0 + t1 + t1) / 3;
      PlotLineAux(t2, t3, lineplotter, color, depth-1);
      PlotLineAux(t3, t1, lineplotter, color, depth-1);
#else
      double t2 = (t0 + t1) / 2;
      PlotLineAux(t0, t2, lineplotter, color, depth-1);
      PlotLineAux(t2, t1, lineplotter, color, depth-1);
#endif
    } else {
#if 0
      fprintf(stderr, "Giving up: %f %f %f %f %f %f\n",
              t0, t1, x0, y0, x1, y1);
#endif
    }
    
  }
}
