// $Revision: 1.3 $
// stars.cpp
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
#include "constants.h"
#include "transform.h"
#include "version.h"

#define LINESIZE 256

const double asecperyear = 2 * pi / 1296000.0;

double topmag = 1.0;
double magfact = 0.666;
bool dobv = false;
double bvfact = 4.0;
double year = 0.0;

// Map a magnitude to a linear level <= 1.0
double convertmag (double mag) {
  mag -= topmag;
  if (mag < 0) mag = 0;
  return pow(magfact, mag);
}

void AddStars(const char* starlist, Image& image,
	      Transform& transform, const TransformParams& params)
{
  int width = image.Width();
  int height = image.Height();

  FILE* starfile;
  if (starlist == NULL) {
    starfile = stdin;
  } else {
    starfile = fopen(starlist, "r");
    if (starfile == NULL) {
      syserror("Can't open starlist");
    }
  }
  char linebuffer[LINESIZE];
  int i = 0;
  while (fgets(linebuffer, LINESIZE, starfile) != NULL) {
    // Read the appropriate entry in the star list
    double mag, bv, ra, dec, raproper, decproper;
    if (sscanf (linebuffer, "%lf %lf %lf %lf %lf %lf", 
		&mag, &bv, &ra, &dec,
		&raproper, &decproper) != 6) {
      syserror("sscanf failure -bad input?");
    }

    // Change right ascension and declination to sensible units (radians)
    ra = pi * ra / 12.0;
    dec = pi * dec / 180.0;

    // And put RA into to the expected range
    if (ra > pi) { ra -= 2 * pi; }

    // Apply proper motion, if required
    if (year != 0.0) {
      // Allegedly, both raproper and dec proper are in arcseconds per year
      ra += year * raproper * asecperyear;
      dec += year * decproper * asecperyear;
      // And put back into range
      if (dec > piovertwo) { dec = pi - dec; ra += pi; }
      if (dec < -piovertwo) { dec = -pi + dec; ra += pi; }
      if (ra > pi) { ra -= 2 * pi; }
      if (ra < -pi) { ra += 2 * pi; }
    }

    double x,y;
    bool show = transform.MapXY(image, params, dec, -ra, x, y);
    if (show) {
      int ix = int(floor(x));
      int iy = int(floor(y));
      if (0 <= ix && ix < width && 0 <= iy && iy < height) {
	int index = iy * width + ix;
	double value = 255 * convertmag(mag);
	Rgb current;
	image.GetRgb(index, current);
	if (value > current.r || value > current.g || value > current.b) {
	  // Adjust for BV color index
	  double r = value;
	  double g = value;
	  double b = value;
	  if (dobv) {
	    if (bv > 0) {
	      if (bv > bvfact) { bv = bvfact; }
	      g *= (bvfact - bv)/bvfact;
	      b *= (bvfact - bv)/bvfact;
	    } else if (bv < 0) {
	      if (bv < -bvfact) { bv = -bvfact; }
	      r *= (bvfact + bv)/bvfact;
	      g *= (bvfact + bv)/bvfact;
	    }
	  }
	  image.SetRgb(index, Rgb(Color(r), Color(g), Color(b)));
	}
      }
    }
    i++;
  }
  fclose(starfile);
}

class SkyMapper : public MapFunction
{
public:
  SkyMapper(Transform& transform_, TransformParams& params_, Rgb& sky_, Rgb& bg_)
    : transform(transform_), params(params_), sky(sky_), bg(bg_) {};
    
  double Scale(int width, int height)
  {
    return transform.BasicScale(width,height)/params.scale;
  }

  void InitY(double y) {
    transform.SetY(y);
  }

  Rgb operator()(double x, double y) const {
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    double phi = 0.0, lambda = 0.0;
    if (transform.Project(params,x,y,x0,y0,z0,phi,lambda)) {
      return bg;
    } else {
      return sky;
    }
  }
private:
  Transform& transform;
  TransformParams& params;
  Rgb& sky;
  Rgb& bg;
};

const char* usage = 
"Usage: %s [-bg r:g:b][-w n][-h n][-adjust][-f filename][-out filename]\n"
"          [-back filename][-year x][-magfact x][-topmag x][-bv][-bvfact x]\n"
"          [-tilt x][-turn x][-dec x][-ra x]\n"
"          [-scale x][-x x][-y x][-z x]\n"
"          [-a x][-date x][-time x][-p x][-aw x]\n"
"          [-conic x][-conicr x][-xoff x][-yoff x]\n"
"          [-grid][-gridx x][-gridy y][-gridcolor r:g:b]\n"
"          latlong|equalarea|sinusoidal|sinusoidal2|mollweide|mercator|azimuthal\n"
"          rectilinear|stereographic|gnomonic|perspective|bonne|hammer|lambert <infile> [<outfile>]\n";

int main(int argc, char* argv[])
{
  const char* outfilename = NULL;
  const char* backfilename = NULL;
  const char* starlist = NULL;
  int width = 800;
  int height = 600;

  bool adjust = false;
  bool addgrid = false;
  TransformParams params;
  params.gridcolor = Rgb(0,0,60);

  Rgb bg(0,0,0);

  const char* type = "perspective";
  bool typeset = false;

  // array of pointers. Dynamic allocation is less typing than
  // stack allocating everything.
  OptSpec *optspecs [] = {
    // General options
    new IntSpec ("-w", width),
    new IntSpec ("-h", height),
    new BoolSpec ("-adjust", adjust),
    new StringSpec ("-f", starlist),
    new StringSpec ("-out", outfilename),
    new StringSpec ("-back", backfilename),
    new RgbSpec ("-bg", bg),
    new DoubleSpec("-year", year),
    new DoubleSpec("-magfact", magfact),
    new DoubleSpec("-topmag", topmag),
    new BoolSpec("-bv", dobv),
    new DoubleSpec("-bvfact", bvfact),

    // First the params options
    new DoubleSpec ("-tilt", params.tilt, 2 * pi / 360.0 ),
    new DoubleSpec ("-turn", params.turn, 2 * pi / 360.0 ),
    new DoubleSpec ("-dec", params.lat, 2 * pi / 360.0 ),
    new DoubleSpec ("-ra", params.lon, -2 * pi / 24.0 ),
    new DoubleSpec ("-scale", params.scale),
    new DoubleSpec ("-x", params.x),
    new DoubleSpec ("-y", params.y),
    new DoubleSpec ("-z", params.z),
    new DoubleSpec ("-a", params.a),
    new DoubleSpec ("-date", params.date),
    new DoubleSpec ("-time", params.time, 2 * pi / 24.0),
    new DoubleSpec ("-p", params.p, 2 * pi / 360.0),
    new DoubleSpec ("-aw", params.aw, 2 * pi / 360.0),
    new DoubleSpec ("-conic", params.conic),
    new DoubleSpec ("-conicr", params.conicr),
    new DoubleSpec ("-xoff", params.xoff),
    new DoubleSpec ("-yoff", params.yoff),

    // Grid options
    new BoolSpec ("-grid", addgrid),
    new IntSpec ("-gridx", params.gridx),
    new IntSpec ("-gridy", params.gridy),
    new RgbSpec ("-gridcolor", params.gridcolor),

    NULL // Finish up
  };

  char *progname = argv[0];
  argc--;argv++;
  if (argc == 1 && strcmp(argv[0],"-v") == 0) {
    fprintf(stderr,"MMPS Version %s\n", VERSION);
    exit(0);
  }
  while (argc > 0) {
    int result = OptSpec::ReadOpts(argc, argv, optspecs);
    if (result == -1) {
      error(usage, progname);
    } else if (result != 1) {
      if (typeset) { error(usage, progname); };
      type = argv[0];
      typeset = true;
    }
    argc--; argv++;
  }

  Transform* transform = Transform::GetTransform(type);
  if (transform == NULL) {
    error(usage, progname);
  }

  if (adjust) {
    transform->AdjustSize(width, height, params);
  }
  Image outImage(width,height);
  transform->Init(params);
  if (backfilename != NULL) {
    Image backImage;
    backImage.Read(backfilename);
    outImage.Copy(backImage);
  } else {
    outImage.Fill(bg);
  }
  if (addgrid) {
    transform->AddLines(outImage, params);
  }
  AddStars(starlist, outImage, *transform, params);
  outImage.Write(outfilename);
  return 0;
}

