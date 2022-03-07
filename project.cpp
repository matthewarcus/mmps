// $Revision: 1.6 $
// project.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "constants.h"
#include "image.h"
#include "transform.h"
#include "version.h"

const char* usage = 
"Usage: %s [-bg r:g:b][-w n][-h n][-adjust][-f filename][-out filename]\n"
"          [-back filename][-tilt x][-turn x][-rotate x][-lat x][-long x]\n"
"          [-scale x][-x x][-y x][-z x][-ox x][-oy x][-oz x]\n"
"          [-a x][-date x][-time x][-p x][-aw x][-radius x]\n"
"          [-conic x][-conicr x][-xoff x][-yoff x][-sun]\n"
"          [-grid][-gridx x][-gridy y][-gridoff x][-gridcolor r:g:b][-i]\n"
"          latlong|equalarea|sinusoidal|sinusoidal2|mollweide|mercator|cylindrical|azimuthal|\n"
"          rectilinear|orthographic|stereographic|gnomonic|perspective|bonne|hammer <infile> [<outfile>]\n";

int main(int argc, char* argv[])
{
  Image inImage;
  const char *infilename = "images/earth.ppm";
  const char *outfilename = NULL;
  const char *backfilename = NULL;
  int width = 800;
  int height = 600;
  bool adjust = false;

  TransformParams params;
  // Multi image stuff
  double xinc = 0.0;
  double yinc = 0.0;
  double zinc = 0.0;
  double turninc = 0.0;
  double tiltinc = 0.0;
  double latinc = 0.0;
  double loninc = 0.0;
  double timeinc = 0.0;
  double dateinc = 0.0;
  int loop = 1;

  bool invert = false; // Apply the inverse mapping

  // array of pointers. Dynamic allocation is less typing than
  // stack allocating everything.
  OptSpec *optspecs [] = {
    // General options
    new RgbSpec ("-bg", params.background),
    new IntSpec ("-w", width),
    new IntSpec ("-h", height),
    new BoolSpec ("-adjust", adjust),
    new StringSpec ("-f", infilename),
    new StringSpec ("-out", outfilename),
    new StringSpec ("-back", backfilename),
    new BoolSpec ("-i", invert),

    // First the projection parameters
    new DoubleSpec ("-tilt", params.tilt, 2 * pi / 360.0 ),
    new DoubleSpec ("-turn", params.turn, 2 * pi / 360.0 ),
    new DoubleSpec ("-rotate", params.rotate, 2 * pi / 360.0 ),
    new DoubleSpec ("-lat", params.lat, 2 * pi / 360.0 ),
    new DoubleSpec ("-long", params.lon, 2 * pi / 360.0 ),
    new DoubleSpec ("-scale", params.scale),
    new DoubleSpec ("-radius", params.radius, 2 * pi / 360.0 ),
    new DoubleSpec ("-x", params.x),
    new DoubleSpec ("-y", params.y),
    new DoubleSpec ("-z", params.z),
    new DoubleSpec ("-ox", params.ox),
    new DoubleSpec ("-oy", params.oy),
    new DoubleSpec ("-oz", params.oz),
    new DoubleSpec ("-a", params.a),
    new DoubleSpec ("-date", params.date),
    new DoubleSpec ("-time", params.time, 2 * pi / 24.0),
    new DoubleSpec ("-p", params.p, 2 * pi / 360.0),
    new DoubleSpec ("-aw", params.aw, 2 * pi / 360.0),
    new DoubleSpec ("-conic", params.conic),
    new DoubleSpec ("-conicr", params.conicr),
    new DoubleSpec ("-xoff", params.xoff),
    new DoubleSpec ("-yoff", params.yoff),
    new BoolSpec ("-sun", params.sun),

    // Grid and dial options
    // These accumulate in a list and are 'executed' in order.
    new CmdSpec ("-color", 1), // Set color
    new CmdSpec ("-gridx", 1), // Set grid coordinates
    new CmdSpec ("-gridy", 1),
    new CmdSpec ("-gridoff", 1), // Angular offset
    new CmdSpec ("-dlat", 1), // Set dial origin
    new CmdSpec ("-dlong", 1),

    new CmdSpec ("-grid", 0),
    new CmdSpec ("-localhours",0),
    new CmdSpec ("-temporaryhours",0),
    new CmdSpec ("-altitudes",0),
    new CmdSpec ("-analemma",0),
    new CmdSpec ("-tropics",0),
    new CmdSpec ("-dateline",1),
    new CmdSpec ("-datetime",1),

    // Backwards compatible - set default initial color
    new RgbSpec ("-gridcolor", params.gridcolor),

    // Loop options
    new IntSpec ("-loop", loop ),
    new DoubleSpec ("-tiltinc", tiltinc, 2 * pi / 360.0 ),
    new DoubleSpec ("-turninc", tiltinc, 2 * pi / 360.0 ),
    new DoubleSpec ("-latinc", latinc, 2 * pi / 360.0 ),
    new DoubleSpec ("-longinc", loninc, 2 * pi / 360.0 ),
    new DoubleSpec ("-xinc", xinc),
    new DoubleSpec ("-yinc", yinc),
    new DoubleSpec ("-zinc", zinc),
    new DoubleSpec ("-dateinc", dateinc, 2 * pi / 365.25),
    new DoubleSpec ("-timeinc", timeinc, 2 * pi / 24.0),

    NULL // Finish up
  };
  
  char *progname = argv[0];
  argc--;argv++;
  if (argc == 1 && strcmp(argv[0],"-v") == 0) {
    fprintf(stderr,"MMPS Version %s " __DATE__ "\n", VERSION);
    exit(0);
  }

  Transform* transform = NULL;

  while (argc > 0) {
    int result = OptSpec::ReadOpts(argc, argv, optspecs);
    if (result == -1) {
      // Assume error has been reported by ReadOpts
      error(usage, progname);
    } else if (result != 1) {
      // Option not recognized
      // Is it a transform>?
      Transform *t = Transform::GetTransform(argv[0]);
      if (t == NULL) {
	fprintf(stderr,"Unrecognized option %s\n",argv[0]);
	error(usage, progname);
      } else {
	if (transform == NULL) {
	  // Found a transform
	  transform = t;
	} else {
	  // Oops, two transforms
	  fprintf(stderr,"Can only specify one projection type\n");
	  error(usage, progname); 
	}
      }
    }
    argc--; argv++;
  }

  if (transform == NULL) {
    fprintf(stderr,"Need to specify a projection type\n");
    error(usage, progname);
  }

  inImage.Read(infilename);

  if (adjust) {
    transform->AdjustSize(width, height, params);
  }
  Image outImage(width, height);
  Image backImage;
  if (backfilename != NULL) {
    backImage.Read(backfilename);
  }
  for (int i = 0; i < loop; i++) {
    if (backfilename != NULL) {
      outImage.Copy(backImage);
      params.noback = true;
    }
    transform->Init(params);
    if (invert) {
      transform->TransformImageInv(inImage, outImage, params);
      // No commands for inverse transformation yet
    } else {
      transform->TransformImage(inImage, outImage, params);
      transform->DoCommands(outImage, params, CmdSpec::getcmdlist());
    }
    if (loop > 1) {
      char fname[256];
      sprintf(fname, "%s%04d.ppm", outfilename, i);
      outImage.Write(fname);
      fprintf(stderr, "Written %s\n", fname);
    } else {
      outImage.Write(outfilename);
    }      
    params.turn += turninc;
    params.tilt += tiltinc;
    params.lat += latinc;
    params.lon += loninc;
    params.x += xinc;
    params.y += yinc;
    params.z += zinc;
    params.time += timeinc;
    params.date += dateinc;
  }
  delete(transform);
  return 0;
}
