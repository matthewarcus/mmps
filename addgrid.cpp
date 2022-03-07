// $Revision: 1.1 $
// addgrid.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "image.h"
#include "transform.h"

const char* usage = 
"Usage: %s [-grid x][-gridx x][-gridy x] file\n";

int main(int argc, char* argv[])
{
  Image inImage;
  char* infilename = NULL;

  double gridx = 10.0;
  double gridy = 10.0;
  Rgb gridcolor(0,0,0);

  char *progname = argv[0];
  argc--;argv++;
  while (argc > 0 && argv[0][0] == '-') {
    if (strcmp(argv[0], "-grid") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      gridx = strtod(argv[0],NULL);
      gridy = gridx;
    } else if (strcmp(argv[0], "-gridx") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      gridx = strtod(argv[0],NULL);
    } else if (strcmp(argv[0], "-gridy") == 0) {
      argc--; argv++;
      if (argc == 0) { error(usage, progname); };
      gridy = strtod(argv[0],NULL);
    } else {
      error(usage, progname);
    }
    argc--; argv++;
  }
  if (argc != 1) {
    error(usage, progname);
  }

  infilename = argv[0];

  inImage.Read(infilename, true);
  inImage.AddGrid(gridx, gridy, gridcolor);
  inImage.Close();
  return 0;
}

