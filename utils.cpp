// $Revision: 1.3 $
// utils.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "utils.h"

using namespace std;

void error(const char *s)
{
  fprintf(stderr, "%s\n", s);
  exit(-1);
}

void syserror(const char *s)
{
  perror(s);
  fprintf(stderr, "Exiting\n");
  exit(-1);
}

OptSpec::OptSpec(const char *optstring_)
  : optstring(optstring_) 
{}

int OptSpec::ReadOpts (int &argc, char** &argv, OptSpec* optspecs[])
{
  int i = 0;
  while (optspecs[i] != NULL) {
    int result = optspecs[i]->ReadOpt(argc, argv);
    if (result != 0) {
      return result;
    }
    i++;
  }
  return 0;
}

BoolSpec::BoolSpec (const char *optstring_, bool &value_)
  : OptSpec(optstring_), value(value_) {}

int BoolSpec::ReadOpt(int &argc, char** &argv) {
  (void(argc));
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    value = true;
    return 1;
  }
}

StringSpec::StringSpec (const char *optstring_, const char *&value_)
  : OptSpec(optstring_), value(value_) 
{}

int StringSpec::ReadOpt(int &argc, char** &argv) {
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    argc--; argv++;
    if (argc == 0) {
      fprintf(stderr,"Need a parameter for %s\n", optstring);
      return -1;
    } else {
      value = argv[0];
      return 1;
    }
  }
}

IntSpec::IntSpec (const char *optstring_, int &value_)
  : OptSpec(optstring_), value(value_) 
{}

int IntSpec::ReadOpt(int &argc, char** &argv) {
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    argc--; argv++;
    if (argc == 0) {
      fprintf(stderr,"Need a parameter for %s\n", optstring);
      return -1;
    } else {
      // Should do some error checking here
      value = strtol(argv[0],NULL,0);
      return 1;
    }
  }
}

DoubleSpec::DoubleSpec (const char *optstring_, double &value_, double conversion_)
  : OptSpec(optstring_), value(value_), conversion(conversion_) 
{};

int DoubleSpec::ReadOpt(int &argc, char** &argv) {
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    argc--; argv++;
    if (argc == 0) {
      fprintf(stderr,"Need a parameter for %s\n", optstring);
      return -1;
    } else {
      // Should do some error checking here
      value = conversion * strtod(argv[0], NULL);
      return 1;
    }
  }
}

vector<vector<string> > CmdSpec::cmdlist;

CmdSpec::CmdSpec (const char *optstring_, int nargs_)
  : OptSpec(optstring_),nargs(nargs_) {}

int CmdSpec::ReadOpt(int &argc, char** &argv) {
  if (strcmp(argv[0], optstring) != 0) {
    return 0;
  } else {
    vector<string>cmd;
    // Skip the initial '-'
    assert(*argv[0] == '-');
    cmd.push_back(argv[0]+1);
    for (int i = 0; i < nargs; i++){
      argc--;argv++;
      if (argc == 0) {
	fprintf(stderr,"Need %d parameter for %s\n", nargs, optstring);
	return -1;
      } else {
	cmd.push_back(argv[0]);
      }
    }
    cmdlist.push_back(cmd);
    return 1;
  }
}
