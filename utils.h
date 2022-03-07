// $Revision: 1.4 $
// utils.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined UTILS_H
#define UTILS_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

static inline double roundint(double x)
{
  return floor (x + 0.5);
}

static const unsigned int MESSAGESIZE = 1024;

void error(const char *s);
void syserror(const char *s);

template <class T>
void error(const char *s, T t)
{
  char message[MESSAGESIZE];
  snprintf(message, MESSAGESIZE, s, t);
  error(message);
}

template <class T>
void syserror(const char *s, T t)
{
  char message[MESSAGESIZE];

  snprintf(message, MESSAGESIZE, s, t);
  syserror(message);
}

class OptSpec 
{
public:
  OptSpec(const char *optstring);
  virtual int ReadOpt(int &argc, char** &argv) = 0;
  const char *optstring;
  static int ReadOpts (int &argc, char** &argv, OptSpec* optspecs[]);
  virtual ~OptSpec() {};
 private:
  OptSpec(const OptSpec&);
  OptSpec &operator=(const OptSpec&);
};

class BoolSpec : public OptSpec
{
 public:
  BoolSpec (const char *optstring, bool &value);
  int ReadOpt(int &argc, char** &argv);
 private:
  bool &value;
};

class StringSpec : public OptSpec
{
 public:
  StringSpec (const char *optstring, const char* &value);
  int ReadOpt(int &argc, char** &argv);
 private:
  const char* &value;
};

class IntSpec : public OptSpec
{
 public:
  IntSpec (const char *optstring, int &value);
  int ReadOpt(int &argc, char** &argv);
 private:
  int &value;
};

class DoubleSpec : public OptSpec
{
 public:
  DoubleSpec (const char *optstring, double &value, double conversion = 1.0);
  int ReadOpt(int &argc, char** &argv);
 private:
  double &value;
  double conversion;
};

class CmdSpec : public OptSpec
{
 public:
  CmdSpec (const char *optstring, int nargs);
  int ReadOpt(int &argc, char** &argv);
  static const std::vector<std::vector<std::string> > &getcmdlist() { return cmdlist; }
 private:
  int nargs;
  static std::vector<std::vector<std::string> > cmdlist;
};
#endif
