// $Revision: 1.3 $
// equations.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined EQUATIONS_H
#define EQUATIONS_H

#include "constants.h"

// 16 minutes half-diameter + 34 minutes refraction

// The earth
const double eccentricity = 0.017;

// It's orbit
const double inclination = radians(23.45);
const double perihelion = radians(12.25);// Angle from winter solstice to perihelion
const double equinoxangle = pi/2.0 - perihelion; // Angle from equinox to perihelion
const double periheliontime = -1.323859; // Time from equinox to perihelion

// Sun declination for sunrise/sunset
const double sunriseangle = radians(-(16.0 + 34)/60.0);
const double daysinyear = 365.242;

inline double days(double t) { return t * daysinyear / (2 * pi); };

double equationoftime(double date);
double sunheight (double date); // sin of the declination at date
double sundec (double date); // sin of the declination at date
double sunrise (double date, double phi);

class Equation
{
 public:
  virtual double operator() (double t) const = 0;
  bool FindRoot (double t0, double t1, double epsilon, double& t) const;
  virtual ~Equation() {};
};

#endif
