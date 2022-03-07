// $Revision: 1.3 $
// equations.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "equations.h"

bool Equation::FindRoot (double t0, double t1, double epsilon, double& t) const
{
  bool result = false;
  double p0 = (*this)(t0);
  double p1 = (*this)(t1);
  if (p0 > p1) {
    double temp;
    temp = t0; t0 = t1; t1 = temp;
    temp = p0; p0 = p1; p1 = temp;
  };
  if (p0 > 0 || p1 < 0) {
    fprintf(stderr, "No root for %f %f %f %f\n", t0, t1, p0, p1);
  } else {
    while (t1 - t0 > epsilon) {
      double t2 = (t0 + t1) / 2.0;
      double p2 = (*this)(t2);
      if (p2 <= 0.0) {
	t0 = t2;
      } else {
	t1 = t2;
      }
    }
    t = t0;
    result = true;
  }
  return result;
}
  
// The angle the sun is round from perihelion at time m, measured from
// perihelion.

double perihelionangle (double m)
{
  double e = eccentricity;
  double r = m + 2*e*sin(m) + 5*sqr(e)*sin(2*m)/4 + 
             cube(e)*(-sin(m)/4 + 13*sin(3*m)/12);
  return r;
}

// Return the time, relative to the spring equinox, of the perihelion
// Since this is constant, the result now in equations.h
double getperiheliontime ()
{
  double theta0 = 0;
  double theta1 = 2 * pi;
  double epsilon = 1e-8;
  while (theta1 - theta0 > epsilon) {
    double theta2 = (theta0 + theta1) / 2.0;
    double p = perihelionangle(theta2) - equinoxangle;
    if (p <= 0.0) {
      theta0 = theta2;
    } else {
      theta1 = theta2;
    }
  }
  printf ("Perihelion time: %f\n", -theta1);
  return -theta1;
}

// The value of the equation of time
// mean time + equation = apparent time
double equationoftime(double date)
{
  double p = perihelion;
  double m = date - periheliontime;
  double s = 
    -591.7 * sin(2 * (m+p))
    -459.6 * sin(m) +
     +19.8 * sin(m + 2 * p)
     -19.8 * sin(3 * m + 2 * p)
     -12.8 * sin(4 * (m + p))
      -4.8 * sin(2 * m)
      +0.9 * sin(3 * m + 4 * p)
      -0.9 * sin(5 * m + 4 * p)
      -0.5 * sin(4 * m + 2 * p)
      -0.4 * sin(6 * (m + p));
    
  // fprintf(stderr, "EOT: %f\n", s / 60.0);
  return s; // seconds
}
  
// The angle the sun is round from the vernal equinox at date,
// also measured from the vernal equinox.
double sunangle (double date)
{
  // Date is measured from vernal equinox
  // Need time from perihelion to equinox

  // Add the time from the perihelion and equinox
  double t = date - periheliontime;
  // r is the angle of the sun, from the perihelion
  double r = perihelionangle(t);
  // So subtract the angle between perihelion and equinox
  return r - equinoxangle;
}

double sunheight (double date)
{
  return sin(sunangle(date)) * sin(inclination);
}

double sundec (double date)
{
  return asin(sunheight(date));
}

double sunaltitude (double date, double phi, double time)
{
  // Compute altitude of sun at the given date and time,
  // on the zero meridian at latitude phi
  double delta = asin(sunheight(date)); // Declination of the sun

  // Rotate lat = phi, long = 0 round by time
  // Find x,y,z coords on the sphere - no need to use an ellipsoid, as the
  // tangents at a location are the same.
  double x = cos(time) * cos(phi);
  double y = sin(time) * cos(phi);
  double z = sin(phi);

  // The sun vector - the sun is at -x
  double sunx = -cos(delta);
  double suny = 0;
  double sunz = sin(delta);

  // Dot product to get angle to normal
  double alt = pi/2.0 - acos(x * sunx + y * suny + z * sunz);
  return alt;
}

double sunrise (double date, double phi)
{
  double t0 = 0;
  double t1 = pi;
  double epsilon = 1e-6;
  while (t1 - t0 > epsilon) {
    double t = (t1 + t0) / 2.0;
    double h = sunaltitude(date + t/daysinyear, phi, t);
    //printf ("%f %f %f %f\n", t0, t1, t, h);
    if (h < sunriseangle) {
      t0 = t;
    } else {
      t1 = t;
    }
  }
  return t0;
}
