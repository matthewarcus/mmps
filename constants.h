// $Revision: 1.1 $
// constants.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined CONSTANTS_H
#define CONSTANTS_H

const long double pi = 3.14159265358979323846;
const long double twopi = 2.0 * pi;
const long double oneoverpi = 1.0 / pi;
const long double twooverpi = 2.0 / pi;
const long double piovertwo = pi / 2.0;
const long double oneovertwopi = 1.0 / (2.0 * pi);

const double secsperday = 86400;

inline double radians(double deg) { return deg * pi / 180.0; };
inline double hourstoradians(double hours) { return hours * pi / 12.0; };
inline double degrees(double rad) { return rad * 180.0 / pi; };
inline double hours(double rad) { return rad * 12.0 / pi; };

inline double sqr(double x) { return x * x; };
inline double cube(double x) { return x * x * x; }
inline double frac(double x) { return x - floor(x); };
#endif
