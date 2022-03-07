// $Revision: 1.1 $
// vector.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined VECTOR_H
#define VECTOR_H

#include <iosfwd>
#include <math.h>

double radians(double deg) {
  return deg * M_PI / 180;
}

double degrees(double rad) {
  return rad * 180 / M_PI;
}

class spherical {
 public:
  double phi;
  double lambda;
  spherical(double phi_, double lambda_)
    : phi(phi_),lambda(lambda_) {}
  friend std::ostream &operator<<(std::ostream &s, const spherical &p) {
    s << "[" << degrees(p.phi) << ", " << degrees(p.lambda) << "]";
    return s;
  }
};

class vector {
 public:
  double x;
  double y;
  double z;
  vector(double x_, double y_, double z_) : x(x_),y(y_),z(z_) {};

  vector(double phi, double lambda) 
    : x(cos(lambda)*cos(phi)),y(sin(lambda)*cos(phi)),z(sin(phi)) {}

  friend double operator* (const vector &v1, const vector &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }
  friend vector operator^ (const vector &v1, const vector &v2) {
    return vector(v1.y * v2.z - v1.z * v2.y,
		  v1.z * v2.x - v1.x * v2.z,
		  v1.x * v2.y - v1.y * v2.x);
  }
  friend vector operator+ (const vector &v1, const vector &v2) {
    return vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
  }
  friend vector operator- (const vector &v1, const vector &v2) {
    return vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }
  friend vector operator *(double d, const vector &v) {
    return vector(v.x * d, v.y *d, v.z * d);
  }
  friend vector operator *(const vector &v, double d) {
    return vector(v.x * d, v.y *d, v.z * d);
  }
  friend vector operator /(const vector &v, double d) {
    return vector(v.x/d, v.y/d, v.z/d);
  }
  friend vector operator -(const vector &v) {
    return -1 * v;
  }
  friend std::ostream &operator<<(std::ostream &s, const vector &v) {
    s << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    return s;
  }
  friend double length(const vector &v) {
    return sqrt(v * v);
  }
  friend vector norm(const vector &v) {
    return v/length(v);
  }
  operator spherical() {
    return spherical(asin(z),atan2(y,x));
  }
};

#endif
