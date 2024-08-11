// $Revision: 1.5 $
// transform.cpp
// (c) 2004-2022 Matthew Arcus

// MIT License

#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "utils.h"
#include "constants.h"
#include "matrix.h"
#include "image.h"
#include "equations.h"
#include "transform.h"

using namespace std;

const double nightdim = 0.5;

// Apply matrix to phi, lambda, and put the resulting
// cartesion coords in x,y,z.
void convertlatlong(double &phi, double &lambda, 
		    double& x, double& y, double& z,
		    const Matrix3& m)
{
  x = cos(lambda) * cos(phi);
  y = sin(lambda) * cos(phi);
  z = sin(phi);
  if (!m.IsIdentity()) {
    m.Apply(x,y,z);
    phi = asin(z);
    lambda = atan2(y,x);
  }
}

double distance(double phi0, double lambda0, double phi1, double lambda1)
{
  double x0 = cos(lambda0) * cos(phi0);
  double y0 = sin(lambda0) * cos(phi0);
  double z0 = sin(phi0);
  double x1 = cos(lambda1) * cos(phi1);
  double y1 = sin(lambda1) * cos(phi1);
  double z1 = sin(phi1);
  double d = acos((x0*x1 + y0*y1 + z0*z1));
  return d;
}

inline double cot(double theta) { return 1.0/tan(theta); }

TransformParams::TransformParams()
  : tilt(0.0), turn(0.0), rotate(0.0), lat(0.0), lon(0.0),
    scale(1.0), radius(0.0), a(1.0), aw(radians(20)),
    x(8.0), y(0.0), z(0.0), ox(1.0), oy(1.0), oz(1.0),
    sun(false), p(0.0), time(0.0), date(0.0), noback(false),
    conic(1.0), conicr(0.0), xoff(0.0), yoff(0.0), background(Rgb(0,0,0)),
    gridx(15), gridy(10), gridcolor(Rgb(0,0,0)), gridoff(0),
    analemma(false), analemmacolor(Rgb(255,0,0)),
    dial(false), dialcolor(Rgb(255,255,0)),
    dlat(0.0),dlon(0.0)
{
}

void Transform::Init(const TransformParams& params)
{
  m.Identity();
  m = m * Matrix3(ZAxis, params.turn);
  m = m * Matrix3(XAxis, params.tilt);
  m = m * Matrix3(YAxis, params.lat);
  m = m * Matrix3(ZAxis, params.lon);

  // Construct the inverse matrix
  minv.Identity();
  minv = minv * Matrix3(ZAxis,-params.lon);
  minv = minv * Matrix3(YAxis,-params.lat);
  minv = minv * Matrix3(XAxis,-params.tilt);
  minv = minv * Matrix3(ZAxis,-params.turn);

  // Construct the sun parameters
  if (!params.sun) {
    sunx = 0.0;
    suny = 0.0;
    sunz = 0.0;
  } else {
    // We will use the astronomical day, starting at midnight (at Greenwich)
    // And start the year at the spring equinox
    // Need the vector to the sun.
    double date = params.date-80;
    while (date < 0) date += 365;
    while (date >= 365) date -= 365;
    date = 2 * pi *date / 365;
    double sheight = sunheight(date);
    double q = sqrt(1-sqr(sheight));
    // params.time is mean time, convert to apparent time
    double eot = equationoftime(date) * 2.0 * pi/(60.0 * 60.0 * 24.0);
    // AT = MT + EOT
    double atime = params.time + eot;
    sunx = -cos(atime) * q;
    suny =  sin(atime) * q;
    sunz =  sheight;
    // fprintf (stderr, "Sun: %f %f\n", degrees(asin(sunz)), hours(atan2(suny,-sunx)));
  }
}

void Transform::SetData(const Image& inImage, Image& outImage,
			const TransformParams& params,
			double x, double y, double z,
			double phi, double lambda,
			unsigned int outIndex)	     
{
  unsigned int iw = inImage.Width();
  unsigned int ih = inImage.Height();

  double hh = ih * 0.5;
  double hw = iw * 0.5;
  double sh = ih * oneoverpi;
  double sw = iw * oneovertwopi;

  // Use unsigned so we don't have to test for negative indices
  unsigned int sx = unsigned(floor(hw + lambda * sw));
  // Clamp in case of rounding errors
  if (sx >= iw) sx = 0;

  // Use unsigned so we don't have to test for negative indices

  unsigned int sy = unsigned(floor(hh - phi * sh));
  if (sy >= ih) sy = ih-1;

  Rgb rgb;
  unsigned int inIndex = sy * iw + sx;
  inImage.GetRgb(inIndex, rgb);
  if (params.sun) {
    // Is where we are sunny?
    // Compute pi/2 - angle of the x-axis with the normal at the relevant point.
    // This should be asin(...) but we are only interested in small angles
    // so assume sin x = x
    double sunaltitude = sunx * x + suny * y + sunz * z;
#if 1
    if (sunaltitude < sunriseangle) {
      rgb.Dim(nightdim);
    } else {
      rgb.Dim(0.8 + 0.2 * sunaltitude);
    }
#else
    rgb.Dim(0.6 + 0.5 * sunaltitude * sunaltitude);
#endif
  }
  outImage.SetRgb(outIndex, rgb);
}

bool Transform::CheckPoint(const TransformParams &params, double phi, double lambda) const
{
  if (params.radius != 0 && 
      distance(phi,lambda,params.lat,params.lon) > params.radius) {
    return false;
  }
  return true;
}

void applyRotation(double r, double &x, double &y)
{
  double x1 =  x * cos(r) + y * sin(r);
  double y1 = -x * sin(r) + y * cos(r);
  //double x1 = -y;
  //double y1 = x;
  //fprintf(stderr,"%g %g:", x, y);
  x = x1;
  y = y1;
  //fprintf(stderr,"%g %g\n", x, y);
}
  
void Transform::TransformImage(const Image& inImage, Image& outImage,
			       const TransformParams& params)
{
  int ow = outImage.Width();
  int oh = outImage.Height();
  double orgx = 0.5 * ow;
  double orgy = 0.5 * oh;
  // Scale so that half width is 1
  double scaleFactor = BasicScale(ow,oh) / params.scale;
  int ox, oy; int outIndex = 0;
  //fprintf(stderr,"Rotation %g %g %g %g\n",
  //        params.rotate, params.xoff, params.yoff, scaleFactor);
  for (oy = 0; oy < oh; oy++) {
    double y = scaleFactor * (orgy - oy - 0.5) + params.yoff;
    if (params.rotate == 0) SetY(y);
    for (ox = 0; ox < ow; ox++, outIndex++) {
      // Really ought to just apply a 3x2 matrix to the point
      double x = scaleFactor * (ox + 0.5 - orgx) + params.xoff;
      double x1 = x;
      double y1 = y;
      if (params.rotate != 0) {
	applyRotation(-params.rotate,x1,y1);
	SetY(y1);
      }
      double x0 = 0.0, y0 = 0.0, z0 = 0.0;
      double phi = 0.0, lambda = 0.0;
      if (Project(params, x1, y1, x0, y0, z0, phi, lambda) &&
	  CheckPoint(params,phi, lambda)) {
	SetData(inImage, outImage, params, x0, y0, z0, phi, lambda, outIndex);
      } else if (!params.noback) {
	outImage.SetRgb(outIndex, params.background);
      }
    }
  }
}

void Transform::TransformImageInv(const Image& inImage, Image& outImage,
				  const TransformParams& params)
{
  int ow = outImage.Width();
  int oh = outImage.Height();
  double orgx = 0.5 * ow;
  double orgy = 0.5 * oh;
  // Now scan across the output image
  int ox, oy; int outIndex = 0;
  for (oy = 0; oy < oh; oy++) {
    for (ox = 0; ox < ow; ox++, outIndex++) {
      // Compute lat and long
      double phi = (orgy - oy - 0.5) * pi / oh;
      double lambda = (ox + 0.5 - orgx) * twopi / ow;
      
      // Compute the scaled x,y coordinates for <phi,lambda>
      double x = 0.0, y = 0.0;
      if (!CheckPoint(params,phi,lambda) ||
	  !ProjectInv(params, phi, lambda, x, y) ||
	  !SetDataInv(inImage, outImage, params, ox, oy, x, y, outIndex)) {
	if (!params.noback) {
	  outImage.SetRgb(outIndex, params.background);
	}
      }
    }
  }
}

bool Transform::SetDataInv(const Image& inImage, Image& outImage,
			   const TransformParams& params,
			   double ox, double oy, // Coordinates in output image
			   double x, double y,   // scaled coordinates in input image
			   unsigned int outIndex)	     
{
  (void(ox));
  (void(oy));
  int iw = inImage.Width();
  int ih = inImage.Height();
  double orgx = 0.5 * iw;
  double orgy = 0.5 * ih;

  // Scale so that half width is 1
  double scaleFactor = BasicScale(iw,ih) / params.scale;

  // This is what we are inverting
  //double x = scaleFactor * (ix - orgx) + params.xoff;
  //double y = scaleFactor * (orgy - iy) + params.yoff;
  //applyRotation(-params.rotate,x,y)
  
  if (params.rotate != 0) {
    applyRotation(params.rotate,x,y);
  }
  int ix = int(floor(orgx + (x - params.xoff) / scaleFactor));
  int iy = int(floor(orgy - (y - params.yoff) / scaleFactor));

  if (ix < 0 || ix >= iw || iy < 0 || iy >= ih) {
    return false;
  } else {
    Rgb rgb;
    unsigned int inIndex = iy * iw + ix;
    inImage.GetRgb(inIndex, rgb);
    outImage.SetRgb(outIndex, rgb);
    return true;
  }
}

bool Transform::MapXY(const Image& outImage, 
		      const TransformParams& params,
		      double phi, double lambda, double& x, double& y) const
{
  bool result = false;
  // Set x and y to where phi and lambda are mapped to
  // x, y are in image coordinates
  // Get projection coordinates for x and y
  if (CheckPoint(params,phi,lambda) &&
      ProjectInv(params, phi,lambda, x, y)) {
    // Now x and y are in 2pi scale
    int ow = outImage.Width();
    int oh = outImage.Height();
    double orgx = 0.5 * ow;
    double orgy = 0.5 * oh;
    double scaleFactor = BasicScale(ow,oh) / params.scale;
    if (params.rotate != 0) {
      applyRotation(params.rotate,x,y);
    }
    x = orgx + (x - params.xoff) / scaleFactor;
    y = orgy + (y - params.yoff) / -scaleFactor;
    result = true;
  }
  return result;
}

Transform* Transform::GetTransform(const char *type)
{
 Transform* transform = NULL;
  if (strcmp(type, "latlong") == 0) {
    transform = new LatLong();
  } else if (strcmp(type, "equalarea") == 0) {
    transform = new EqualArea();
  } else if (strcmp(type, "sinusoidal") == 0) {
    transform = new Sinusoidal();
  } else if (strcmp(type, "sinusoidal2") == 0) {
    transform = new Sinusoidal2();
  } else if (strcmp(type, "mollweide") == 0) {
    transform = new Mollweide();
  } else if (strcmp(type, "mercator") == 0) {
    transform = new Mercator();
  } else if (strcmp(type, "cylindrical") == 0) {
    transform = new Cylindrical();
  } else if (strcmp(type, "azimuthal") == 0) {
    transform = new Azimuthal();
  } else if (strcmp(type, "orthographic") == 0) {
    transform = new Rectilinear();
  } else if (strcmp(type, "rectilinear") == 0) {
    transform = new Rectilinear();
  } else if (strcmp(type, "stereographic") == 0) {
    transform = new Stereographic();
  } else if (strcmp(type, "gnomonic") == 0) {
    transform = new Gnomonic();
  } else if (strcmp(type, "perspective") == 0) {
    transform = new Perspective();
  } else if (strcmp(type, "bonne") == 0) {
    transform = new Bonne();
  } else if (strcmp(type, "hammer") == 0) {
    transform = new Hammer();
  }
  return transform;
}

void PolarTransform::Init (const TransformParams& params)
{
  Transform::Init(params);

  // Polar transforms have the north pole as the centre.
  // For consistency with other transforms, we rotate to set lat = 0
  m = Matrix3(YAxis,-piovertwo) * m;
  minv = minv * Matrix3(YAxis,piovertwo);

  k = params.conic;
}

bool PolarTransform::Project(const TransformParams& params,
			     double x0, double y0, 
			     double& x, double& y, double& z,
			     double& phi, double& lambda) const
{
  bool result = false;
  double r = sqrt(sqr(x0) + sqr(y0)) - params.conicr;
  lambda = params.conic * atan2(x0, -y0);
  if (r > 0.0 && lambda >= -pi && lambda <= pi) {
    GetPhi(r, phi);
    if (phi >= -piovertwo && phi <= piovertwo) {
      convertlatlong(phi,lambda,x,y,z,m);
      result = true;
    }
  }
  return result;
}

bool PolarTransform::ProjectInv(const TransformParams& params,
				double phi, double lambda,
				double& x, double& y) const
{
  bool result = false;
  // Set x and y to where phi and lambda are mapped to
  // x, y are in image coordinates
  // Get projection coordinates for x and y
  double x0,y0,z0;
  convertlatlong(phi,lambda,x0,y0,z0,minv);
  double r = 0.0;
  if (GetR(phi, r)) {
    x = (r + params.conicr) * sin(lambda / params.conic);
    y = -(r + params.conicr) * cos(lambda / params.conic);
    result = true;
  }
  return result;
}

bool CylindricalTransform::Project(const TransformParams& params,
				   double x0, double y0, 
				   double& x, double& y, double& z,
				   double& phi, double& lambda) const
{
  (void(params));
  (void(y0));
  bool result = false;

  phi = phi0;
  lambda = GetLong(x0);

  if (lambda >= -pi && lambda <= pi && 
      phi >= -piovertwo && phi <= piovertwo) {
    // Transform to new lat and long.
    // cartesian coords from latlong
    convertlatlong(phi,lambda,x,y,z,m);
    result = true;
  }
  return result;
}

bool CylindricalTransform::ProjectInv(const TransformParams& params,
				      double phi, double lambda, 
				      double& x, double& y) const
{
  (void(params));
  minv.ApplyLatLong(phi,lambda);
  return GetXY(phi,lambda, x, y);
}

bool SimpleTransform::Project(const TransformParams& params,
			      double x, double y, 
			      double& x1, double& y1, double& z1,
			      double& phi, double& lambda) const
{
  (void(params));
  bool result = ProjectSimple(x,y,phi,lambda);
  convertlatlong(phi,lambda,x1,y1,z1,m);
  return result;
}

bool SimpleTransform::ProjectInv(const TransformParams& params,
				 double phi, double lambda, double& x, double& y) const
{
  (void(params));
  minv.ApplyLatLong(phi,lambda);
  return ProjectInvSimple(phi,lambda,x,y);
}

double LatLong::GetMaxHeight(const TransformParams& params) 
{
  (void(params));
  return piovertwo;
}

double LatLong::GetLat(double y)
{
  return y;
}

double LatLong::GetLong(double x) const
{
  return x;
}

bool LatLong::GetXY(double phi, double lambda, double& x, double& y) const
{
  x = lambda;
  y = phi;
  return true;
}

// f(phi) = sin(phi)/k, g(phi) = 1, where k = sqr(cos(p))
double EqualArea::GetMaxHeight(const TransformParams& params)
{
  return 1 / (cos(params.p) * cos(params.p));
}

void EqualArea::Init(const TransformParams& params) 
{
  CylindricalTransform::Init(params);
  k = cos(params.p) * cos(params.p);
} 

double EqualArea::GetLat(double y)
{
  return asin(y * k);
}

double EqualArea::GetLong(double x) const
{
  return x;
}

bool EqualArea::GetXY(double phi, double lambda, double& x, double& y) const
{
  x = lambda;
  y = sin(phi)/k;
  return true;
}

// f(phi) = phi, g(phi) = cos(phi)
double Sinusoidal::GetMaxHeight (const TransformParams& params)
{
  (void(params));
  return piovertwo;
}

double Sinusoidal::GetLat(double y)
{
  // Latitude is just y
  m = 1/cos(y); // Cache the cos
  return y;
}

double Sinusoidal::GetLong(double x) const
{
  return x * m;
}

bool Sinusoidal::GetXY(double phi, double lambda, double& x, double& y) const
{
  x = lambda * cos(phi);
  y = phi;
  return true;
}

double Mercator::GetLat(double y)
{
  double k = exp(fabs(y));
  double phi = acos(2*k/(k*k+1));
  if (y < 0) phi = -phi;
  return phi;
}

double Mercator::GetLong(double x) const
{
  return x;
}

bool Mercator::GetXY(double phi, double lambda, double& x, double& y) const
{
  x = lambda;
  y = log((1 + sin(fabs(phi))) / cos(fabs(phi)));
  if (phi < 0) {
    y = -y;
  }
  return true;
}

double Cylindrical::GetLat(double y)
{
  return atan(y);
}

double Cylindrical::GetLong(double x) const
{
  return x;
}

bool Cylindrical::GetXY(double phi, double lambda, double& x, double& y) const
{
  x = lambda;
  y = tan(phi);
  return true;
}

void Azimuthal::AdjustSize (int& w, int& h, TransformParams& params)
{
  int w1 = (int (params.scale * h));
  if (w1 < w) { w = w1; }
}

void Azimuthal::GetPhi(double r, double& phi) const
{
  phi = pi * (0.5 - r);
}

bool Azimuthal::GetR(double phi, double& r) const
{
  r = 0.5 - phi / pi;
  return true;
}

void Rectilinear::GetPhi(double r, double& phi) const
{
  phi = acos(r);
}

bool Rectilinear::GetR(double phi, double& r) const
{
  bool result = false;
  if (phi >= 0.0) {
    r = cos(phi);
    result = true;
  }
  return result;
}

void Stereographic::GetPhi(double r, double& phi) const
{
  double sr = pow(r,k);
  phi = asin((1-sqr(sr)) / (1+sqr(sr)));
}

bool Stereographic::GetR(double phi, double& r) const
{
  r = cos(phi) / (1 + sin(phi));  
  r = pow(r,1/k);
  return true;
}

void Gnomonic::GetPhi(double r, double& phi) const
{
  phi = atan(1 / (2 * r));
}

bool Gnomonic::GetR(double phi, double& r) const
{
  bool result = false;
  if (phi > 0) {
    r = 1 / (2 * tan(phi));
    result = true;
  }
  return result;
}

void Hammer::Init(const TransformParams& params)
{
  Transform::Init(params);
}

bool Hammer::ProjectSimple(double x, double y, 
			   double& phi, double& lambda) const
{
  bool result = false;
  double z2 = 2-sqr(x)-sqr(2*y);
  double z = sqrt(z2);
  double t1 = 2*y*z;
  if (t1 >= -1.0 && t1 <= 1.0) {
    phi = asin(t1);
    lambda = 2 * atan2(x*z,z2 - 1);
    if (lambda >= -pi && lambda <= pi) {
      result = true;
    }
  }
  return result;
}

bool Hammer::ProjectInvSimple(double phi, double lambda, double& x, double& y) const
{
  double z = sqrt(1 + cos(phi)*cos(lambda/2));
  x = cos(phi)*sin(lambda/2)/z;
  y = sin(phi)/(2 * z);
  return true;
}

void Bonne::Init(const TransformParams& params)
{
  Transform::Init(params);
  p = params.p;
  cotp = cot(p);
}

bool Bonne::ProjectSimple(double x, double y, 
			  double& phi, double& lambda) const
{
  bool result = false;
  double rho = sqrt(sqr(x) + sqr(cotp - y));
  if (p > 0) {
    phi = cotp + p - rho;
    lambda = rho * atan2(x, cotp - y)/cos(phi);
  } else if (p < 0) {
    phi = cotp + p + rho;
    lambda = rho * atan2(x, y - cotp)/cos(phi);
  } else {
    // Degenerate case - the sinusoidal projection
    phi = y;
    lambda = x/cos(phi);
  }
  if (phi >= -piovertwo && phi <= +piovertwo &&
      lambda >= -pi && lambda <= pi) {
    result = true;
  }
  return result;
}

bool Bonne::ProjectInvSimple(double phi, double lambda, double& x, double& y) const
{
  if (p == 0.0) {
    x = lambda * cos(phi);
    y = phi;
  } else {
    double rho = cotp + p - phi;
    double e = lambda * cos(phi)/rho;
    x = rho * sin(e);
    y = cotp - rho * cos(e);
  }
  return true;
}

void Perspective::Init(const TransformParams& params) 
{
  Transform::Init(params);
  vx = params.x; vy = params.y; vz = params.z;
  rx = vx; ry = vy; rz = vz;
  m.Apply(rx,ry,rz);
  scalefactor = (vx + 1) * tan(params.aw/2);
  iscalefactor = 1/scalefactor;

  a = params.ox;
  b = params.oy;
  c = params.oz;
  a2 = sqr(a);
  b2 = sqr(b);
  c2 = sqr(c);
  ia2 = 1/a2;
  ib2 = 1/b2;
  ic2 = 1/c2;
}

bool Perspective::Project(const TransformParams& params,
			  double x0, double y0, 
			  double& x, double& y, double& z,
			  double& phi, double& lambda) const
{
  (void(params));
  bool result = false;
  // Apply our rotation to the point we are projecting to
  double x1 = -1;
  double y1 = scalefactor * x0 + vy;
  double z1 = scalefactor * y0 + vz; 
  m.Apply(x1,y1,z1);

  // Projecting from (rx, ry, rz) to point (x1, y1, z1)
  // Solve a quadratic obtained from equating line equation with r = 1
  double qa = a2 * sqr(rx-x1) + b2 * sqr(ry-y1) + c2 * sqr(rz-z1);
  double qb = 2 * (a2 * x1 * (rx-x1) + b2 * y1 * (ry-y1) + c2 * z1 * (rz-z1));
  double qc = a2 * sqr(x1) + b2 * sqr(y1) + c2 * sqr(z1) - 1;
  double qm = qb * qb - 4 * qa * qc;
  if (qm >=  0) {
    // Since qa is always positive, the + solution is nearest to the point
    // of view, so we always use that one.
    double k = (-qb + sqrt(qm)) / (2*qa);
    x = k * rx + (1-k) * x1;
    y = k * ry + (1-k) * y1;
    z = k * rz + (1-k) * z1;
    if (a == 1.0 && b == 1.0 && c == 1.0) {
      phi = asin(z);
      lambda = atan2(y, x);
    } else {
      // This is a point on the ellipsoid, so convert to lat long
      double r = sqrt(sqr(a2 * x) + sqr(b2 * y));
      phi = atan(c2 * z / r);
      lambda = atan2(b2*y, a2*x);
      // Now return the spherical x, y, z corresponding to phi, lambda
      x = cos(lambda) * cos(phi);
      y = sin(lambda) * cos(phi);
      z = sin(phi);
    }
    result = true;
  }
  return result;
}

bool Perspective::ProjectInv(const TransformParams& params, 
			     double phi, double lambda, double& x, double& y) const
{
  (void(params));
  bool result = false;
  // Find where phi, lambda project to on the ellipsoid
  double x0,y0,z0;
  double nx,ny,nz;
  if (a == 1.0 && b == 1.0 && c == 1.0) {
    // Just a sphere
    x0 = cos(lambda) * cos(phi);
    y0 = sin(lambda) * cos(phi);
    z0 = sin(phi);
    minv.Apply(x0,y0,z0);
    nx = x0; ny = y0; nz = z0;
  } else {
    // The normal vector
    nx = cos(lambda);
    ny = sin(lambda);
    nz = tan(phi);

    // The actual vector
    double r = sqrt(1/(ia2 * sqr(nx) + ib2 * sqr(ny) + ic2 * sqr(nz)));
    x0 = nx * r * ia2;
    y0 = ny * r * ib2;
    z0 = nz * r * ic2;
    // And apply rotations
    minv.Apply(x0,y0,z0);
    minv.Apply(nx,ny,nz);
  }

  // Test for visibility - dot product of nx,ny,nz and the line
  // towards to viewpoint must be positive
  double p = (vx-x0) * nx + (vy-y0) * ny + (vz-z0) * nz;
  if (p >= 0) {
    // Project from (vx, vy, vz) through (x0,y0,z0) to (-1, y, z)
    double t = (1 + x0)/(x0 - vx);
    x = t * vy + (1-t) * y0; // Note change of axes
    y = t * vz + (1-t) * z0;
    x = (x - vy) * iscalefactor;
    y = (y - vz) * iscalefactor;
    result = true;
  }
  return result;
}

void CylindricalTransform::AdjustSize (int& w, int& h, 
				       TransformParams& params)
{
  if (params.scale < 1.0) {
    w = int(w * params.scale);
    params.scale = 1.0;     
  }
  // If a parallel has been given, check it's valid
  if (params.p < 0 || params.p >= piovertwo) {
    error ("Invalid parallel");
  }
  double h0 = GetMaxHeight(params);
  if (h0 > 0.0) { // Return <= 0 for don't adjust
    int newh = int(h0 * w * params.scale / pi);
    if (newh < h) { h = newh; };
  }
}

// f(phi) = bt, g(phi) = (cos(t) + a)/(a+1), where sin(t) + at = k(a+1)sin(phi)
// and b + (a * pi)/2 = k*(a+1)

double Sinusoidal2::GetMaxHeight(const TransformParams& params)
{
  double k0, b0;
  CalcParams(params.a, params.p, k0, b0);
  return b0 * piovertwo;
}

void Sinusoidal2::CalcParams(double a, double phi, double& k, double& b)
{
  // Find the t value for phi
  double t1 = 0;
  double t2 = pi/2.0;
  double epsilon = 1e-6;
  double k0 = 0.0, b0 = 0.0;
  while (fabs(t2 - t1) > epsilon) {
    double t0 = (t1 + t2)/2.0;
    k0 = sqr((cos(t0) + a) / ((a+1)*cos(phi)));
    b0 = (k0 * (a+1))/(1 + a * pi / 2.0);
    double p0 = b0 * sin(t0) + a * t0 - k0 * (a+1) * sin(phi);
    if (p0 <= 0.0) t1 = t0;
    else t2 = t0;
  }
  k = k0; b = b0;
}

void Sinusoidal2::Init(const TransformParams& params)
{
  CylindricalTransform::Init(params);
  a = params.a;
  CalcParams(params.a, params.p, k, b);
}

double Sinusoidal2::GetLat(double y)
{
  double t = y/b;
  double phi = asin((b * (sin(t) + a * t)) / (k * (a+1)));
  m = (a+1) / (cos(t) + a);
  return phi;
}

double Sinusoidal2::GetLong(double x) const
{
  return m * x;
}

class Sinusoidal2Equation : public Equation
{
public: 
  double operator() (double t) const {
    // b(sin(t) + at) / k(a+1)-sin(phi);
    return b * (sin(t) + a * t) / (k * (a+1)) - sinphi;
  }
public:
  double b;
  double a;
  double k;
  double sinphi;
};
    
bool Sinusoidal2::GetXY(double phi, double lambda, double& x, double& y) const
{
  bool result = false;
  Sinusoidal2Equation equation;
  equation.b = b; equation.a = a; equation.k = k; equation.sinphi = sin(phi);
  double t = 0.0;
  if (equation.FindRoot(-pi/2, pi/2, 1e-5, t)) {
    y = b * t;
    x = lambda * (cos(t) + a) / (a+1);
    result = true;
  }
  return result;
}

// f(phi) = 1/a sin t, g(phi) = cos t, where t + 0.5 sin(2t) = 0.5 pi sin(phi)
double Mollweide::GetMaxHeight(const TransformParams& params)
{
  double a0;
  CalcParams(params.p, a0);
  return 1.0/a0;
}

void Mollweide::CalcParams(double phi, double& a)
{
  // Now we need to find a suitable t
  // Good old binary search - note we cavalierly fail to check if there's
  // a solution in the interval. Hell, we're going to terminate with something.
  double t = 0;     // phi = 0;
  double t1 = pi/2;  // phi = pi/2;
  double b = 0.5 * pi * sin(phi);
  double epsilon = 1e-6;
  while (fabs(t1 - t) > epsilon) {
    double t2 = (t + t1)/2.0;
    double p = t2 + 0.5 * sin (2 * t2) - b;
    if (p <= 0.0) t = t2; 
    else t1 = t2;
  }
  // t is our selected value
  a = pi / (4 * sqr(cos (t)/cos(phi)));
}

void Mollweide::Init(const TransformParams& params)
{
  CylindricalTransform::Init(params);
  double phi = radians(40.73333); //params.p;
  CalcParams(phi, a);
}

double Mollweide::GetLat(double y)
{
  // find t;
  double t = asin(y*a);
  m = 1/cos(t);
  return asin((t + 0.5 * sin(2 * t)) * twooverpi);
}
double Mollweide:: GetLong(double x) const
{
  return m * x;
}

class MollweideEquation : public Equation
{
public: 
  double operator() (double t) const {
    return (t + 0.5 * sin(2 * t)) * twooverpi - sinphi;
  }
public:
  double sinphi;
};
    
bool Mollweide::GetXY(double phi, double lambda, double& x, double& y) const
{
  bool result = false;
  MollweideEquation equation;
  equation.sinphi = sin(phi);
  double t = 0.0;
  if (equation.FindRoot(-pi/2, pi/2, 1e-5, t)) {
    y = sin(t)/a;
    x = lambda * cos(t);
    result = true;
  }
  return result;
}

class TransformPlotter : public LinePlotter
{
public:
  TransformPlotter(const Image& image_, const TransformParams& params_, 
	  const Transform& transform_)
    : image(image_), params(params_), transform(transform_) {}
 protected:
  const Image& image;
  const TransformParams& params;
  const Transform& transform;
};

// Draw a line of latitude at longitude lambda
class LatPlotter : public TransformPlotter
{
 public:
  LatPlotter(const Image& image_, const TransformParams& params_, 
	      const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), lambda(0) {}
  bool GetXY(double t, double&x, double&y) const
  {
    return transform.MapXY(image, params, t, lambda, x, y);
  }

 public:
  double lambda;
};

// Draw a temporary hour line
class TempPlotter : public TransformPlotter
{
 public:
  TempPlotter(const Image& image_, const TransformParams& params_, 
	      const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), time(0), phi(0), lambda(0) {}
  // t is the solar declination
  // phi,lambda is origin
  bool GetXY(double t, double&x, double&y) const
  {
    double delta = t;
    // time/angle between sunrise and noon and noon and sunset
    double tau = acos(-tan(phi)*tan(delta));
    //fprintf(stderr,"%g %g %g %g\n", t, phi, lambda, tau);
    double h = tau * (time - 12) / 6.0 + lambda;
    return transform.MapXY(image, params, delta, h, x, y);
  }

 public:
  double time;
  double phi;
  double lambda;
};

class AnalemmaPlotter : public TransformPlotter
{
 public:
  AnalemmaPlotter(const Image& image_, const TransformParams& params_, 
	      const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), time(0) {}
  bool GetXY(double t, double&x, double&y) const
  {
    // AT = MT + EOT
    // lambda = -AT (measure longitude to the east)
    double eot = twopi * equationoftime(t) / (24 * 60 * 60);
    double delta = sundec(t);
    
    return transform.MapXY(image, params, delta, -(time+eot), x, y);
  }

 public:
  double time;
};

// Draw a line of longitude at latitude phi
class LongPlotter : public TransformPlotter
{
 public:
  LongPlotter(const Image& image_, 
	      const TransformParams& params_, 
	      const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), phi(0) {}
  bool GetXY(double t, double&x, double&y) const
  {
    return transform.MapXY(image, params, phi, t, x, y);
  }
 public:
  double phi;
};

inline void RotateX(double theta, double&x, double&y, double& z)
{
  double x0 = x;
  double y0 = cos(theta) * y - sin(theta) * z;
  double z0 = sin(theta) * y + cos(theta) * z;
  x = x0; y = y0; z = z0;
}

inline void RotateY(double theta, double&x, double&y, double& z)
{
  double x0 = cos(theta) * x - sin(theta) * z;
  double y0 = y;
  double z0 = sin(theta) * x + cos(theta) * z;
  x = x0; y = y0; z = z0;
}

inline void RotateZ(double theta, double&x, double&y, double& z)
{
  double x0 = cos(theta) * x - sin(theta) * y;
  double y0 = sin(theta) * x + cos(theta) * y;
  double z0 = z;
  x = x0; y = y0; z = z0;
}

class CirclePlotter : public TransformPlotter
{
 public:
  CirclePlotter(const Image& image_, const TransformParams& params_, 
		const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), theta(0), phi(0){}
  bool GetXY(double t, double&x, double&y) const
  {
    double x0 = cos(t);
    double y0 = sin(t);
    double z0 = 0;

    RotateX(phi,x0,y0,z0);
    RotateY(phi,x0,y0,z0);
    RotateZ(theta,x0,y0,z0);

    double phi0 = asin(z0);
    double lambda0 = atan2(y0,x0);
    return transform.MapXY(image, params, phi0, lambda0, x, y);
  }
 public:
  double theta;
  double phi;
};

class FooPlotter : public TransformPlotter
{
 public:
  FooPlotter(const Image& image_, const TransformParams& params_, 
	     const Transform& transform_)
    : TransformPlotter (image_, params_, transform_), 
      theta(0), lambda(0), phi(0){}
  bool GetXY(double t, double&x, double&y) const
  {
    // Angle up from the pole
    double x0 = cos(t)*cos(theta);
    double y0 = sin(t)*cos(theta);
    double z0 = sin(theta);
    // Now rotate about y axis
    double x1 = sin(phi)*x0 + cos(phi)*z0;
    double y1 = y0;
    double z1 = -cos(phi)*x0 + sin(phi)*z0;
    double phi0 = asin(z1);
    if (phi0 < -inclination || phi0 > inclination) {
      return false;
    } else {
      double lambda0 = atan2(y1,x1);
      return transform.MapXY(image, params, phi0, lambda0+lambda, x, y);
    }
  }
 public:
  double theta;
  double lambda;
  double phi;
};

void Transform::AddLines (Image& image, const TransformParams& params) const
{
  CmdParams cmdparams;
  cmdparams.color = params.gridcolor;
  cmdparams.gridx = params.gridx;
  cmdparams.gridy = params.gridy;
  DrawGrid(image,params,cmdparams);
  //DrawDial(image,params,cmdparams);
}

void Transform::DoCommands (Image& image, 
			    const TransformParams& params,
			    const vector<vector<string> > cmdlist) const
{
  CmdParams cmdparams;
  cmdparams.color = params.gridcolor;
  for (size_t i = 0; i < cmdlist.size(); i++) {
    string cmd = cmdlist[i][0];
    if (cmd == "grid") {
      DrawGrid(image,params,cmdparams);
    } else if (cmd == "gridx") {
      cmdparams.gridx = atoi(cmdlist[i][1].c_str());
    } else if (cmd == "gridy") {
      cmdparams.gridy = atoi(cmdlist[i][1].c_str());
    } else if (cmd == "gridoff") {
      cmdparams.gridoff = atoi(cmdlist[i][1].c_str());
    } else if (cmd == "temporaryhours") {
      DrawTemporaryHours(image,params,cmdparams);
    } else if (cmd == "localhours") {
      DrawLocalHours(image,params,cmdparams);
    } else if (cmd == "altitudes") {
      DrawAltitudes(image,params,cmdparams);
    } else if (cmd == "analemma") {
      DrawAnalemma(image,params,cmdparams);
    } else if (cmd == "tropics") {
      DrawTropics(image,params,cmdparams);
    } else if (cmd == "dateline") {
      double p = strtod(cmdlist[i][1].c_str(),NULL);
      DrawDateline(image,params,cmdparams, p);
    } else if (cmd == "datetime") {
      double p = strtod(cmdlist[i][1].c_str(),NULL);
      DrawDatetime(image,params,cmdparams, p);
    } else if (cmd == "dlat") {
      double p = strtod(cmdlist[i][1].c_str(),NULL);
      cmdparams.dlat = radians(p);
    } else if (cmd == "dlong") {
      double p = strtod(cmdlist[i][1].c_str(),NULL);
      cmdparams.dlon = radians(p);
    } else if (cmd == "color") {
      cmdparams.color = Rgb(cmdlist[i][1].c_str());
    } else {
      fprintf(stderr,"Unknown command: %s\n", cmd.c_str());
    }
  }
}

void Transform::DrawLocalHours (Image& image, 
				const TransformParams& params,
				const CmdParams &cmdparams) const
{
  LatPlotter latplotter(image, params, *this);
  int nx = 24;
  for (int i = 0; i < nx; i++) {
    double lambda = (i - nx/2) * 2 * pi / nx + cmdparams.dlon;
    latplotter.lambda = lambda;
    // Omit the last section of the lines of latitude.
    image.PlotLine(-inclination, +inclination, latplotter, cmdparams.color, 4);
  }
}

void Transform::DrawGrid (Image& image, 
			  const TransformParams& params,
			  const CmdParams &cmdparams) const
{
  int gridx = cmdparams.gridx;
  int gridy = cmdparams.gridy;
  int gridoff = cmdparams.gridoff;
  double loff = 2*pi*gridoff/360.0;
  //loff = 0.1;
  fprintf(stderr,"%g\n",loff);
  int nx = 360/gridx;
  int ny = 180/gridy;
  for (int i = 0; i < nx; i++) {
    LatPlotter latplotter(image, params, *this);
    double lambda = (i - nx/2) * 2 * pi / nx + loff;
    latplotter.lambda = lambda;
    // Omit the last section of the lines of latitude.
    //image.PlotLine(-pi/2+radians(gridy), pi/2-radians(gridy), latplotter, cmdparams.color, 16);
    image.PlotLine(-pi/2, pi/2, latplotter, cmdparams.color, 16);
  }
  LongPlotter longplotter(image, params, *this);
  for (int i = 0; i <= ny; i++) {
    longplotter.phi = (i - ny/2) * pi / ny;
    image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
  }
}

void Transform::DrawAnalemma (Image& image, 
                              const TransformParams& params,
			      const CmdParams &cmdparams) const
{
  int gridx = cmdparams.gridx;
  AnalemmaPlotter analemmaplotter(image, params, *this);
  int nx = 360/gridx;
  for (int i = 0; i < nx; i++) {
    double time = 2 * pi * (i - nx/2) / nx;
    analemmaplotter.time = time;
    image.PlotLine(0, twopi, analemmaplotter, cmdparams.color,16);
  }
}

void Transform::DrawAltitudes(Image& image, 
			 const TransformParams& params,
			 const CmdParams &cmdparams) const
{
  FooPlotter fooplotter(image,params,*this);
  fooplotter.lambda = cmdparams.dlon;
  fooplotter.phi = cmdparams.dlat;
  for (int i = 10; i <= 80; i += 10) {
    fooplotter.theta = radians(i);
    image.PlotLine(0,twopi,fooplotter,cmdparams.color,16);
  }
}

void Transform::DrawTemporaryHours(Image& image, 
				   const TransformParams& params,
				   const CmdParams &cmdparams) const
{
  TempPlotter tempplotter(image, params, *this);
  tempplotter.phi = cmdparams.dlat;
  tempplotter.lambda = cmdparams.dlon;
  for (int i = 6; i <= 18; i++) {
    tempplotter.time  = i;
    image.PlotLine(-inclination,+inclination,tempplotter,cmdparams.color,4);
  }
}

void Transform::DrawTropics(Image& image, 
			    const TransformParams& params,
			    const CmdParams &cmdparams) const
{
  LongPlotter longplotter(image, params, *this);
  longplotter.phi = inclination;
  image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
  longplotter.phi = -inclination;
  image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
  longplotter.phi = 0.5*pi-inclination;
  image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
  longplotter.phi = inclination-0.5*pi;
  image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
}

void Transform::DrawDateline (Image& image, 
			      const TransformParams& params,
			      const CmdParams &cmdparams,
			      double day) const
{
  LongPlotter longplotter(image, params, *this);
  // Spring equinox is date 0, and is day 80 of a normal year.
  day = day - 80;
  while (day < 0) day += 365;
  while (day >= 365) day -= 365;
  double date = 2 * pi * day / 365;
  longplotter.phi = sundec(date);
  image.PlotLine(-pi, pi, longplotter, cmdparams.color,16);
}

void Transform::DrawDatetime (Image& image, 
			      const TransformParams& params,
			      const CmdParams &cmdparams,
			      double day) const
{
  LongPlotter longplotter(image, params, *this);
  // Spring equinox is day 0
  day = day - 80;
  while (day < 0) day += 365;
  while (day >= 365) day -= 365;
  double time = fmod(day, 1.0)-0.5; // Time relative to noon
  double date = 2 * pi * day / 365;
  double phi = sundec(date);
  double eot = equationoftime(date) / (24 * 60 * 60);
  double apparenttime = time + eot; // AT = MT + EOT
  double lambda = -(2 * pi * apparenttime);
  double x,y;
  MapXY(image, params, phi, lambda, x, y); 
  //fprintf(stderr,"%g %g %g %g %g %g\n", day, time,date,eot,phi,lambda);
  image.PlotPoint(x,y,1,cmdparams.color);
}
