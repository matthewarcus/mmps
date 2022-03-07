// $Revision: 1.5 $
// transform.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined TRANSFORM_H
#define TRANSFORM_H

#include "constants.h"
#include "matrix.h"
#include "image.h"

class TransformParams {
public:
  TransformParams();
  double tilt;
  double turn;
  double rotate;
  double lat;
  double lon;
  double scale;
  double radius;
  double a;
  double aw;
  double x;
  double y;
  double z;
  double ox;
  double oy;
  double oz;
  bool   sun;
  double p;
  double time;
  double date;
  bool   noback;
  double conic;
  double conicr;
  double xoff;
  double yoff;
  Rgb background;

  // grid
  int gridx;
  int gridy;
  Rgb gridcolor;
  double gridoff;

  // dial
  bool analemma;
  Rgb analemmacolor;
  bool dial;
  Rgb dialcolor;
  double dlat;
  double dlon;
};
  
class Transform {
 public:
  void TransformImage(const Image& inImage, Image& outImage,
		      const TransformParams& params);
  void TransformImageInv(const Image& inImage, Image& outImage,
			 const TransformParams& params);
  bool MapXY(const Image& outImage, const TransformParams& params,
	     double phi, double lambda, double& x, double& y) const;
  void AddLines (Image& image, const TransformParams& params) const;
  void DoCommands (Image& image, 
		   const TransformParams& params,
		   const std::vector<std::vector<std::string> > cmdlist) const;
  virtual void SetY(double y) { (void(y)); };
  static Transform* GetTransform(const char *type);
  virtual void Init(const TransformParams& params);
  virtual void AdjustSize (int& w, int& h, TransformParams& params) { 
    (void(w)); (void(h)); (void(params)); 
  };
  virtual double BasicScale(int width, int height) const = 0;
  virtual bool Project(const TransformParams& params,
		       double x, double y, 
		       double& x1, double& y1, double& z1,
		       double& phi, double& lambda) const  = 0;
  virtual bool ProjectInv(const TransformParams& params,
			  double phi, double lambda, 
			  double& x, double& y) const = 0;
  virtual bool CheckPoint(const TransformParams &params, double phi, double lambda) const;
  virtual ~Transform() {};
 private:
  struct CmdParams {
    CmdParams() : color(),gridx(15),gridy(10),gridoff(0), dlat(0), dlon(0) {}
    Rgb color;
    int gridx;
    int gridy;
    double gridoff;
    double dlat;
    double dlon;
  };
  void DrawGrid (Image& image, 
		 const TransformParams& params,
		 const CmdParams &cmdparams) const;
  void DrawLocalHours (Image& image, 
		       const TransformParams& params,
		       const CmdParams &cmdparams) const;
  void DrawTemporaryHours (Image& image, 
			   const TransformParams& params,
			   const CmdParams &cmdparams) const;
  void DrawAltitudes (Image& image, 
		      const TransformParams& params,
		      const CmdParams &cmdparams) const;
  void DrawAnalemma (Image& image, 
		     const TransformParams& params,
		     const CmdParams &cmdparams) const;
  void DrawTropics (Image& image, 
		    const TransformParams& params,
		    const CmdParams &cmdparams) const;
  void DrawDateline (Image& image, 
		     const TransformParams& params,
		     const CmdParams &cmdparams,
		     double day) const;
  void DrawDatetime (Image& image, 
		     const TransformParams& params,
		     const CmdParams &cmdparams,
		     double daytime) const;
  void SetData(const Image& inImage, Image& outImage,
	       const TransformParams& params,
	       double x, double y, double z,
	       double phi, double lambda,
	       unsigned int outIndex);
  bool SetDataInv(const Image& inImage, Image& outImage,
		  const TransformParams& params,
		  double ox, double oy, // Coordinates in output image
		  double x, double y,   // scaled coordinates in input image
		  unsigned int outIndex);
 protected:
  Matrix3 m;
  Matrix3 minv;
  double sunx;
  double suny;
  double sunz;
};

class PolarTransform : public Transform
{
public:
  void Init(const TransformParams& params);
  bool Project(const TransformParams& params,
	       double x, double y, 
	       double& x1, double& y1, double& z1,
	       double& phi, double& lambda) const;
  bool ProjectInv (const TransformParams& params,
		   double phi, double lambda, double& x, double& y) const;
  double BasicScale(int width, int height) const { (void(width)); return 2.0 / height; };
private:
  virtual void GetPhi(double r, double& phi) const = 0;
  virtual bool GetR(double phi, double& r) const = 0;
protected:
  double k;
};
  
class CylindricalTransform : public Transform
{
  // x is g(phi)*lambda
  // y is f(phi)
  // The drawing area is 2pi units across. The scale factor therefore
  // is w/2pi. The north pole is at w/2pi * f(pi/2)
public:
  void AdjustSize (int& w, int& h, TransformParams& params);
  bool Project(const TransformParams& params,
	       double x, double y, 
	       double& x1, double& y1, double& z1,
	       double& phi, double& lambda) const;
  bool ProjectInv(const TransformParams& params,
		  double phi, double lambda, double& x, double& y) const;
  double BasicScale(int width, int height) const { (void(height)); return 2.0 * pi / width; };
  void SetY(double y) { phi0 = GetLat(y); };
private:
  virtual bool GetXY(double phi, double lambda, double& x, double& y) const = 0;
  virtual double GetLat(double y) = 0;
  virtual double GetLong(double x) const = 0;
  virtual double GetMaxHeight(const TransformParams& params) { (void(params)); return 0.0; };
 private:
  double phi0;
};

class SimpleTransform: public Transform
{
 public:
  bool Project(const TransformParams& params,
	       double x0, double y0, 
	       double& x, double& y, double& z,
	       double& phi, double& lambda) const;
  bool ProjectInv(const TransformParams& params,
		  double phi, double lambda, double& x, double& y) const;
  virtual bool ProjectSimple(double x, double y, 
			     double& phi, double& lambda) const = 0;
  virtual bool ProjectInvSimple(double phi, double lambda,
				double &x, double &y) const = 0;
};

// Now the specific transforms

class LatLong : public CylindricalTransform {
public:
  double GetMaxHeight(const TransformParams& params);
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
private:
  double k;
};

// f(phi) = sin(phi)/k, g(phi) = 1, where k = sqr(cos(p))
class EqualArea : public CylindricalTransform {
private:
  double GetMaxHeight(const TransformParams& params);
  void Init(const TransformParams& params);
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
private:
  double k;
};

class Mercator : public CylindricalTransform {
private:
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
};

class Cylindrical : public CylindricalTransform {
private:
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
};

class Sinusoidal : public CylindricalTransform {
  double GetMaxHeight (const TransformParams& params);
private:
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
private:
  double m;
};

class Sinusoidal2 : public CylindricalTransform {
private:
  double GetMaxHeight(const TransformParams& params);
  void Init(const TransformParams& params);
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
  static void CalcParams(double a, double phi, double& k, double& b);
private:
  double a;
  double b;
  double k;
  double m;
};

class Mollweide : public CylindricalTransform {
private:
  double GetMaxHeight(const TransformParams& params);
  void Init(const TransformParams& params);
  double GetLat(double y);
  double GetLong(double x) const;
  bool GetXY(double phi, double lambda, double& x, double& y) const;
  static void CalcParams(double phi, double& a);
private:
  double a;
  double m;
};

class Azimuthal : public PolarTransform
{
 public:
  void AdjustSize (int& w, int& h, TransformParams& params);
 private:
  void GetPhi(double r, double& phi) const;
  bool GetR(double phi, double& r) const;
};
  
class Rectilinear : public PolarTransform
{
private:
  void GetPhi(double r, double& phi) const;
  bool GetR(double phi, double& r) const;
};

// This is a conformal map.  
class Stereographic : public PolarTransform
{
private:
  void GetPhi(double r, double& phi) const;
  bool GetR(double phi, double& r) const;
};
  
class Gnomonic : public PolarTransform
{
private:
  void GetPhi(double r, double& phi) const;
  bool GetR(double phi, double& r) const;
};

class Bonne : public SimpleTransform
{
 public:
  void Init(const TransformParams& params);
  double BasicScale (int width, int height) const { (void(height)); return 2.0 * pi / width; };
  bool ProjectSimple(double x, double y, 
		     double& phi, double& lambda) const;
  bool ProjectInvSimple(double phi, double lambda,
			double &x, double &y) const;
 private:
  double p;
  double cotp;
};

class Hammer : public SimpleTransform
{
 public:
  void Init(const TransformParams& params);
  double BasicScale (int width, int height) const { (void(height)); return 2.0 / width; };
  void AdjustSize (int& w, int& h, TransformParams& params) { (void(params)); h = w/2; };
  bool ProjectSimple(double x, double y, double& phi, double& lambda) const;
  bool ProjectInvSimple(double phi, double lambda, double& x, double& y) const;
};

class Lambert : public SimpleTransform
{
 public:
  void Init(const TransformParams& params);
  double BasicScale (int width, int height) const { (void(height)); return 2.0 / width; };
  void AdjustSize (int& w, int& h, TransformParams& params) { (void(params)); h = w/2; };
  bool ProjectSimple(double x, double y, double& phi, double& lambda) const;
  bool ProjectInvSimple(double phi, double lambda, double& x, double& y) const;
};

class Perspective : public Transform
{
 public:
  void Init(const TransformParams& params);
  bool Project(const TransformParams& params,
	       double x0, double y0, 
	       double& x, double& y, double& z,
	       double& phi, double& lambda) const;
  double BasicScale (int width, int height) const { (void(width)); return 2.0 / height; };
  bool ProjectInv (const TransformParams& params,
		   double phi, double lambda, double& x, double& y) const;
 private:
  double scalefactor;
  double iscalefactor;
  // Viewpoint position
  double vx;
  double vy;
  double vz;
  // Rotated viewpoint position
  double rx;
  double ry;
  double rz;
  // Cached oblateness factors
  double a; double b; double c;
  // And their squares
  double a2; double b2; double c2;
  // And their inverse squares
  double ia2; double ib2; double ic2;
};

#endif
