// $Revision: 1.1 $
// matrix.cpp
// (C) 2004 by Matthew Arcus

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include "matrix.h"

Matrix3::Matrix3()
{ 
  Identity();
}

Matrix3::Matrix3(RotationAxis r, double theta)
{ 
  switch (r) {
  case XAxis: 
    {
      RotX(theta);
      break;
    }
  case YAxis: 
    {
      RotY(theta);
      break;
    }
  case ZAxis: 
    {
      RotZ(theta);
      break;
    }
  }
}

Matrix3::Matrix3(const Matrix3& m) 
{ 
  Copy(m);
}

void Matrix3::Copy(const Matrix3& m)
{
  for (int i = 0; i < 9; i++) {
      Set(i, m.Get(i));
  }
  identity = m.identity;
}

void Matrix3::Identity () {
  int i;
  for (i = 0; i < 9; i++) {
    Set(i,0.0);
  }
  for (i = 0; i < 9; i+=4){
    Set(i,1.0);
  }
  identity = true;
}

Matrix3 Matrix3::operator* (const Matrix3& m2) const
{
  Matrix3 result;
  result.Mult(*this, m2);
  return result;
}

// Multiply m1 by m2 and put the result in this
void Matrix3::Mult (const Matrix3& m1, const Matrix3& m2) 
{ 
  if (m1.identity) {
    if (m2.identity) {
      Identity();
    } else {
      Copy(m2);
    }
  } else if (m2.identity) {
    Copy(m1);
  } else {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	double x = 0.0;
	x += m1.Get(i,0) * m2.Get(0,j);
	x += m1.Get(i,1) * m2.Get(1,j);
	x += m1.Get(i,2) * m2.Get(2,j);
	Set(i,j,x);
      }
    }
    identity = false; // We could check the indices at this point
  }
}

// Set m to be a rotation of theta about the x-axis
void Matrix3::RotX (double theta) {
  if (theta == 0.0) {
    Identity();
  } else {
    Set(0, 1.0); Set(1, 0.0);        Set(2, 0.0);
    Set(3, 0.0); Set(4, cos(theta)); Set(5, -sin(theta));
    Set(6, 0.0); Set(7, sin(theta)); Set(8,  cos(theta));
    identity = false;
  }
}

// Set m to be a rotation of theta about the y-axis
void Matrix3::RotY (double theta) {
  if (theta == 0.0) {
    Identity();
  } else {
    Set(0, cos(theta)); Set(1, 0.0); Set(2, -sin(theta));
    Set(3, 0.0);        Set(4, 1.0); Set(5, 0.0);
    Set(6, sin(theta)); Set(7, 0.0); Set(8, cos(theta));
    identity = false;
  }
}

// Set m to be a rotation of theta about the z-axis
void Matrix3::RotZ (double theta) {
  if (theta == 0.0) {
    Identity();
  } else {
    Set(0, cos(theta)); Set(1, -sin(theta)); Set(2, 0.0);
    Set(3, sin(theta)); Set(4,  cos(theta)); Set(5, 0.0);
    Set(6, 0.0);        Set(7, 0.0);         Set(8, 1.0);
    identity = false;
  }
}

void Matrix3::Apply (double& x, double& y, double& z) const
{
  if (!identity) {
    double x1 = Get(0) * x + Get(1) * y + Get(2) * z;
    double y1 = Get(3) * x + Get(4) * y + Get(5) * z;
    double z1 = Get(6) * x + Get(7) * y + Get(8) * z;
    x = x1; y = y1; z = z1;
  }
}

void Matrix3::ApplyLatLong (double& phi, double& lambda) const
{
  if (!identity) {
    double x = cos(lambda) * cos(phi);
    double y = sin(lambda) * cos(phi);
    double z = sin(phi);
    Apply(x,y,z);
    phi = asin(z);
    lambda = atan2(y,x);
  }
}
