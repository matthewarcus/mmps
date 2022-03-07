// $Revision: 1.1 $
// matrix.h
// (c) 2004-2022 Matthew Arcus

// MIT License

#if !defined MATRIX_H
#define MATRIX_H

enum RotationAxis {
  XAxis,
  YAxis,
  ZAxis,
};

class Matrix3 {
public:
  Matrix3();
  Matrix3(RotationAxis r, double theta);
  Matrix3(const Matrix3& m);
  void Apply (double& x, double& y, double& z) const;
  void ApplyLatLong (double& phi, double& lambda) const;
  Matrix3 operator *(const Matrix3& m2) const;
  inline double Get(int i, int j) const { return data[i + j + j + j]; };
  void Identity ();
  bool IsIdentity() const { return identity; };
private:
  void Copy(const Matrix3& m);
  void Mult (const Matrix3& m1, const Matrix3& m2);
  void RotX (double theta);
  void RotY (double theta);
  void RotZ (double theta);
private:
  inline double Get(int i) const { return data[i]; };
  inline void Set(int i, int j, double x) { data[i + j + j + j] = x; };
  inline void Set(int i, double x) { data[i] = x; };
  double data[9];
  bool identity;
};
#endif
