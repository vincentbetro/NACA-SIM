#include <stdio.h>
#include <math.h>
#include "Point.h"

#ifndef Vector_h
#define Vector_h
class Vector
{
 double vec[2];
 public:
  inline Vector(double dx=0.0, double dy=0.0)
  { vec[0] = dx; vec[1] = dy; }
  inline Vector(const Point &from, const Point &to);
  void print(FILE *outf)
  { fprintf(outf,"\nVector (x,y)= (%.12g,%.12g)",vec[0],vec[1]); }
  inline Vector operator + ( const Vector &v);
  inline Vector operator - ( const Vector &v);
  inline Vector &operator += ( const Vector &v);
  inline Vector &operator -= ( const Vector &v);
  inline Vector &operator += ( const double &s);
  inline Vector &operator -= ( const double &s);
  inline Vector &operator *= ( const double &s);
  inline Vector &operator /= ( const double &s);
  inline double operator * ( const Vector &v)
  { return double(vec[0]*v.vec[0]+vec[1]*v.vec[1]); }
  inline double operator % ( const Vector &v);
  inline Vector operator * ( const double &s);
  inline Vector operator / ( const double &s);
  inline double operator () (int i) const;
  inline double &operator () (int i);
  inline double &operator [] (int i);
  inline double magnitude();
  inline void normalize();
  ~Vector(){};

};

inline Vector::Vector(const Point &from, const Point &to)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i] = to(i) - from(i);
          }
}
inline double &Vector::operator () (int i)
{ return vec[i]; }
inline double Vector::operator () (int i) const
{ return vec[i]; }
inline double &Vector::operator [] (int i)
{ return vec[i]; }
inline double Vector::operator % ( const Vector &v)
{
  return (vec[0]*v.vec[1]-vec[1]*v.vec[0]);
}
inline Vector Vector::operator +( const Vector &v)
{ return Vector(vec[0]+v.vec[0], vec[1]+v.vec[1]); }
inline Vector Vector::operator -( const Vector &v)
{ return Vector(vec[0]-v.vec[0], vec[1]-v.vec[1]); }
inline Vector Vector::operator *( const double &s)
{ return Vector(s*vec[0], s*vec[1]); }
inline Vector Vector::operator /( const double &s)
{ return Vector(vec[0]/s, vec[1]/s); }
inline double Vector::magnitude()
{ return double(sqrt(vec[0]*vec[0] + vec[1]*vec[1])); }
inline Vector &Vector::operator += (const Vector &v)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]+=v.vec[i];
          }
  return *this;
}
inline Vector &Vector::operator -= (const Vector &v)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]-=v.vec[i];
          }
  return *this;
}
inline Vector &Vector::operator += (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]+=s;
          }
  return *this;
}
inline Vector &Vector::operator -= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]-=s;
          }
  return *this;
}
inline Vector &Vector::operator *= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]*=s;
          }
  return *this;
}
inline Vector &Vector::operator /= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]/=s;
          }
  return *this;
}

// we do not need to worry if mag less than machine eps
// if so, vector is 0 anyways, and if we divide nearly 0 by zero,
// we get a huge number and not a normalized vector like we want 
inline void Vector::normalize()
{
  double mag = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
  if (mag > 1.0e-20) { vec[0] = vec[0] / mag; vec[1] = vec[1] / mag; }
}

#endif
