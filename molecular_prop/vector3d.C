/*! \file vector3d.C
\brief Stores x,y,z of a 3D vector 
\ingroup DATA_CLASSES
*/

#include "vector3d.h"

Vector3D& Vector3D::print(const char* nameTag)
{
  std::cout << nameTag;
  std::cout << "("<< coord[0] << ","  << coord[1] << ","  << coord[2] << ")\n";
  return *this;
}

double Vector3D::getNorm()
{
  return sqrt(coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2]);
}

Vector3D& Vector3D::operator=(const Vector3D& other)
{
  if (this != &other)
    {
      for (int i=0; i<CARTDIM; i++) 
	coord[i]=other.coord[i];
    }
  return *this;
}

Vector3D& Vector3D::operator+=(Vector3D& other)
{
  for (int i=0; i<3; i++)
    coord[i]+=other.getCoord(0);
  return *this;
}

Vector3D& Vector3D::operator-=(Vector3D& other)
{
  for (int i=0; i<3; i++)
    coord[i]-=other.getCoord(0);
  return *this;
}

Vector3D& Vector3D::operator*=(const double dbl)
{
  for (int i=0; i<3; i++)
    coord[i]*=dbl;
  return *this;
}

void Vector3D::reset()
{
  for (int i=0; i<3; i++)
    coord[i]=0;
}

Vector3D& Vector3D::applyThreshold(const double threshold)
{
  for (int i=0; i<3; i++)
    if ( fabs(coord[i]) < threshold)
      coord[i]=0.0;
  return *this;
}

