#ifndef _vector3d_h_
#define _vector3d_h_

/*! \file vector3d.h
\brief Stores x,y,z of a 3D vector 
\ingroup DATA_CLASSES
*/

#include "genincludes.h"
#include "constants.h"

//! Vector3D: stores x,y,z.                 
class Vector3D
{
  double coord[CARTDIM];

 public:
  Vector3D(){coord[0]=coord[1]=coord[2]=0;}
  Vector3D(const double x,const double y,const double z){coord[0]=x; coord[1]=y; coord[2]=z;}

  //  double& Coord(int i) {return coord[i];}
  double& getCoord(const int i) {return coord[i];}
  Vector3D& print(const char* nameTag);

  double getNorm();

  // coord should be const..., as whel as 'other':
  Vector3D& operator=(const Vector3D& other);
  Vector3D& operator+=(Vector3D& other);
  Vector3D& operator-=(Vector3D& other);
  Vector3D& operator*=(const double dbl);
  double operator[](const int index) {return getCoord(index); };

  
 Vector3D& applyThreshold(const double threshold);

  void reset();
};

#endif

