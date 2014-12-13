/*! \file atom.C
\brief Stores x,y,z mass, and atom name 
\ingroup DATA_CLASSES
*/

#include "atom.h"

double Atom::getR()
{
  double R=0;
  for (int i=0; i<CARTDIM; i++)
    R+=coord[i]*coord[i];
  return sqrt(R);
}

double Atom::getMomentumProj(const int axis)
{
  double tmpDbl;
  tmpDbl=coord[axis]*mass;
  return tmpDbl;
}


void Atom::shiftCoordinates(Vector3D& vector)
{
  for (int i=0; i<CARTDIM; i++)
      coord[i]-=vector[i];
}

void Atom::transformCoordinates(const KMatrix& matrix_3x3)
{
  double coord_old[CARTDIM];

  for (int i=0; i<CARTDIM; i++)
    {
      coord_old[i]=coord[i];
      coord[i]=0;
    }

  //vector*matrix multiplication:
  for (int i=0; i<CARTDIM; i++) // columns #, i.e. i=0 is x, i=1 is y, i=2 is z;
    for (int j=0; j<CARTDIM; j++) // rows #, i.e. eigen vectors, which are in rows;
      coord[i]+=coord_old[j]*matrix_3x3.Elem2(i,j);
}

void Atom::applyCoordinateThreshold(const double threshold)
{
  for (int i=0; i<CARTDIM; i++)
      if (fabs(coord[i])< threshold)
	coord[i]=0.0;
}

void Atom::rotateX_90deg()
{
  // x'=x; y'=z; z'=-y;
  double coord_tmp;
  coord_tmp=coord[1];// =y
  coord[1]=coord[2];  // y'=z
  coord[2]=-coord_tmp; // z'=-y
}

void Atom::rotateY_90deg()
{
  // x'=-z; y'=y; z'=x;
  double coord_tmp;
  coord_tmp=coord[0];// =x
  coord[0]=-coord[2]; // x'=-z
  coord[2]=coord_tmp; // z'=x
}

void Atom::rotateZ_90deg()
{
  // x'=y; y'=-x; z'=z;
  double coord_tmp;
  coord_tmp=coord[0];// =x
  coord[0]=coord[1];  // x'=y
  coord[1]=-coord_tmp; // y'=-x
}

