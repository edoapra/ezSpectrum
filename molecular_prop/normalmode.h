#ifndef _normalmode_h_
#define _normalmode_h_

/*! \file normalmode.h
\brief Normal mode: frequency and 3N displacements (N-number of atoms)
\ingroup DATA_CLASSES
*/

//ACHTUNG! different constructor from the previous versions -- now one should supply
//         it with number of atoms, NOT nAtoms*CartDim.

#include "genincludes.h"
#include "kmatrix.h"
#include "vector3d.h"

class NormalMode
{
  // 1D structure AtomNumber*CARTDIM+CoordinateNumber:
  KMatrix displacement;
  double freq;

  int nAtoms;

  //!  Re-allocates memory to size NOfAtoms_x_CARTDIM. 
  //  void AllocForNAtoms(int NOfAtoms) 
  //    {
  //      displacement.Adjust(NOfAtoms*CARTDIM,1); 
  //      nAtoms=NOfAtoms;
  //    }
  
 public:
  NormalMode(int n, double fr) : displacement(n*CARTDIM,1), freq(fr), nAtoms(n) {};
  NormalMode(const NormalMode& other);
  NormalMode& operator=(const NormalMode& other);
  
  //! Returns frequency
  double& getFreq() { return freq; }
  // 3N displacements (x,y,z for N atoms) as a (3Nx1) KMatrix
  KMatrix& getDisplacement() { return displacement; }

  // Shifts coordinate's origin by "vector" 
  // IMPORTANT: there is no need to "shift" normal coordinates in general. (i.e. never use this function!)
  void shiftCoordinates(Vector3D& vector);
  // matrix multiolication coordinates*matrix_3x3; "coordinates" matrix is of 3 columns: x,y,z. matrix_3x3 has eigen vectors of the transformation in rows.
  void transformCoordinates(const KMatrix& matrix_3x3);
  void applyCoordinateThreshold(const double threshold);
  // three rotations around x, y, z by PI/2:
  void rotateX_90deg();
  void rotateY_90deg();
  void rotateZ_90deg();
 
};
#endif


