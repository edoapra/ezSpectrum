#ifndef _molstate_h_
#define _molstate_h_

/*! \file molstate.h
\brief Molecular state: Stores Geometry, Normal Modes & Frequencies. Also reads input data from the XML file (vm'06) 
\ingroup MOLECULAR_PROP
*/

#include "genincludes.h"
#include "simple_xml_parser.h"
#include "atom.h"
#include "normalmode.h"
#include <vector>
#include "vector3d.h"
#include "kmatrix.h"
#include <algorithm>


class MolState
{
  //! N Atoms 
  std::vector<Atom> atoms;
  //! M Normal modes and frequencies (3N-5 or 3N-6) Stored in mass unweighted format in Angstoms (i.e as in ACESII)
  std::vector<NormalMode> normModes;
  //! Normal modes order (relative to the input file's order)
  std::vector<int> normModesOrder;
  //! may be removed later
  bool ifLinear;
  //! IP
  double energy;
  //! Input file
  simpleXMLparser xmlF; 
  //! "AtomName<->Mass" table
  simpleXMLparser amuF;  

  //!move this to functions...
  bool ifLetterOrNumber(char Ch);

  Vector3D centerOfMass;
  KMatrix momentOfInertiaTensor;

  //! reduced masses
  KMatrix reduced_masses;
  
  //!if geometry transformation was performed manually
  bool if_aligned_manually;
  bool if_nm_reordered_manually;

 public:

  MolState (const char* xmlFileName, const char* amuFileName);
  MolState (const MolState& other);
  //  ~MolState();

  MolState& operator=(const MolState& other);
  //! Read the state data from file
  bool Read(); 
  void Print(); 
  void printGeometry(); 
  void printNormalModes(); 
  // returns true if the overlap matrix is diagonal; makes a list of normal modes which form a non-diagonal minor.
  bool getNormalModeOverlapWithOtherState(MolState& other, KMatrix& overlap, std::vector<int>& normal_modes_list);

  //--- interface ---------------------------------------------------

  //! Returns Atom #i (this atoms can be different/isotops from the normal mode section)
  Atom& getAtom(int i){ return atoms[i]; } 
  //! Returns NormalMode #i
  NormalMode& getNormMode(int i) { return normModes[i]; } 
  //! Returns NormalMode order index
  int getNormModeIndex(int i) { return normModesOrder[i]; } 
  //! Linear?
  bool IfLinear () const { return ifLinear; }
  //! Energy (IP)
  double Energy () const { return energy; }
  //! Returns number of atoms
  int NAtoms() const { return atoms.size(); }
  //! Returns number of normal modes
  int NNormModes() const { return normModes.size(); }

  //--- alignment ------------- ------------------------------------
  //! align each state: center of mass in the coordinates origin, moment of ineretia principal axes along the coordinate axes:
  void align();
  //! align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  void align(MolState& other);
  //! check that states are similar (same number of atoms and same "liniarity", i.e. same number of normal modes); 
  bool ifSimilar(MolState& other);

  //--- Geometry transformations ------------------------------------
  // Center eof mass vector:
  Vector3D& getCenterOfMass();
  // Moment of inertia tensor:
  KMatrix& getMomentOfInertiaTensor();
  // Shifts coordinate's origin by "vector" 
  void shiftCoordinates(Vector3D& vector);
  // matrix multiolication coordinates*matrix_3x3; "coordinates" matrix is of 3 columns: x,y,z. matrix_3x3 has eigen vectors of the transformation in rows.
  void transformCoordinates(const KMatrix& matrix_3x3);
  void applyCoordinateThreshold(const double threshold);
  double getGeomDifference(MolState& other);
  // three rotations around x, y, z by PI/2:
  void rotateX_90deg();
  void rotateY_90deg();
  void rotateZ_90deg();
  // arbitrary rotation around x ,y, z axes by angles alpha_x, alpha_y,alpha_z respectevly:
  void rotate(const double alpha_x,const double alpha_y,const double alpha_z);
  bool ifAlignedManually();
  bool ifNMReorderedManually();

  //-- XML parser ---------------------------------------------------
  //! 2DO (should be from the parent XML_simple_parser)
  MolState& node(const char* elementName, int elementNumber) {  xmlF.node(elementName,elementNumber); return *this; };
  //! 2DO (should be from the parent XML_simple_parser)
  MolState& node(const char* elementName) {  xmlF.node(elementName,1); return *this; };
  //! 2DO (should be from the parent XML_simple_parser)
  MolState& reset() {  xmlF.reset(); return *this; };
  //! 2DO (should be from the parent XML_simple_parser)
  bool Check() {return xmlF.Check(); };
  // ----------------------------------------------------------------

};


#endif

/*!
2DO:
(should be from the parent XML_simple_parser)
// after read state (2DO check if was read) get
//DONE if number_of_atoms is greater than the actual number in the "text" ...
//hek -- check the parcer -- if specified file does not exist -- do something 
*/
