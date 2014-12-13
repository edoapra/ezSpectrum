#include "molstate.h"

/*! \file molstate.C
\brief Molecular state: Stores Geometry, Normal Modes & Frequencies. Also reads input data from the XML file (vm'06) 
\ingroup DATA_CLASSES
*/

//------------------------------
MolState::MolState (const char* xmlFileName, const char* amuFileName)
{
  xmlF.assignFile(xmlFileName);
  amuF.assignFile(amuFileName);
  momentOfInertiaTensor.Adjust(CARTDIM,CARTDIM);
  if_aligned_manually=false;
  if_nm_reordered_manually=false;
  normModesOrder.clear();
}


//------------------------------
MolState::MolState (const MolState& other)
{
  atoms=other.atoms;
  normModes=other.normModes;
  ifLinear=other.ifLinear;
  energy=other.energy;
  xmlF=other.xmlF;
  amuF=other.amuF;
  centerOfMass = other.centerOfMass;
  momentOfInertiaTensor = other.momentOfInertiaTensor;
  if_aligned_manually = other.if_aligned_manually;
  if_nm_reordered_manually =other.if_nm_reordered_manually;
  normModesOrder = other.normModesOrder;
  reduced_masses=other.reduced_masses;
}

//------------------------------
MolState& MolState::operator=(const MolState& other)
{
  if (this != &other)
    {
      atoms=other.atoms;
      normModes=other.normModes;
      ifLinear=other.ifLinear;
      energy=other.energy;
      xmlF=other.xmlF;
      amuF=other.amuF;
      centerOfMass = other.centerOfMass;
      momentOfInertiaTensor = other.momentOfInertiaTensor;
      if_aligned_manually = other.if_aligned_manually;
      if_nm_reordered_manually =other.if_nm_reordered_manually;
      normModesOrder = other.normModesOrder;
      reduced_masses=other.reduced_masses;
    }
  return *this;
}


//------------------------------
void MolState::align()
{
  // Shift center of mass to the origin:
  shiftCoordinates( getCenterOfMass() );

  // Aligh moment of inertia principal axes:
  KMatrix MOI_tensor(3,3), MOI_eigenValues(3,1), MOI_eigenVectors(3,3);
  MOI_tensor = getMomentOfInertiaTensor();
  //MOI_tensor.Print("Moment of inertia tensor");
  MOI_eigenVectors=MOI_tensor.Herm_Diag(MOI_eigenValues, true); // true=>eigen vals in the descending order; eigen vectors in ROWS.
  // MOI_eigenVectors.Print("MOI eigen vectors");
  // MOI_eigenValues.Print("MOI eigen values");
  
  // compute the determinant of MOI matrix:
  double det_tmp =  MOI_eigenVectors.Determinant();
  // if determinant is = 1 than it is a proper rotation; 
  // if it is = -1 than it is rotation+reflection (switch from left handed to the right handed coord.system); swap x&y axes:
  if (det_tmp < 0)
    MOI_eigenVectors.SwapRows(0,1);

  transformCoordinates(MOI_eigenVectors); // rotate coordinates using matrix of the MOI eigen vectors 

  std::cout << "Coordinate axes were aligned with MOI principal axes; center mass was shifted to the origin.\n";
}


//------------------------------
void MolState::align(MolState& other)

{
  //align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  // by minimizing the "sum" of distances between the same atoms: SUM[i=1..Natoms](deltaRi)
  // DeltaRi: distance between i'th atoms of the initial and target states

  double min_diff=DBL_MAX; // minimum of the "sum"
  int min_x_rot, min_y_rot, min_z_rot; // rotations, that correspond to the minimum of the "sum"

  // x_rot, y_rot, z_rot = {0, pi/2, pi, 3pi/2} rotation angle around x, y, and z axes ;
  // rotation in 3D is not commutative, so first rotate around X, than Y, than Z:
  for (int z_rot=0;  z_rot<4; z_rot++)
    {
      for (int y_rot=0;  y_rot<4; y_rot++)
	{
	  for (int x_rot=0;  x_rot<4; x_rot++)
	    {
	      double diff=getGeomDifference(other);
	      if (diff<=min_diff)
		{
		  min_diff=diff;
		  min_x_rot=x_rot;
		  min_y_rot=y_rot;
		  min_z_rot=z_rot;
		}
	      rotateX_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
	    }
	  rotateY_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
	}
      rotateZ_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry.
    }// at the end the molecule rotated by 4*pi/2 around x and y, i.e. it is non-rotated
  
  // rotate to the angles with the min norm of the deltaRi differences (min_x_rot, min_y_rot, min_z_rot):
  
  rotate(min_x_rot*PI/2, min_y_rot*PI/2, min_z_rot*PI/2);

  std::cout << "Also rotated by " << min_x_rot <<"/2*pi, " 
	    << min_y_rot <<"/2*pi, and " << min_z_rot <<"/2*pi CCW around x, y, and z.\n";
  std::cout << "The norm of the geometry difference from the initial state is " << sqrt(min_diff)<<"\n";
}


//------------------------------
bool MolState::ifSimilar(MolState& other)
{
  bool return_bool=true;
  std::string error_str="";
  if ( NAtoms()!=other.NAtoms() )
    {
      std::cout << "\nError: different number of atoms in the initial and the target states\n\n";
      return_bool=false;	
    }
  else if ( NNormModes()!=other.NNormModes() )
    {
      std::cout << "\nError: different number of normal modes in the initial and the target state\n\n";
      return_bool=false;	
    }

  return return_bool;
}


//------------------------------
Vector3D& MolState::getCenterOfMass()
{
  centerOfMass.reset();
  double totalMass=0;
  for (int i=0; i<atoms.size();i++)
    {
      for (int axis=0; axis<3; axis++)
	centerOfMass.getCoord(axis)+=atoms[i].getMomentumProj(axis);
      totalMass+=atoms[i].getMass();
    }
  centerOfMass*=1/totalMass;
  return centerOfMass;
}

//------------------------------
KMatrix& MolState::getMomentOfInertiaTensor()
{
  //  I_ij=SUM_k[m_k*(r_k^2*delta_ij-r_ki*r_kj)]
  for (int i=0; i< CARTDIM; i++)
    for (int j=0; j<CARTDIM; j++)
      {
	momentOfInertiaTensor.Elem2(i,j)=0;
	for (int atom=0; atom<NAtoms(); atom++)
	  { 
	    momentOfInertiaTensor.Elem2(i,j)-=atoms[atom].getCoord(i)*atoms[atom].getCoord(j)*atoms[atom].getMass();

	    // add R^2, i.e. I_zz=m*(R^2-r_z^2)=m*(r_x^2+r_y^2)
	    if (i==j)
	      momentOfInertiaTensor.Elem2(i,j)+=atoms[atom].getR()*atoms[atom].getR()*atoms[atom].getMass();
	  }
      }
  return momentOfInertiaTensor;
}


//------------------------------
void MolState::shiftCoordinates(Vector3D& vector)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].shiftCoordinates(vector);
  // there is no need to shift normal coordinates!!
}

//------------------------------
void MolState::transformCoordinates(const KMatrix& matrix_3x3)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].transformCoordinates(matrix_3x3);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].transformCoordinates(matrix_3x3);
}

//------------------------------
void MolState::rotateX_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateX_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateX_90deg();
}

//------------------------------
void MolState::rotateY_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateY_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateY_90deg();
}

//------------------------------
void MolState::rotateZ_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateZ_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateZ_90deg();
}

//------------------------------
void MolState::rotate(const double alpha_x, const double alpha_y, const double alpha_z)
{
  // three rotation matrices (instead of making one matrix) arouns x, y, z axes
  KMatrix R(CARTDIM,CARTDIM), Rx(CARTDIM,CARTDIM), Ry(CARTDIM,CARTDIM), Rz(CARTDIM,CARTDIM);

  Rx.Set(0);
  Rx.Elem2(0,0)=1;
  Rx.Elem2(1,1)=cos(alpha_x);
  Rx.Elem2(2,2)=cos(alpha_x);
  Rx.Elem2(2,1)=-sin(alpha_x);
  Rx.Elem2(1,2)=sin(alpha_x);

  Ry.Set(0);
  Ry.Elem2(1,1)=1;
  Ry.Elem2(0,0)=cos(alpha_y);
  Ry.Elem2(2,2)=cos(alpha_y);
  Ry.Elem2(2,0)=sin(alpha_y);
  Ry.Elem2(0,2)=-sin(alpha_y);

  Rz.Set(0);
  Rz.Elem2(2,2)=1;
  Rz.Elem2(0,0)=cos(alpha_z);
  Rz.Elem2(1,1)=cos(alpha_z);
  Rz.Elem2(1,0)=-sin(alpha_z);
  Rz.Elem2(0,1)=sin(alpha_z);

  // overall rotation R=Rx*Ry*Rz
  R.SetDiagonal(1);
  R*=Rx;
  R*=Ry;
  R*=Rz;
  
  // and now rotates using matix R:
  transformCoordinates(R);
}
//------------------------------
bool MolState::ifAlignedManually()
{
  if (if_aligned_manually)
    return true;
  else
    return false;
}
//------------------------------
bool MolState::ifNMReorderedManually()
{
  if (if_nm_reordered_manually)
    return true;
  else
    return false;
}

//------------------------------
double MolState::getGeomDifference(MolState& other)
{
  double diff=0;
  for (int i=0; i<atoms.size(); i++)
    for (int j=0; j<CARTDIM; j++) 
      // diff+=DeltaR^2[=DetaX^2+DeltaY^2+DeltaZ^2]
      diff+= (getAtom(i).getCoord(j)-other.getAtom(i).getCoord(j)) * (getAtom(i).getCoord(j)-other.getAtom(i).getCoord(j));
  return diff;
}

//------------------------------
void MolState::applyCoordinateThreshold(const double threshold)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].applyCoordinateThreshold(threshold);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].applyCoordinateThreshold(threshold);
}


//------------------------------
bool MolState::Read()
{
  int i,j,k,l;
  char tmpCh;
  std::string tmp_str, tmp_atomName;
  std::istringstream tmp_iStr;
  Atom tmp_atom;

  //------------ Read IP (if provided) ----------------------------------

  if ( xmlF.CheckSubNode("ip") )
    {
      xmlF.node("ip");
      energy = xmlF.getDoubleValue();
      xmlF.stepBack();
    }
  //------------ Read the Geometry --------------------------------------
  
  xmlF.node("geometry");
  
  int tmp_nAtoms, tmp_nNormMds;
  tmp_nAtoms=xmlF.getIntValue("number_of_atoms");
 
  std::string units;
  units=xmlF.value("units");

  ifLinear= xmlF.getBoolValue("linear"); 
  if (ifLinear) 
    tmp_nNormMds = (3*tmp_nAtoms - 5);
  else
    tmp_nNormMds = (3*tmp_nAtoms - 6);

  tmp_iStr.str(xmlF.value("text")); 
  tmp_iStr.clear();
  atoms.clear();

  for (i=0; i<tmp_nAtoms; i++)
    {
      tmp_atomName="";
      tmpCh=' ';

      while ( !ifLetterOrNumber(tmpCh) and !(tmp_iStr.fail()) )
	tmp_iStr.get(tmpCh);  

      if (tmp_iStr.fail())
	{
	  std::cout << "\nError: less then "<< tmp_nAtoms <<" atoms found or the format error, \"text\"-s line nuumber "<<i+1<<"\n";
	  xmlF.exitOnFormatError(true);
	}
      
      do {  
	tmp_atomName+=tmpCh; 
	tmp_iStr.get( tmpCh );  
      } while( ifLetterOrNumber(tmpCh) );
      tmp_atom.Name() = tmp_atomName;

      for (k=0; k<CARTDIM; k++)
	{
	  tmp_iStr >> tmp_atom.Coord(k);
	  if (units=="au")
	    tmp_atom.Coord(k)*=AU2ANGSTROM;
	}

      if (tmp_iStr.fail())
	{
	  std::cout << "\nError: less then "<< tmp_nAtoms <<" atoms found or the format error in the atomic coordinates,\n"
		    << "       \"text\"-s line number "<<i+1<<"\n";
	  xmlF.exitOnFormatError(true);
	}
      
      atoms.push_back(tmp_atom);
    }


  
  //------------ Convert atomic names to masses --------------------------

   for (i=0; i<NAtoms(); i++)
    {
      tmp_iStr.str( amuF.reset().node("masses").node(  getAtom(i).Name().c_str()   ).value() ); 
      tmp_iStr.clear();
      tmp_iStr >> getAtom(i).Mass(); 
      amuF.exitOnFormatError(tmp_iStr.fail());
    }

  //------------ Read Normal Modes ---------------------------------------
  NormalMode tmp_normMode(NAtoms(),0); // one temp. norm mode

  normModes.clear();

  xmlF.stepBack();
  xmlF.node("normal_modes");

  tmp_iStr.str(xmlF.value("text")); 
  tmp_iStr.clear();

  for (i=0; i < tmp_nNormMds; i++)
    normModes.push_back( tmp_normMode );
 
  int nModesPerLine=3;  // three number of vib. modes per Line

  int nLines;
  nLines = tmp_nNormMds / nModesPerLine;
  if ( tmp_nNormMds % nModesPerLine != 0)
    nLines++;

  for (k = 0; k < nLines; k++)  // number of blocks with 3 norm.modes. ("lines")
    {
      int current_nModesPerLine = nModesPerLine;
      // for the last entree, nModesPerString may differ from 3
      if (nLines - 1 == k)
  	if ( tmp_nNormMds % nModesPerLine != 0 )
	  current_nModesPerLine = tmp_nNormMds % nModesPerLine;
 
      for (i=0; i < NAtoms(); i++)   
  	for (j=0; j < current_nModesPerLine; j++) 
  	  for (l=0; l < CARTDIM; l++)
	    {
	      tmp_iStr >> normModes[k*nModesPerLine+j].getDisplacement()[i*CARTDIM+l]; 
	      if (tmp_iStr.fail())
		{
		  std::cout << "\nError: less normal modes found than expected or the format error\n";
		  xmlF.exitOnFormatError(true);
		}
	    }
    }
   
  //------------ Read Frequencies ----------------------------------------
  xmlF.stepBack();
  xmlF.node("frequencies");

  tmp_iStr.str(xmlF.value("text")); 
  tmp_iStr.clear();
  
  for (i=0; i < tmp_nNormMds; i++)
    {
      tmp_iStr >> getNormMode(i).getFreq();
      if (tmp_iStr.fail())
	{
	  std::cout << "\nError: format error at frequency #"<<i<< " or less than " << tmp_nNormMds <<" frequencies found\n";
	  xmlF.exitOnFormatError(true);
	}
      if (getNormMode(i).getFreq()<=0)
	{
	  std::cout <<"\nError. The frequency ["<<getNormMode(i).getFreq() <<"] is negative\n";
	  xmlF.exitOnFormatError(true);
	}
    }

  // Now MolState is in a good shape, and some transformations can be performed

  //------------ 1. Un-mass-weight normal modes----------------------------

  xmlF.stepBack();
  xmlF.node("normal_modes");
  bool if_massweighted;
  if_massweighted=xmlF.getBoolValue("if_mass_weighted");

  //------------ 1. mass un-weight normal modes, if needed (QChem-->ACES format; )----------------
  // qchem if_massweighted="true"; aces if_massweighted="false";
  reduced_masses.Adjust(NNormModes(),1);
  reduced_masses.Set(1);
  if (if_massweighted)
    {

      // Read atomic names from "...->normal_modes->atoms"
      std::vector<Atom> normalModeAtoms;

      tmp_iStr.str(xmlF.value("atoms")); 
      tmp_iStr.clear();
      normalModeAtoms.clear();

      for (i=0; i<NAtoms(); i++)
	{
	  tmp_atomName="";
	  tmpCh=' ';

	  while ( !ifLetterOrNumber(tmpCh) and !(tmp_iStr.fail()) )
	    tmp_iStr.get(tmpCh);  

	  xmlF.exitOnFormatError(tmp_iStr.fail());
      
	  do {  
	    tmp_atomName+=tmpCh; 
	    tmp_iStr.get( tmpCh );  
	  } while( ifLetterOrNumber(tmpCh) );
	  tmp_atom.Name() = tmp_atomName;
	  normalModeAtoms.push_back(tmp_atom);
	}
      // Get masses for each atomic name:
      for (i=0; i<NAtoms(); i++)
	{
	  tmp_iStr.str( amuF.reset().node("masses").node(  normalModeAtoms[i].Name().c_str()   ).value() ); 
	  tmp_iStr.clear();
	  tmp_iStr >> normalModeAtoms[i].Mass(); 
	  amuF.exitOnFormatError(tmp_iStr.fail());
	}
      // Mass-un-weight normal modes:
      for (int nm=0; nm<NNormModes(); nm++)
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i) *= sqrt(normalModeAtoms[a].Mass());

      // normolize each normal mode (/sqrt(norm) which is also /sqrt(reduced mass)):
      // KMatrix reduced_masses(NNormModes(),1);
      for (int nm=0; nm<NNormModes(); nm++)
	{
	  reduced_masses[nm] = 0;
	  for (int a=0; a<NAtoms(); a++)
	    for (int i=0; i<CARTDIM; i++ )
	      reduced_masses[nm]+= getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i) * getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i);
	}
      // reduced_masses.Print("Reduced masses:");
      // Normalize:
      for (int nm=0; nm<NNormModes(); nm++)
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i)/=sqrt(reduced_masses[nm]);

    }
  xmlF.stepBack();

  //------------ 2. Align geometry if requested ----------------------------
  if_aligned_manually=false;
  double man_rot_x, man_rot_y, man_rot_z;
  Vector3D man_shift;
  if ( xmlF.CheckSubNode("manual_coordinates_transformation") )
     {
       xmlF.node("manual_coordinates_transformation");
       
       man_rot_x=xmlF.getDoubleValue("rotate_around_x"); 
       man_rot_y=xmlF.getDoubleValue("rotate_around_y"); 
       man_rot_z=xmlF.getDoubleValue("rotate_around_z"); 

       man_shift.getCoord(0)=-xmlF.getDoubleValue("shift_along_x"); 
       man_shift.getCoord(1)=-xmlF.getDoubleValue("shift_along_y"); 
       man_shift.getCoord(2)=-xmlF.getDoubleValue("shift_along_z"); 
       
       std::cout << "Molecular structure and normal modes of this electronic state\nwill be transformed as requested in the input.\n";

       shiftCoordinates(man_shift);
       rotate(man_rot_x*PI,man_rot_y*PI,man_rot_z*PI);
       applyCoordinateThreshold(COORDINATE_THRESHOLD);
       if_aligned_manually=true;

       xmlF.stepBack();
     }


  //------------ 3. Reorder normal modes if requested --------------------
   if ( xmlF.CheckSubNode("manual_normal_modes_reordering") )
    {
      xmlF.node("manual_normal_modes_reordering");
      if_nm_reordered_manually=false;
      normModesOrder.clear();

      std::cout << "New normal modes order was requested:\n" << xmlF.value("new_order") <<"\n";
      tmp_iStr.str(xmlF.value("new_order")); 
      int tmpInt;
      for (int nm=0; nm < NNormModes(); nm++)
	{
	  tmp_iStr >> tmpInt;
	  //input error check:
	  if (tmp_iStr.fail())
	    {
	      std::cout << "\nFormat error: non numeric symbol or less entries then the number of normal modes\n\n";
	      xmlF.exitOnFormatError(true);
	    }
	  if ( (tmpInt<0) or (tmpInt>=NNormModes()) )
	    {
	      std::cout << "\nError: normal mode number ["<< tmpInt<<"] is out of range [0.."<<NNormModes()-1<<"].\n\n";
	      xmlF.exitOnFormatError(true);
	    }
	  normModesOrder.push_back(tmpInt);
	}
      
      // check if there are duplicates in the list:
      std::vector<int> tmpIntVector, tmpIntVector2;
      tmpIntVector = normModesOrder;
      std::sort( tmpIntVector.begin(), tmpIntVector.end() );
      tmpIntVector2 = tmpIntVector;
      std::vector<int>::const_iterator intVec_iter;
      intVec_iter= unique( tmpIntVector.begin(), tmpIntVector.end() );
      if (intVec_iter != tmpIntVector.end())
	{
	  std::cout << "\nFormat error: there are non unique entries. Check the sorted list:\n";
	  for (std::vector<int>::const_iterator tmp_iter=tmpIntVector2.begin(); tmp_iter!=tmpIntVector2.end(); tmp_iter++)
	    std::cout << ' ' << *tmp_iter;
	  std::cout<<'\n';
	  xmlF.exitOnFormatError(true);
	}
      
      // backup normal modes
      std::vector<NormalMode> oldNormModes;
      oldNormModes.clear();
      for (int nm=0; nm < NNormModes(); nm++)
	{
	  tmp_normMode = getNormMode(nm);
	  oldNormModes.push_back( tmp_normMode );
	}
      
      // copy normal modes using new order
      for (int nm=0; nm < NNormModes(); nm++)
	getNormMode(nm) = oldNormModes[  normModesOrder[nm] ];
      
      std::cout << "Normal modes were reordered accordingly.\n";
      if_nm_reordered_manually=true;

      xmlF.stepBack();
    }
   else
     for (int nm=0; nm < NNormModes(); nm++)
       normModesOrder.push_back(nm);
    
  //------------ 4. Reorder atoms if requested --------------------

  Atom tmp_Atom; // one temp. norm mode

  if ( xmlF.CheckSubNode("manual_atoms_reordering") )
    {
      xmlF.node("manual_atoms_reordering");
      std::vector<int> atomsOrder;
      atomsOrder.clear();

      std::cout << "New order of atoms was requested:\n" << xmlF.value("new_order") <<"\n";
      tmp_iStr.str(xmlF.value("new_order")); 
      int tmpInt;
      for (int nm=0; nm < NAtoms(); nm++)
	{
	  tmp_iStr >> tmpInt;
	  //input error check:
	  if (tmp_iStr.fail())
	    {
	      std::cout << "\nFormat error: non numeric symbol or less entries then the number of atoms\n\n";
	      xmlF.exitOnFormatError(true);
	    }
	  if ( (tmpInt<0) or (tmpInt>=NAtoms()) )
	    {
	      std::cout << "\nError: stom number ["<< tmpInt<<"] is out of range [0.."<<NAtoms()-1<<"].\n";
	      xmlF.exitOnFormatError(true);
	    }
	  atomsOrder.push_back(tmpInt);
	}
      
      // check if there are duplicates in the list:
      std::vector<int> tmpIntVector, tmpIntVector2;
      tmpIntVector = atomsOrder;
      std::sort( tmpIntVector.begin(), tmpIntVector.end() );
      tmpIntVector2 = tmpIntVector;
      std::vector<int>::const_iterator intVec_iter;
      intVec_iter= unique( tmpIntVector.begin(), tmpIntVector.end() );
      if (intVec_iter != tmpIntVector.end())
	{
	  std::cout << "\nFormat error: there are non unique entries. Check the sorted list:\n";
	  for (std::vector<int>::const_iterator tmp_iter=tmpIntVector2.begin(); tmp_iter!=tmpIntVector2.end(); tmp_iter++)
	    std::cout << ' ' << *tmp_iter;
	  std::cout<<'\n';
	  xmlF.exitOnFormatError(true);
	}
      
      // backup molecular geometry and normal modes
      std::vector<Atom> oldAtoms;
      oldAtoms.clear();
      for (int a=0; a < NAtoms(); a++)
	{
	  tmp_Atom = getAtom(a);
	  oldAtoms.push_back( tmp_Atom );
	}
      std::vector<NormalMode> oldNormModes;
      oldNormModes.clear();
      for (int nm=0; nm < NNormModes(); nm++)
	{
	  tmp_normMode = getNormMode(nm);
	  oldNormModes.push_back( tmp_normMode );
	}
      
      // copy molecular geometry using the new order of atoms
      for (int a=0; a < NAtoms(); a++)
	getAtom(a) = oldAtoms[  atomsOrder[a] ];
      
      // copy normal modes using the new order of atoms
      for (int nm=0; nm < NNormModes(); nm++)
	for (int a=0; a < NAtoms(); a++) 
	  for (int k=0; k < CARTDIM; k++)
	    getNormMode(nm).getDisplacement()[a*CARTDIM+k]=oldNormModes[nm].getDisplacement()[  atomsOrder[a]*CARTDIM + k  ];
	  
      std::cout << "Atoms were reordered accordingly.\n";

      xmlF.stepBack();
    }
  //----------------------------------------------------------------------
  // print normal modes in original format (i.e. mass weighted or not -- as in the input):
  /*
  std::cout << "Normal modes:\n";
  for (int nm=0; nm < NNormModes(); nm++)
    {
      std::cout <<"Normal mode #"<< nm<< ":\n";
      for (int a=0; a < NAtoms(); a++) 
	{
	  for (int k=0; k < CARTDIM; k++)
	    std::cout << getNormMode(nm).getDisplacement()[a*CARTDIM+k]*sqrt( reduced_masses[nm]/getAtom(a).Mass() )<<' ';
	  std::cout << '\n';
	}
      std::cout << '\n';
    }
  */
 
  //----------------------------------------------------------------------
  return true;
}

bool MolState::getNormalModeOverlapWithOtherState(MolState& other, KMatrix& overlap, std::vector<int>& normal_modes_list)
{
  overlap; //#NM x #NM
  overlap.Adjust(NNormModes(),NNormModes());
  overlap.Set(0.0);

  //normalization of initial and target normal modes (stored as un-mass-weighted)
  double norm_ini, norm_targ;
  // norm of the L
  for (int nm1=0; nm1<NNormModes(); nm1++)
    for (int nm2=0; nm2<NNormModes(); nm2++)
      {
	norm_ini = 0;
	norm_targ = 0;
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    {
	      //x1*x1+y1*y1+..
	      norm_ini+= getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i) * getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i);
	      //x2*x2+y2*y2+..
	      norm_targ+= other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i) * other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i);
	      //x1*x2+y1*y2+...
	      overlap.Elem2(nm1,nm2)+= getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i) * other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i);
	    }
	overlap.Elem2(nm1,nm2)/= sqrt(norm_ini) * sqrt(norm_targ);
      }
  bool return_bool=true;


  // scan the diagonal; if the diagonal element is NOT the maximum element in the raw and the colum -- than ADD this nm to the list;
  double max;
  for (int nm=0; nm<NNormModes(); nm++)
    {
      max=fabs(overlap.Elem2(nm,nm)); // maximum should be at the diagonal
      for (int nm_scan=0; nm_scan<NNormModes(); nm_scan++) // scan row&column
	if ( (fabs(overlap.Elem2(nm,nm_scan))>max) or (fabs(overlap.Elem2(nm_scan,nm))>max) ) // check nm-th row & column
	  {
	    return_bool=false;
	    normal_modes_list.push_back(nm);
	    normal_modes_list.push_back(nm_scan); // the column/row with larger element should be also included
	  }
    }

  // Now Sort and Remove duplicates from the normal_modes_list:
  std::sort( normal_modes_list.begin(), normal_modes_list.end() );
  std::vector<int>::iterator new_end_pos;
  new_end_pos = std::unique( normal_modes_list.begin(), normal_modes_list.end() );
  normal_modes_list.erase( new_end_pos, normal_modes_list.end() );

  return return_bool;
 }



void MolState::Print()
{
  std::cout << "Geometry=\n";
  printGeometry();    
  std::cout << "Normal modes=\n";
  for (int k=0; k<NNormModes(); k++)
    {
      std::cout << "   Frequency="; 
      std::cout << getNormMode(k).getFreq() << '\n';
      std::cout << "   Displacement=\n";
      for (int i=0; i<NAtoms(); i++)
	{
	  for (int l=0; l<CARTDIM; l++)
	    std::cout << getNormMode(k).getDisplacement()[i*CARTDIM+l] << ' ';
	  std::cout << '\n';
	}
    }
  std::cout <<" end of the electronic state \n";   
}

void MolState::printGeometry()
{
  for (int i=0; i<NAtoms(); i++)
    {
      std::cout << std::setw(4) << std::right  << getAtom(i).Name();
      for (int k=0; k<CARTDIM; k++)
	std::cout << std::setw(12) << std::right << std::fixed << std::setprecision(4) << std::showpoint << getAtom(i).Coord(k) << ' '; 
      std::cout << '\n';
    }
}


void MolState::printNormalModes()
{
  // print in qchem format (3 per line)
  int nModesPerLine=3;  // three number of vib. modes per Line
  int nLines = NNormModes() / nModesPerLine;
  if ( NNormModes() % nModesPerLine != 0)
    nLines++;

  for (int n = 0; n < nLines; n++)  // number of blocks with 3 norm.modes. ("lines")
    {
      int current_nModesPerLine = nModesPerLine;
      // for the last entree, nModesPerString may differ from 3
      if (nLines - 1 == n)
  	if ( NNormModes() % nModesPerLine != 0 )
	  current_nModesPerLine = NNormModes() % nModesPerLine;
 
      for (int a=0; a < NAtoms(); a++)   
	{
	  for (int j=0; j < current_nModesPerLine; j++) 
	    {
	      for (int k=0; k < CARTDIM; k++)
	      
		std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(3)
			  << getNormMode(n*nModesPerLine+j).getDisplacement()[a*CARTDIM+k]*sqrt( reduced_masses[n*nModesPerLine+j]/getAtom(a).Mass() )<<' ';
	      std::cout <<  "  ";
	    }
	  std::cout << "\n";
	}
      std::cout << "\n";
    }

  /*

  for (int nm=0; nm < NNormModes(); nm++)
    {
      for (int a=0; a < NAtoms(); a++) 
	{
	  for (int k=0; k < CARTDIM; k++)
	    std::cout << std::setw(10) << std::right << std::fixed << std::setprecision(4)
		      << getNormMode(nm).getDisplacement()[a*CARTDIM+k]*sqrt( reduced_masses[nm]/getAtom(a).Mass() )<<' ';
	  std::cout << '\n';
	}
      std::cout << '\n';
    }

  */


}


//------------------------------
bool MolState::ifLetterOrNumber(char Ch)
{
if ( ((int(Ch)<=int('Z'))&&(int(Ch)>=int('A')))  
  || ((int(Ch)<=int('z'))&&(int(Ch)>=int('a'))) 
  || ((int(Ch)<=int('9'))&&(int(Ch)>=int('0'))) )
  return true;
else return false;
}



