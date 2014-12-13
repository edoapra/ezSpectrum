#ifndef _harmonic_pes_main_
#define _harmonic_pes_main_

/*! \file harmonic_pes__main.h
\brief  the "main" for analytic-harmonic Franck-Condon Factors calculation
*/


//////////////////////////////////////////////////////////////////////////////
//                                                                          // 
//  calculates photoelectron spectrum in 1D Franck-Condon approximation      // 
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

//2DO: 
// 1 ADD DEGENERATE NORMAL MODES
// 2 several target states
// cut into nice picies

//! .xml file name with atomName<->atomMass table
#define ATOMIC_MASSES_FILE ("atomicMasses.xml")

#include "genincludes.h"
#include <vector>
#include <set>

#include "molstate.h"
#include "simple_xml_parser.h"
#include "kmatrix.h"
#include "vector3d.h"

#include "parallel_approximation.h"
#include "dushinsky.h"
#include "vibrational_indexing.h"

#include <limits>

//! program itself
bool harmonic_pes_main (const char* xmlFileName);


//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string& ex_str, simpleXMLparser& xmlF, int& qnt, int& nm );
//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state (i.e. vector of integers)
void fillVibrState(My_istringstream& vibr_str, VibronicState& v_state, simpleXMLparser& xmlF, const int nm_max );

#endif


