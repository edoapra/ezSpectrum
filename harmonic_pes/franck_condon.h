#ifndef _franck_condon_h_
#define _franck_condon_h_

/*! \file franck_condon.h 
\ingroup METHODS
\brief Franck-Condon factors (vm'06)
*/

#include "genincludes.h"
#include "mathutil.h"
#include "constants.h"
#include "kmatrix.h"

/*!
This functions follows the algorithm described by Elmer Hutchisson in Phys. Rev. 36: 410-420 (1930) \n
FCFs are correct in absolute values (i.e. overlap of normolized vawefunctions)
Mass - is a reduced mass in amu, \n
dQ - displacement in angstoms (or in angstroms/sqrt(amu) -- in this case Mass should be =1amu) \n
Nu_ini/Nu_targ -- initial and target state's frequencies in cm-1  
*/

//! Analytic harmonic 1D Franck-Condon Factors.
void harmonic_FCf (KMatrix& FCf, double reducedMass, double deltaQ, 
                   double Nu_initial, double Nu_target);


#endif
