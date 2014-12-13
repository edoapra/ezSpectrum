#ifndef _constants_h_
#define _constants_h_

/*! \file constants.h
\brief Simple math functions and constants 
\ingroup METHODS
*/

//! PI CHECK THIS
#define PI 3.141592654
//! Dimension of Cartesian 3D world
#define CARTDIM 3
//! Normalized random numbers
#define RANDMAX ((double)pow(2.,31.)-1.)
//! random number in [0,1] 
#define RANDOM ((double)random()/RANDMAX) 

//! Pi
#define C_Pi  3.141592654
//! c (speed of light) cm/second
#define SPEEDOFLIGHT_CM_SEC  2.99792458E10
//! c (speed of light) bohr/a.time.u
#define SPEEDOFLIGHT_AU 137.0359996

#define BOHR_2_CM 0.00000000529177249
//! h, erg-seconds 
#define PLANKCONSTANT_ERGxSEC  6.626176E-27

#define KELVINS2EV 8.61735E-05

#define AMU2GRAM 1.6605655E-24
#define AU2ANGSTROM  0.52918
#define ANGSTROM2AU 1.889716
#define WAVENUMBERS2EV (1.0/8065.479)

#define HBAR 1.054572669E-34
#define C_IN_SI  2.99792458E8

#define AMU_2_ELECTRONMASS (1822.8884843)

//! if less than threshold, considered to be zero
#define COORDINATE_THRESHOLD 0.00000001 //Angstrom

#endif
