#ifndef _spectral_point_h_
#define _spectral_point_h_

/*! \file spectralpoint.h
\brief Spectral point I(E)
stores pair Intensity-Energy and the transition (vibrational quanta of the initial and target states)

\ingroup MOLECULAR_PROP
*/

#include "constants.h"
#include "vibronic_state.h"

class SpectralPoint
{
  double I;
  double E;

  //energy of the hot bands (vibr. excited states of the initial electronic state)
  double Epp; //E''

  // Franck-Condon factor;
  double FCF;

  bool if_print;

  VibronicState State1, State2;
  
 public:
  SpectralPoint(){ if_print=true; }; // by default the point is active(visible)

  double& getIntensity() { return I; };
  double& getEnergy() { return E; };

  double& getE_prime_prime() { return Epp; };
  double& getFCF() { return FCF; };

  bool getIfPrint() {return if_print; };
  void setIfPrint(const bool flag) {if_print=flag; };

  VibronicState& getVibrState1() { return State1; };
  VibronicState& getVibrState2() { return State2; };

  //void ini(const double E, const double I, const double fcf, const double Epp, VibronicState state_ini, VibronicState state_targ );

  void print();
};

#endif
