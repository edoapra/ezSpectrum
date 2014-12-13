#include "spectralpoint.h"

void SpectralPoint::print()
{
  std::cout << std::fixed << std::setprecision(4) << std::setw(7) << E;
  std::cout << std::scientific << std::setprecision(6) << std::setw(18) << I << "  ";
  std::cout << std::scientific << std::setprecision(6) << std::setw(11) << std::showpos << FCF << std::noshowpos ;
  std::cout << std::fixed << std::setprecision(3) << std::setw(10) << Epp/KELVINS2EV << "  ";
  State1.print();
  std::cout << "->";
  State2.print();
  std::cout << "\n";
}


/*
void SpectralPoint::ini(const double E, const double I, const double fcf, const double Epp, VibronicState state_ini, VibronicState state_targ )
{
  getEnergy() = E;
  getIntensity() = I; 
  getFCF()= fcf;
  getE_prime_prime()= Epp;
  getVibrState1().reset();
  getVibrState1().setElStateIndex(state_ini.getElStateIndex());
  getVibrState2().reset();
  getVibrState2().setElStateIndex(state_targ.getElStateIndex());


  for (int nm=0; nm<state_ini.getVibrQuantaSize(); nm++)
    getVibrState1().setVibrQuanta(nm, state_ini.getVibrQuanta(nm));
  for (int nm=0; nm<state_targ.getVibrQuantaSize(); nm++)
    getVibrState2().setVibrQuanta(nm, state_targ.getVibrQuanta(nm));

}
*/
