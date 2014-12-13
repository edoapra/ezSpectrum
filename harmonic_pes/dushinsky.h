#ifndef _dushinsky_h_
#define _dushinsky_h_

/*! \file dushinsky.h 
\ingroup METHODS
\brief calculates and stores full dimentional Franck-Condon factors (i.e. including Dushinsky rotations) 

the notation and equations are from [Berger et al. JPCA 102:7157(1998)]
To initialize provide:
1. vector of molStates; ([0]-initial, [1],[2],... --target states)
2. vector of int -- list of normal modes in the subspace;
3. threshold for FCFs to accept (in_fcf_threshold)
4. target state index =1..Nmax for spectral lines assignment (in_targN)

*/

#include "genincludes.h"
#include "mathutil.h"
#include "constants.h"
#include "kmatrix.h"
#include "vibronic_state.h"
#include "vibrational_indexing.h"
#include <vector>
#include "spectrum.h"
#include "molstate.h"


// IMPORTANT! All matrices and zero-zero integral are evaluated for the FULL space of the normal modes of the molecule
// than space is shrinked to the space described by nm_list, so excitations will be added only to those normal modes.

class Dushinsky
{
  //! number of normal modes (dimentionality) 
  int N;
  //! frequently used matrices in the "excite subspace"
  KMatrix ompd, tpmo, rd, tqmo, tr;
  //! frequently used matrices in the full space; requred for single excitation (outside the "excite subspace") and hot bands recursive calculations;
  KMatrix ompd_full, tpmo_full, rd_full, tqmo_full, tr_full;
  //! K' -- maximum layer which was stored ( maximum total number of quanta in the target state):
  int Kp_max_saved;
  //! K' -- maximum layer which was evaluated( maximum total number of quanta in the target state):
  int Kp_max;
  //! <zero|zero> integral -- initital one in the recurrent procedure
  double zero_zero;
  //! threshold for FCFs to be included in the spectrum:
  double fcf_threshold;
  //! target state index
  int targN;

  //! initial('constant' for the no hot bands case) and target(variable) states:
  VibronicState state0, state;

  //! "layer" #K' (K'=0...Kp_max)contains all states with the total number of excitations = K' in the target state only
  // each layer is a linear vector; to get between "vector's index" and "vibrational state" 
  // functions convVibrState2Index() and convIndex2VibrState are used
  std::vector< std::vector<double>* > layers;

  //! spectrum which stores all the points
  Spectrum spectrum;


  // matrix C: combinations C(n,k)=C_n^k=n!/(k!*(n-k)!); 
  // n=0..N+K-1; k=0..N-1=0..min(n,N-1); matrix (N+K)x(N) -- N is known everywhere;
  // Combination(n,k)=C(n,k)=C[n*N+k]
  // N-total number of normal modes; K -- max number of quanta in target state
  // ZZZ  4/11/2012 removed, and the combinations are calculated now on the fly
  // unsigned long *C;

  // array of sqrt: sqrtArray[i]=sqrt(i); 
  double *sqrtArray;

 public:
  Dushinsky(std::vector <MolState>& molStates, std::vector<int>& nm_list, double in_fcf_threshold, int in_targN, int max_quanta_target, int max_quanta_initial);
  ~Dushinsky();

  //! K'=0 is zero_zero; so it starts with K'=1 and increments Kp_max; Also updates the spectrum and returns number of points below the threahold added;
  int evalNextLayer(const bool if_save=true);

  //! add hot bands to the spectrum
  int addHotBands(std::vector <MolState>& molStates, std::vector<int>& nm_list, 
		   double fcf_threshold, double temperature, 
		   int max_n_initial, int max_n_target, 
		   double energy_threshold_initial,  double energy_threshold_target);

  //! if total number excitations in the stete is larger than Kp_max (stored), evaluates FCF recursively; 
  double evalSingleFCF(VibronicState& state_ini, int K, VibronicState& state_targ, int Kp);
  //! a separate copy of evalSingleFCF() for faster execution; if there are excitations in normal modes outside "excite subspace", than use this one;
  double evalSingleFCF_full_space(VibronicState& state_ini, int K, VibronicState& state_targ, int Kp);


  //! returns the up to date spectrum
  Spectrum& getSpectrum(){return spectrum;};


  //! add a single point to the spectrum
  void addSpectralPoint(const double fcf, VibronicState state_ini, VibronicState state_targ );
  //! returns Kp_max
  int getLmax(){return Kp_max;};
  //! Prints estimated sizes of each layer up to Kp_max (works before the memory is actually allocated)
  void  printLayersSizes(const int uptoKp);

};

#endif
