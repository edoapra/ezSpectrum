#ifndef _vibronic_state_h_
#define _vibronic_state_h_

/*! \file vibronic_state.h
\brief Vibronic state, 
i.e. stores electronic state index 
and excitation quantum number along each normal mode;

\ingroup MOLECULAR_PROP
*/

#include "genincludes.h"
#include <vector>

class VibronicState
{
  //! electronic state index: 0=initial state, 1,2... target states
  int elStateIndex; 
  //! Vibrational state quantum numbers, e.g. (0,2,1) for three normal modes
  std::vector<int> vibrQuanta;
  //! nombers of normal modes; e.g. if the "excite subspace" is {v3,v7,v12}, than (0,2,1) excitaion will be (2v7, 1v12) in the full space
  std::vector<int> excite_subspace;

 public:
  VibronicState() {elStateIndex=0; vibrQuanta.clear(); };

  //??
  void reset() {elStateIndex=0; vibrQuanta.assign(vibrQuanta.size(),0);};
  void clear() {elStateIndex=0; vibrQuanta.clear(); excite_subspace.clear(); };

  int getElStateIndex() { return elStateIndex; };
  void setElStateIndex(const int index) { elStateIndex = index; };

  int getNNormModes() { return vibrQuanta.size(); };
  int getVibrQuantaSize() { return vibrQuanta.size(); };

  // adds a quanta the normal mode number
  void addVibrQuanta(const int quantNumber, const int nm) { vibrQuanta.push_back(quantNumber); excite_subspace.push_back(nm); };
  int getVibrQuanta(const int i) { return vibrQuanta[i]; };
  void setVibrQuanta(const int i, const int quantNumber) { vibrQuanta[i]=quantNumber; };

  //! checks if nm is in the "excite subspace", otherwise returns zero
  int getV_full_dim(const int nm);

  //!returns the number of the normal mode nm (nm in the full space) in the "excite subspace"
  //int getEx_full_dim(const int nm);

  int getTotalQuantaCount();

  void incrIndex(const int index, const int add){vibrQuanta[index]+=add;};

  //avoid this:
  std::vector<int>& getV(){return vibrQuanta;};
  std::vector<int>& getEx(){return excite_subspace;};


  bool if_equal(VibronicState& other);

  void print();
  //void print_full();
};

#endif
