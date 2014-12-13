#include "vibronic_state.h"

bool  VibronicState::if_equal(VibronicState& other)
{
  bool return_bool = true;
  if(elStateIndex!=other.elStateIndex)
    return_bool = false;

  if (getVibrQuantaSize()!=other.getVibrQuantaSize())
    return_bool = false;
  else
    for (int nm=0; nm<getVibrQuantaSize(); nm++)
      if (getVibrQuanta(nm)!=other.getVibrQuanta(nm))
	return_bool=false;

  return return_bool;
}

int VibronicState::getTotalQuantaCount()
{
  int count=0;
  for (int nm=0; nm<getV().size(); nm++)
    count+=getV()[nm];
  return count;
}

int VibronicState::getV_full_dim(const int nm)
{
  int return_quanta = 0;
  // check if nm in the "excite space"
  for (int i=0; i<getEx().size(); i++)
    if (nm==getEx()[i])
      return_quanta = getV()[ i ];

  return return_quanta;
}


/*
//!returns the number of the normal mode nm (nm in the full space) in the "excite subspace"
int VibronicState::getEx_full_dim(const int nm)
{
  int return_number = -1; // -1 -- error, not found;
  // check if nm in the "excite space"
  for (int i=0; i<getEx().size(); i++)
    if (nm==getEx()[i])
      return_number = nm;

  return return_number;
}
*/

//NEW format
void VibronicState::print()
{
  int quanta_printed=0;
  // print in format N(a,b,c,d..), where N is the electronic state index, a,b,c..--vibrational quanta
  std::cout<<elStateIndex<<'(';

  for (int nm=0; nm<excite_subspace.size(); nm++)
    {
      if (vibrQuanta[nm]>0)
	{
	  if (quanta_printed>0) // not the first normal mode in the vector
	    std::cout << ',';
	  quanta_printed++;
	  std::cout<<vibrQuanta[nm]<<'v'<<excite_subspace[nm];
	}
    }

  // if now excitations, than it is the ground vibbr.  state, i.e. "(0)"    
  if (quanta_printed==0)
    std::cout<<'0';
  std::cout << ")";
}

/*
//OLD format ( NOT UPDATED FOR SUBSPACE USE )
void VibronicState::print_full()
{
  // print in format N(a,b,c,d..), where N is the electronic state index, a,b,c..--vibrational quanta
  std::cout<<elStateIndex<<'(';
  // if full printout:
  if (excite_subspace.size()==0)
    {
      for (int i=0; i<vibrQuanta.size(); i++)
	{
	  std::cout<<vibrQuanta[i];
	  if (i<vibrQuanta.size()-1) // not the last normal mode in the vector
	    std::cout << ",";
	}
    }
  else
    {
      for (int i=0; i<excite_subspace.size(); i++)
	{
	  std::cout<<vibrQuanta[excite_subspace[i]];
	  if (i<excite_subspace.size()-1) // not the last normal mode in the vector
	    std::cout << ",";
	}
    }
  std::cout << ")";
}

*/
