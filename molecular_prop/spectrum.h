#ifndef _spectrum_h_
#define _spectrum_h_

/*! \file spectrum.h
\brief Spectrum object
Implemented:
1) stores stick spectrum (vector of "SpectralPoint"s {I,E}) ) 
2) sorts it by energy
3) prints as a table

Sould be done soon:
1) prints as a spectrum with some bradening..
2)

\ingroup MOLECULAR_PROP
*/
#include "genincludes.h"
#include "spectralpoint.h"
#include <vector>
#include <algorithm>

class Spectrum
{
  std::vector<SpectralPoint> spectralPoints;

  struct SortByEnergy
  {
    bool operator()(SpectralPoint p1, SpectralPoint p2)
    {
      return (p1.getEnergy() > p2.getEnergy());
    };
  };

 public:
  Spectrum(){};
  
  SpectralPoint& getSpectralPoint(int i) { return spectralPoints[i]; };
  int getNSpectralPoints() {return spectralPoints.size(); };

  void AddSpectralPoint( SpectralPoint& PointToAdd);
  void AddSpectralPoint(const double E, const double I, const double fcf, const double Epp, 
			VibronicState state_ini, VibronicState state_targ, const bool if_print_flag=true );


  void Sort() { sort( spectralPoints.begin(), spectralPoints.end(), SortByEnergy() );};

  void PrintStickTable();
  void PrintStickTable(const char* spectrumFileName);
 
};

#endif

