#include "spectrum.h"


void Spectrum::PrintStickTable()
{
  std::cout << "Energy,eV    Intensity     FC factor         E\",K   Transition\n\n" ;
  for (int i=0; i<spectralPoints.size(); i++)
    if ( spectralPoints[i].getIfPrint() )
      getSpectralPoint(i).print();
}


void Spectrum::PrintStickTable( const char* spectrumFileName )
{
  int quanta_printed;
  std::ofstream spectrumF;     
  spectrumF.open(spectrumFileName, std::ios::out);
  for (int i=0; i<spectralPoints.size(); i++)
    if ( spectralPoints[i].getIfPrint() )
      {
	spectrumF << std::fixed << std::setprecision(4) << std::setw(7) << getSpectralPoint(i).getEnergy();
	spectrumF << std::scientific << std::setprecision(6) << std::setw(18) << getSpectralPoint(i).getIntensity() << "  ";
	spectrumF << std::scientific << std::setprecision(6) << std::setw(11) << std::showpos << getSpectralPoint(i).getFCF()<< std::noshowpos;
	spectrumF << std::fixed << std::setprecision(3) << std::setw(10) << getSpectralPoint(i).getE_prime_prime()/KELVINS2EV << "  ";
	
	// print initial & target states
	spectrumF << getSpectralPoint(i).getVibrState1().getElStateIndex() << '(';
	quanta_printed=0;
	for (int n=0; n<getSpectralPoint(i).getVibrState1().getVibrQuantaSize(); n++)
	  if (getSpectralPoint(i).getVibrState1().getVibrQuanta(n)>0)
	    {
	      if (quanta_printed>0) // not the first normal mode in the vector
		spectrumF << ',';
	      quanta_printed++;
	      spectrumF << getSpectralPoint(i).getVibrState1().getVibrQuanta(n)<<'v'<<getSpectralPoint(i).getVibrState1().getEx()[n];
	    }
	if (quanta_printed==0)
	  spectrumF << '0';
	spectrumF << ")->" << getSpectralPoint(i).getVibrState2().getElStateIndex() << '(';
	quanta_printed=0;
	for (int n=0; n<getSpectralPoint(i).getVibrState2().getVibrQuantaSize(); n++)
	  if (getSpectralPoint(i).getVibrState2().getVibrQuanta(n)>0)
	    {
	      if (quanta_printed>0) // not the first normal mode in the vector
		spectrumF << ',';
	      quanta_printed++;
	      spectrumF << getSpectralPoint(i).getVibrState2().getVibrQuanta(n)<<'v'<<getSpectralPoint(i).getVibrState2().getEx()[n];
	    }
	if (quanta_printed==0)
	  spectrumF << '0';
	spectrumF << ")\n";
      }
  spectrumF.close();
}


void Spectrum::AddSpectralPoint( SpectralPoint& PointToAdd)  
{ 
  spectralPoints.push_back(PointToAdd); 
};


void Spectrum::AddSpectralPoint(const double E, const double I, const double fcf, const double Epp, 
				VibronicState state_ini, VibronicState state_targ, const bool if_print_flag)
{
  SpectralPoint tmpPoint;
  for (int nm=0; nm<state_ini.getVibrQuantaSize(); nm++)
    tmpPoint.getVibrState1().addVibrQuanta( state_ini.getV()[nm], state_ini.getEx()[nm] );
  for (int nm=0; nm<state_targ.getVibrQuantaSize(); nm++)
    tmpPoint.getVibrState2().addVibrQuanta( state_targ.getV()[nm], state_targ.getEx()[nm] );
  tmpPoint.getVibrState1().setElStateIndex(state_ini.getElStateIndex());
  tmpPoint.getVibrState2().setElStateIndex(state_targ.getElStateIndex());

  tmpPoint.getEnergy() = E;
  tmpPoint.getIntensity() = I; 
  tmpPoint.getFCF()= fcf;
  tmpPoint.getE_prime_prime()= Epp;

  tmpPoint.setIfPrint(if_print_flag);

  AddSpectralPoint(tmpPoint);
}
