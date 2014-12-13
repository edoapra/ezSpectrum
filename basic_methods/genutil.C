/*! \file genutil.C
\ingroup BASIC_FUNCTIONS
*/

#include "genutil.h"
#include "genincludes.h"

// Return string with the current time:
std::string GetTime(){
  time_t rawtime;
  time ( &rawtime );
  return  ctime ( &rawtime );
}

time_t GetRawTime(){
  time_t rawtime;
  time ( &rawtime );
  return  rawtime;
}


void error(const char *const msg)
{
  std::cout << "\nError! " << msg <<"\n\n";
  exit(1);
}

//! Swap two integers
void Swap(int& i, int& j)
{
  static int tmp;
  tmp=i;
  i=j;
  j=tmp;
}

//! Swap two doubles
void Swap(double& i, double& j)
{
  static double tmp;
  tmp=i;
  i=j;
  j=tmp;
}
