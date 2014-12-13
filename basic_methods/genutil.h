#ifndef _genutil_h_
#define _genutil_h_

/*! \file genutil.h
\brief General utilities (get Time; error messages )  
\ingroup BASIC_METHODS
*/

#include <time.h>
#include "genincludes.h"


//! Error handling (print msg, exit(1))
void error(const char *const msg);

//! if expr -> error(msg);
inline void check(const bool expr, const char *const msg)
{
  if (expr) error(msg);
}

//! Returns string with the current time:
std::string GetTime();
time_t GetRawTime();

//! Swap two integers
void Swap(int& i, int& j);

//! Swap two doubles
void Swap(double& i, double& j);


#endif
