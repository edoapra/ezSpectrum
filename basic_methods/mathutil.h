#ifndef _mathutil_h_
#define _mathutil_h_

/*! \file mathutil.h
\brief General math utilities: fact, MIN, MAX, SIGN, RANDOM etc.
*/

#include "genincludes.h"
#include <vector>

//! Factorial; N=(-inf..0..+inf)  
double Factorial(const int N); 

//! Factorial function INTEGER
unsigned long FactorialInt(const int N);

//! Factorial ratio function n1!/(n1-n2)!=(n1-n2+1)*(n1-n2+2)*...*(n1)
unsigned long  FactorialRatioInt(const int n1, const int n2);

//! Combination, C_n^k  
unsigned long Combination(const int n, const int k); 

//! provide next combination once previous is given in j (total number is C_n^k) -- see comments in  .C
bool enumerateCombinations(int n, int k, std::vector <int>& j);

//Max and Min f-s
//! A macro that returns the maximum of a and b
#define MAX(a,b) ( (a>=b) ? (a) : (b) )
//! A macro that returns the minimum of a and b
#define MIN(a,b) ( (a<b) ? (a) : (b) ) 
//! A macro that returns the sign of as a (double)+/-1.
#define SIGN(a)  ( (a<0.) ? (-1.) : (1.) )
//! random number in [0,1] 

#define RANDOM ((double)random()/RANDMAX) 
//! Randomizer (use system time)
void Randomize(void );


//! Making sure bool, TRUE and FALSE are defined
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (1)
#endif


#endif


