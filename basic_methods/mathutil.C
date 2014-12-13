/*! \file mathutil.C
\ingroup BASIC_METHODS
*/

#include "mathutil.h"

//! Factorial function
double Factorial(const int N)
{
  double tmp=double(abs(N));
  
  if(tmp==0)
    return 1.0;
  
  for(int j=abs(N)-1;j>1;j--)
    tmp*=j;
  
  return tmp;
}

//! Factorial function INTEGER
unsigned long FactorialInt(const int N)
{
  unsigned long tmp=abs(N);
  
  if(tmp==0)
    return 1.0;
  
  for(int j=abs(N)-1;j>1;j--)
    tmp*=j;
  
  return tmp;
}

//! Factorial ratio function n1!/(n1-n2)!=(n1-n2+1)*(n1-n2+2)*...*(n1)
unsigned long FactorialRatioInt(const int n1, const int n2)
{
  if ((n1<0)or((n1-n2)<0))
    {
      std::cout << "\n DEBUG message: FactorialRatio() call error.\n\n";
      exit(2);
    }

  if (n1==0)
    return 1;

  unsigned long tmp=1;
 
  for(int j=n1;j>n2;j--)
    tmp*=j;
  
  return tmp;
}


//! Combination function
//long Combination(const int n, const int k)
//{
//  return  long(Factorial(n)/(Factorial(k)*Factorial(n-k)));
//}

//! Combination function (both parts should be less than max-unsigned-int)
unsigned long Combination(const int n, const int k)
//  return  long(Factorial(n)/(Factorial(k)*Factorial(n-k)));
{
  if (k>(n-k))
      return  FactorialRatioInt(n,k)/FactorialInt(n-k);
  else
      return  FactorialRatioInt(n,n-k)/FactorialInt(k);
}



//! provide next combination once previous is given in j (total number is C_n^k)
bool enumerateCombinations(int n, int k, std::vector <int>& j)
/*-------------------------------------------------------------------
NOTE: k>0!!

 Description:
   Enumerates all possible combinations of choosing k objects from n
   distint objects.  Initialize the enumeration by setting j[0] to a
   negative value.  Then, each call to enumerateCombinations will
   generate the next combinations and place it in j[0..k-1].  A
   return value of false indicates there are no more combinations to
   generate.  j needs to be allocated with a size for at least k
   elements.

 Usage:
   int comb[10] = {-1};
   while(enumerateCombinations(10, 5, comb)){
       // do something with comb[0..4]
   }

 Reference:
   Tucker, Allen.  Applied Combinatorics.  3rd Ed.  1994.
-------------------------------------------------------------------*/
{
    int i;
    if(j[0] < 0){
        for(i = 0; i < k; ++i){
            j[i] = i;
        }
        return true;
    } else {
        for(i = k - 1; i >= 0 && j[i] >= n - k + i; --i){}
        if(i >= 0){
            j[i]++;
            int m;
            for(m = i + 1; m < k; ++m){
                j[m] = j[m-1] + 1;
            }
            return true;
        } else {
            return false;
        }
    }
}
