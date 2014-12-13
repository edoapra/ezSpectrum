#include "franck_condon.h"

void harmonic_FCf(KMatrix& FCf, double Mass, double dQ, 
		  double Nu_ini, double Nu_targ)
{

  //Zero the FCF matrix
  FCf.Set(0);

  double tmpProduct;
  
  // alpha and delta in the (E.Hutchisson,1930) paper
  double a, d;

  dQ *= 1E-8;  // convert from Angstroms/XXX to centimiters/XXX
  Mass *= AMU2GRAM;  // convert amu to grams

  a = sqrt(Nu_ini/Nu_targ);  // alpha
  d = dQ * 2 * PI * sqrt (SPEEDOFLIGHT_CM_SEC*Nu_ini*Mass/PLANKCONSTANT_ERGxSEC ); //delta

  for (int n_ini=0; n_ini < FCf.Dim1(); n_ini++ )  // for each vibr level of the initial state
    {
      for (int n_targ=0; n_targ < FCf.Dim2(); n_targ++ ) // for each vibr level of the target state
    	{
  	  FCf.Elem2(n_ini,n_targ)= 0;
	  for (int L=0; L <= MIN (n_ini, n_targ); L++)
	    {
	      for (int I=0; I <= ( ((n_targ-L)%2==0) ? (n_targ-L)/2 : (n_targ-L-1)/2 ) ; I++)
		{
		  for (int J=0; J <=  ( ((n_ini-L)%2==0) ? (n_ini-L)/2 : (n_ini-L-1)/2) ; J++)
		    {
		      tmpProduct=1;
		      tmpProduct *= (1/Factorial(L))*pow( 4*a / (1+a*a), L);                         // * a_2L
		      tmpProduct *= (1/Factorial(I))*pow( (1-a*a) / (1+a*a), I);                   // * b_2I
		      tmpProduct *= (1/Factorial(J))*pow(-(1-a*a)/(1+a*a), J);                      //* c_2J
		      tmpProduct *= (1/Factorial(n_targ-2*I-L))*pow( -2*a*d/(1+a*a), n_targ-2*I-L ); // * d_xxx
		      tmpProduct *= (1/Factorial(n_ini -2*J-L))*pow(  2  *d/(1+a*a), n_ini -2*J-L ); // * e_xxx
		      FCf.Elem2(n_ini,n_targ)+=tmpProduct ;
		    }
		}
	    }
	  FCf.Elem2(n_ini,n_targ)*= sqrt(Factorial(n_ini)*Factorial(n_targ)) / pow( 2, (n_ini+n_targ)/2.0 ); // * factorial fraction
	  // multiplication by C3 (next two lines):
	  FCf.Elem2(n_ini,n_targ)*= exp( -d*d / (2*(a*a+1))  );
	  // that's like in the book:
	  // FCf.Elem2(n_ini,n_targ)*= sqrt (PLANKCONSTANT_ERGxSEC / (2*PI*PI*Mass*Nu_targ*a*(a*a+1)) ); 

	  // here what one gets for 0-0 transition (overlap of two ground state vibr.w.f., i.e gaussians, with different frequencies)
	  FCf.Elem2(n_ini,n_targ)*=sqrt(a*2/(a*a+1));
	}
    }
}


