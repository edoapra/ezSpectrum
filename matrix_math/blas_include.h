#ifndef _blas_include_h
#define _blas_include_h

/*! 
\file blas_include.h
\brief This file contains our C-wrappers for BLAS/LAPACK-function calls.
\ingroup (MATRIX_MATH)

  Wrappers are in the file "blas_calls.C". All documentation in that 
  file is from the FORTRAN descriptions: to call these functions from C 
  take care about rows<->columns.

  To emphasize that these are literal C wrappers, and not C translations
  of the Fortran calls, we have added a prefix "CL_" for "C Literal".
  The prefix also helps avoid name conflicts on CRAY machines where
  the Fortran library uses all-caps, like DAXPY).
  David Sherrill, 17 June 1999

 AIK, 1998.
*/

void CL_DAXPY(int ntot, double coeff, double *copy_from, int inc1,
              double *copy_to, int inc2);
void CL_DCOPY(int ntot,double *copy_from,int inc1,double *copy_to,int inc2);
int CL_IDAMAX(int ntot, double *data, int inc);
//int CL_IDAMIN(int ntot, double *data, int inc);
double CL_DDOT(int ntot, double *x, int incx, double *y, int incy);
void CL_DGEMM(char transa, char transb, int m, int n, int k, double alpha, 
              double *A, int lda, double *B, int ldb, double beta, double *C, 
              int ldc);
void CL_DGESVD(char jobu, char jobvt, int m, int n, double *A, int lda,
               double *S, double *U, int ldu, double *VT, int ldvt, 
               double *work, int lwork, int& info);
void CL_DSCAL(int size, double coeff, double *data, int inc);	    
void CL_DSYEV(char jobs, char type, int dim, double *matrix, int ld,
	      double *diag, double *work, int lwork, int& info);
void CL_DGEEV(char jobvl, char jobvr, int dim, double *A, double *wr, 
              double *wi, double *vl, double *vr, double *work, 
              int lwork, int& info);
void CL_DSYSV(char type,int dim, int nrhs, double *A, int lda, int *ipiv,
              double *B, int ldb, double *work, int lwork, int& info);	   
void CL_DROT(int ntot, double *x, int incx, double *y, int incy, 
             double costheta, double sintheta);
void CL_DGESV(int dim, int nrhs, double *A, int lda, int *ipiv,
              double *B, int ldb, int& info);
void CL_DSYGV(int itype, char jobz, char uplo, int dim, double *A, int lda,
	      double *B, int ldb, double *ene, double *work, int lwork,
	      int& info);

void CL_DGETRF(int M, int N, double *A, int lda, int *ipiv,  int& info);
void CL_DGETRI(int N, double *A, int lda, int *ipiv, double *work, int lwork, int& info);


#endif

#if 0 // uncomment to inline have blas calls 

//#include "blas_include.h"
#ifdef QCHEM
#define EXTREF
#include "qcmacdep.h"
#endif

//Perhaps in the future we should avoid using this file and
//use the version in libqt directly for PSI.  How to do that
//and keep this around for Q-Chem? --CDS 3/00
#ifdef PSI
#if FCLINK==3
  #define CAPS
#elif FCLINK==1
  #define UNDERSCORE
#endif
#endif

#if defined CAPS          // CRAY, etc.
  #define F_DAXPY DAXPY
  #define F_DCOPY DCOPY
  #define F_IDAMAX IDAMAX
  #define F_IDAMIN IDAMIN
  #define F_DDOT DDOT
  #define F_DGEMM DGEMM
  #define F_DGESVD DGESVD
  #define F_DSCAL DSCAL
  #define F_DSYEV DSYEV
  #define F_DGEEV DGEEV
  #define F_DSYSV DSYSV
  #define F_DROT DROT
  #define F_DGESV DGESV
  #define F_DSYGV DSYGV
  #define F_DGETRF DGETRF
  #define F_DGETRI DGETRI
#elif defined UNDERSCORE  // GNU, most systems
  #define F_DAXPY daxpy_
  #define F_DCOPY dcopy_
  #define F_IDAMAX idamax_
  #define F_IDAMIN idamin_
  #define F_DDOT ddot_
  #define F_DGEMM dgemm_
  #define F_DGESVD dgesvd_
  #define F_DSCAL dscal_
  #define F_DSYEV dsyev_
  #define F_DGEEV dgeev_
  #define F_DSYSV dsysv_
  #define F_DROT drot_
  #define F_DGESV dgesv_
  #define F_DSYGV dsygv_
  #define F_DGETRF dgetrf_
  #define F_DGETRI dgetri_
#else                     // AIX, etc.
  #define F_DAXPY daxpy
  #define F_DCOPY dcopy
  #define F_IDAMAX idamax
  #define F_IDAMIN idamin
  #define F_DDOT ddot
  #define F_DGEMM dgemm
  #define F_DGESVD dgesvd
  #define F_DSCAL dscal
  #define F_DSYEV dsyev
  #define F_DGEEV dgeev
  #define F_DSYSV dsysv
  #define F_DROT drot
  #define F_DGESV dgesv
  #define F_DSYGV dsygv
  #define F_DGETRF dgetrf
  #define F_DGETRI dgetri
#endif


extern "C"
{
  void F_DAXPY(int *ntot, double *coeff, double *copy_from, int *inc1,
	        double *copy_to, int *inc2);
  void F_DCOPY(int *ntot,double *copy_from,int *inc1, double *copy_to,
                int *inc2);
  int F_IDAMAX(int *ntot, double *data, int *inc);
  int F_IDAMIN(int *ntot, double *data, int *inc);
  double F_DDOT(int *ntot, double *x, int *incx, double *y, int *incy);
  void F_DGEMM(char *transa, char *transb, int *m, int *n, int *k, 
                double *alpha, double *A, int *lda, double *B, int *ldb, 
                double *beta, double *C, int *ldc);
  void F_DGESVD(char *jobu, char *jobvt, int *m, int *n, double *A, int *lda,
	         double *S, double *U, int *ldu, double *VT, int *ldvt,
	         double *work, int *lwork, int *info);
  void F_DSCAL(int *size,double *coeff,double *matrix,int *inc);
  void F_DSYEV(char *jobs, char *type, int *dim, double *matrix, int *ld,
               double *diag, double *work, int *lwork, int *info);
  void F_DGEEV(char *jobvl, char *jobvr, int *n, double *da, int *lda,
               double *wr, double *wi, double *vl, int *ldvl, double *vr,
               int *ldvr, double *work, int *lwork, int *info);
  void F_DSYSV(char *type, int *dim, int *nrhs, double *A, int *lda, 
               int *ipiv, double *B, int *ldb, double *work, int *lwork, 
               int *info);	     
  void F_DROT(int *ntot,double *x, int *incx,double *y, int *incy,
              double *cotheta,double *sintheta);
  void F_DGESV(int *dim, int *nrhs, double *A, int *lda,
               int *ipiv, double *B, int *ldb, int *info);
  void F_DSYGV(int *itype, char *jobz, char *uplo, int *dima, double *A,
	       int *lda, double *B, int *ldb, double *W, double *work,
	       int *lwork, int *info);
  void F_DGETRF(int *M, int *N, double *A, int *lda, int *ipiv, int *info);
  void F_DGETRI(int *N, double *A, int *lda, int *ipiv, double *work, int* lwork, int *info);
}



/*!
   void CL_DAXPY(int ntot, double coeff, double *copy_from, int inc1,
                 double *copy_to, int inc2);

   This function adds data from copy_from to copy_to.
 
   int ntot: length of data.

   int inc1,inc2: increments for copy_from, copy_to.
*/
inline void CL_DAXPY(int ntot, double coeff, double *copy_from, int inc1,
              double *copy_to, int inc2)
{
  F_DAXPY(&ntot,&coeff,copy_from,&inc1,copy_to,&inc2);
}


/*!
   void CL_DCOPY(int ntot, double *x, int incx, double *y, int *incy);

   This function copies x to y.
   
   int ntot: length of x,y;

   int incx,incy: increments for x,y.
*/
inline void CL_DCOPY(int ntot,double *copy_from,int inc1,double *copy_to,int inc2)
{
  F_DCOPY(&ntot,copy_from,&inc1,copy_to,&inc2);
}


/*!
   int CL_IDAMAX(int ntot, double *data, int inc);

   Returns adress of maximum value in data.
   
   int ntot: length of data.

   int inc: increments for data.
*/
int CL_IDAMAX(int ntot, double *data, int inc)
{
  return F_IDAMAX(&ntot,data,&inc);
}

/*!
   void CL_IDAMIN(int ntot, double *data, int inc);
  
   This function finds minimum element in *data.
   
   int ntot: length of xdata to scan.

   int inc: increments for data.
*/
/*
   Blas on IBM does not have it for some reason
   int CL_IDAMIN(int ntot, double *data, int inc)
   {
   return F_IDAMIN(&ntot,data,&inc);
   }
*/


/*!
   void CL_DDOT(int ntot, double *x, int incx, double *y, int *incy);
  
   This function calculates dot product (x,y).
   
   int ntot: length of x,y.

   int incx,incy: increments for x,y.
*/
double CL_DDOT(int ntot, double *x, int incx, double *y, int incy)
{
  return F_DDOT(&ntot,x,&incx,y,&incy);
}


/*!
   void CL_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
  	         double *A, int lda, double *B, int ldb, double beta, 
                 double *C, int ldc);
  
   This is a C-wrapper for the literal DGEMM Fortran call; note that the
      C caller must be aware of how to swap rows/columns to make the proper
      call to DGEMM; this can be confusing.
  
   This function calculates C(m,n)=alpha*(opT)A(m,k)*(opT)B(k,n)+ beta*C(m,n).
   
   char transa:       On entry, specifies the form of (op)A used in the
                      matrix multiplication:

                      If transa = 'N' or 'n', (op)A = A

                      If transa = 'T' or 't', (op)A = transp(A)

                      If transa = 'R' or 'r', (op)A = conjugate(A)

                      If transa = 'C' or 'c', (op)A = conjug_transp(A)

                      On exit, transa is unchanged.

   char transb:       On entry, specifies the form of (op)B used in the
                      matrix multiplication:

                      If transb = 'N' or 'n', (op)B = B

                      If transb = 'T' or 't', (op)B = transp(B)

                      If transb = 'R' or 'r', (op)B = conjugate(B)

   int m:             On entry, the number of rows of the matrix (op)A and of
                      the matrix C; m >= 0. On exit, m is unchanged.

   int n:             On entry, the number of columns of the matrix (op)B and
                      of the matrix C; n >= 0. On exit, n is unchanged.

   int k:             On entry, the number of columns of the matrix (op)A and
                      the number of rows of the matrix (op)B; k >= 0. On exit,
  		      k is unchanged.

   double alpha:      On entry, specifies the scalar alpha. On exit, alpha is
                      unchanged.

   double *A:         On entry, a two-dimensional array A with dimensions lda
                      by ka. For (op)A = A  or  conjugate(A), ka >= k and the 
                      leading m by k portion of the array A contains the matrix
                      A. For (op)A = transp(A) or conjug_transp(A), ka >= m
                      and the leading k by m part of the array A contains the
                      matrix A. On exit, a is unchanged.

   int lda:           On entry, the first dimension of array A.
                      For (op)A = A  or conjugate(A), lda >= MAX(1,m).
                      For (op)A=transp(A) or conjug_transp(A), lda >= MAX(1,k).
                      On exit, lda is unchanged.

   double *B:         On entry, a two-dimensional array B with dimensions ldb
                      by kb. For (op)B = B or conjugate(B), kb >= n and the
                      leading k by n portion of the array contains the matrix
                      B. For (op)B = transp(B) or conjug_transp(B), kb >= k and
                      the leading n by k part of the array contains the matrix
  		      B. On exit, B is unchanged.

   int ldb:           On entry, the first dimension of array B.
                      For (op)B = B or <conjugate(B), ldb >= MAX(1,k).
                      For (op)B = transp(B) or conjug_transp(B), ldb >=
                      MAX(1,n). On exit, ldb is unchanged.

   double *beta:      On entry, specifies the scalar beta. On exit, beta is
                      unchanged.

   double C:          On entry, a two-dimensional array with the dimension
                      ldc by at least n. On exit,  the leading  m by n part of
                      array C is overwritten by the matrix alpha*(op)A*(op)B +
                      beta*C.

   int ldc:           On entry, the first dimension  of array C; ldc >=MAX(1,m)
                      On exit, ldc is unchanged.
*/
inline void CL_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
              double *A, int lda, double *B, int ldb, double beta, double *C,
              int ldc)
{
  F_DGEMM(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}


/*!
   void CL_DGESVD(char jobu, char jobvt, int m, int n, double *A, int lda,
                  double *S, double *U, int ldu, double *VT, int ldvt, 
                  double *work, int lwork, int info)
   
   This function performs SVD-decomposition: A(m,n)=U*S*VT, FORTRAN notations.

   char jobu:  Specifies options for computing all or part of the matrix U:

                 = 'A':  all M columns of U are returned in array U:

                 = 'S':  the first min(m,n) columns of U (the left singular
  	                 vectors) are returned in the array U;

                 = 'O':  the first min(m,n) columns of U (the left singular
                         vectors) are overwritten on the array A;

  		 = 'N':  no columns of U (no left singular vectors) are
                         computed

   char jobvt:  Specifies options for computing all or part of the matrix V**T:

                 = 'A':  all N rows of V**T are returned in the array VT;

                 = 'S':  the first min(m,n) rows of V**T (the right singular
                         vectors) are returned in the array VT;

                 = 'O':  the first min(m,n) rows of V**T (the right singular
                         vectors) are overwritten on the array A;

                 = 'N':  no rows of V**T (no right singular vectors) are com-
                         puted.

   int m:       The number of rows of the input matrix A.  M >= 0.

   int n:       The number of columns of the input matrix A.  N >= 0.

   double *A:   Array, dimension (LDA,N). On entry, the M-by-N matrix A.
                On exit, if JOBU = 'O',  A is overwritten with the first
                min(m,n) columns of U (the left singular vectors, stored
                columnwise); if JOBVT = 'O', A is overwritten with the first
                min(m,n) rows of V**T (the right singular vectors, stored
                rowwise); if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of
                A are destroyed.

   int lda:     The leading dimension of the array A.  LDA >= max(1,M).

   double *S    Array, dimension (min(M,N)). The singular values of A, sorted
                so that S(i) >= S(i+1).

   double *U    Array, dimension (LDU,UCOL). (LDU,M) if JOBU = 'A' or
                (LDU,min(M,N)) if JOBU = 'S'.  If JOBU = 'A', U contains the
                M-by-M orthogonal matrix U; if JOBU = 'S', U contains the
                first min(m,n) columns of U (the left singular vectors,
                stored columnwise); if JOBU = 'N' or 'O', U is not referenced.

   int ldu:     The leading dimension of the array U.  LDU >= 1; if JOBU = 'S'
                or 'A', LDU >= M.

   double *VT:  Array, dimension (LDVT,N). If JOBVT = 'A', VT contains the
                N-by-N orthogonal matrix V**T; if JOBVT = 'S', VT contains
                the first min(m,n) rows of V**T (the right singular vectors,
                stored rowwise); if JOBVT = 'N' or 'O', VT is not referenced.

   int ldvt:    The leading dimension of the array VT.  LDVT >= 1; if JOBVT =
                'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).

   double *work: Array, dimension (LWORK). On exit, if INFO = 0, WORK(1)
                returns the optimal LWORK; if INFO > 0, WORK(2:MIN(M,N))
                contains the unconverged superdiagonal elements of an upper
                bidiagonal matrix B whose diagonal is in S (not necessarily
                sorted). B satisfies A = U * B * VT, so it has the same
                singular values as A, and singular vectors related by U and VT.

   int lwork    The dimension of the array WORK. LWORK >= 1.  LWORK >=
                MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N)-4).  For good performance,
  		LWORK should generally be larger.

   int INFO     = 0:  successful exit.
                < 0:  if INFO = -i, the i-th argument had an illegal value.
                > 0:  if DBDSQR did not converge, INFO specifies how many
                superdiagonals of an intermediate bidiagonal form B did not
                converge to zero. See the description of WORK above for detail
*/

inline void CL_DGESVD(char jobu, char jobvt, int m, int n, double *A, int lda,
               double *S, double *U, int ldu, double *VT, int ldvt, 
               double *work, int lwork, int& info)
{
  F_DGESVD(&jobu,&jobvt,&m,&n,A,&lda,S,U,&ldu,VT,&ldvt,work,&lwork,&info);
}


/*!
    void CL_DSCAL(int ntot, double coeff, double *data, int inc);

** This function scales ntot elements from array data with coeff and increment.

   inc : data*=coeff.
   
   int ntot:      length of data to scale.

   double *data:  data to scale.

   int inc:       increments for data.
*/

inline void CL_DSCAL(int size, double coeff, double *data, int inc)
{
  if( 1.!=coeff )
    F_DSCAL(&size,&coeff,data,&inc);
}


/*!
   void CL_DSYEV(char jobs, char type, int dim, double *A, int lda,
                 double *diag, double *work, int lwork, int info);
  
   This program compute all eigenvalues and, optionally, eigenvectors of a real
   symmetric matrix A.

   char jobs:      = 'N':  Compute eigenvalues only;

                   = 'V':  Compute eigenvalues and eigenvectors.

   char uplo:      = 'U':  Upper triangle of A is stored;

                   = 'L':  Lower triangle of A is stored.

   int dim:        The order of the matrix A.  N >= 0.

   double *A:      Array, dimension (LDA, N). On entry, the symmetric matrix A.
                   If uplo = 'U', the leading N-by-N upper triangular part of
                   A contains the upper triangular part of the matrix A.
                   If uplo = 'L', the leading N-by-N lower triangular part of
                   A contains the lower triangular part of the matrix A.
                   On exit, if jobs = 'V', then if INFO = 0, A contains the
                   orthonormal eigenvectors of the matrix A.  If jobs = 'N',
                   then on exit the lower triangle (if uplo='L') or the upper
                   triangle (if uplo='U') of A, including the diagonal, is
                   destroyed.

   int lda:        The leading dimension of the array A.  LDA >= max(1,N).

   double work:    Array, dimension (LWORK). On exit, if INFO = 0, WORK(1)
                   returns the optimal LWORK.

   int lwork:      The length of the array WORK.  LWORK >= max(1,3*N-1).
                   For optimal efficiency, LWORK >= (NB+2)*N, where NB is the
                   blocksize for DSYTRD returned by ILAENV.

   int info:       = 0:  successful exit

                   < 0:  if INFO = -i, the i-th argument had an illegal value

                   > 0:  if INFO = i, the algorithm failed to converge; i
  		   off-diagonal elements of an intermediate tridiagal form did
  		   not coerge to zero.
*/
void CL_DSYEV(char jobs, char type, int dim, double *matrix, int ld,
              double *diag, double *work, int lwork, int& info)
{
  F_DSYEV(&jobs, &type, &dim, matrix, &ld, diag, work, &lwork, &info);
}

/*!
  All explanation below -- in Fortran notations. Must worry about
  permuting everithing...

  void CL_DSYGV(int itype, char jobz, char uplo, int dim, double *A, int lda,
	      double *B, int ldb, double *diag, double *work, int lwork,
	      int& info);

  compute all the eigenvalues, and optionally, the eigenvectors for real
  generalized symmetric-definite  eigenproblem,  of  the form:
  A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
  Here A and B are assumed to be symmetric and B is also positive definite.

  int itype (input):  Specifies the problem type to be solved:

                 = 1:  A*x = (lambda)*B*x

		 = 2:  A*B*x = (lambda)*x

		 = 3:  B*A*x = (lambda)*x

  char jobz (input):
  
               = 'N':  Compute eigenvalues only;

	       = 'V':  Compute eigenvalues and eigenvectors.

  char uplo (input):
  
               = 'U':  Upper triangles of A and B are stored;
               = 'L':  Lower triangles of A and B are stored.

  int dim (input): The order of the matrices A and B.  dim >= 0.

  double *A (input/output):  dimension lda*dim. On entry, the symmetric
               matrix A.  If uplo =  'U', the leading dim-by-dim upper
	       triangular part of A contains the upper triangular part of
	       the  matrix  A. If uplo = 'L', the leading dim-by-dim lower
	       triangular part of A contains the lower  triangular  part  of
               the matrix A.

               On  exit,  if jobz = 'V', then if info = 0, A contains the
	       matrix Z of eigenvectors.  The eigenvectors are normalized
	       as follows:

	       if itype = 1 or 2, Z**T*B*Z = I;

	       if itype = 3, Z**T*inv(B)*Z = I.

	       If jobz  =  'N',  then on exit the upper triangle (if uplo='U')
	       or the lower triangle (if  uplo='L')  of A, including the
	       diagonal, is destroyed.

  int lda (input): The  leading  dimension  of  the  array A. lda >=max(1,dim).

  double *B (input/output):  dimension ldb*dim. On entry, the symmetric
                positive  definite  matrix B.  If uplo = 'U', the leading
		dim-by-dim upper triangular part of B contains the upper
		triangular part of the matrix B.  If uplo = 'L', the
		leading dim-by-dim lower triangular part of B  contains the
		lower triangular part of the matrix B.

		On  exit,  if  info<=dim, the part of B containing the matrix
		is overwritten by the triangular factor U  or L from the
		Cholesky factorization B = U**T*U or B = L*L**T.

  int ldb (input): The leading dimension of the array B. ldb >=  max(1,dim).

  double *ene (output): dimension dim. If INFO = 0, the eigenvalues in
                        ascending order.

  double *work (workspace/output): dimension (lwork).
               On  exit, if info = 0, work[0] returns the optimal
               lwork. 

  int lwork (input): The length of the array work. lwork >= max(1,3*dim-1).
                     For  optimal  efficiency, lwork >= (NB+2)*dim, where NB is
		     the  blocksize  for  DSYTRD returned by ILAENV.

		     If  lwork = -1, then a workspace query is assumed;
		     the routine only calculates the  optimal  size  of
		     the  work  array,  returns this value as the first
		     entry of the work array,  and  no  error  message
		     related to lwork is issued by XERBLA.

  int info (output):

               = 0:  successful exit
	       
               < 0:  if info = -i, the i-th argument had an illegal value

               > 0:  DPOTRF or DSYEV returned an error code:
	       
               <= N:  if INFO = i, DSYEV failed  to  converge;  i
               off-diagonal elements of an intermediate tridiago­
               nal form did not converge to zero; > N:   if  INFO
	       =  N  + i, for 1 <= i <= N, then the leading minor
               of order i of B is  not  positive  definite.   The
               factorization  of  B could not be completed and no
               eigenvalues or eigenvectors were computed.
*/  
void CL_DSYGV(int itype, char jobz, char uplo, int dim, double *A, int lda,
	      double *B, int ldb, double *ene, double *work, int lwork,
	      int& info)
{
  F_DSYGV(&itype,&jobz,&uplo,&dim,A,&lda,B,&ldb,ene,work,&lwork,&info);
}

/*!
     All explanation below -- in Fortran notations. Must worry about
     permuting everithing...

     void CL_DGEEV(char jobvl, char jobvr, int dim, double *A, double *wr,
                   double *wi, double *vl, double *vr,double *work, int lwork,
                   int& info);
  
   PURPOSE:
     DGEEV computes for an N-by-N real nonsymmetric matrix A, the eigenvalues
     and, optionally, the left and/or right eigenvectors.

                U^H A V=L
   
    The right eigenvector v(j) of A satisfies

  		   A * v(j) = lambda(j)	* v(j) (AV=VL)

    where lambda(j) is its eigenvalue.

    The left eigenvector u(j) of A satisfies

  		u(j)**H	* A = lambda(j)	* u(j)**H

    where u(j)**H denotes the conjugate transpose of u(j).
  
    The computed eigenvectors are normalized to have Euclidean norm equal to 1
    and largest component real.
  
   ARGUMENTS:

          char jobvl (input)
  
      	          = 'N': left eigenvectors of A	are not	computed;
	  
	          = 'V': left eigenvectors of A	are computed.
  
          char jobvr (input) 

       	          = 'N': right eigenvectors of A are not computed;

	          = 'V': right eigenvectors of A are computed.

        int dim (input)

                The order of the matrix A. N >= 0.

          double *A (input/output) matrix, dimension (dim,dim)

  	  On entry, the	N-by-N matrix A.  On exit, A has been overwritten.
  
          double *wr (output), dimension (dim)

          double *wi (output), dimension (dim)
  
          wr and wi contain the real and imaginary parts, respectively, of the
          computed eigenvalues.	Complex	conjugate pairs	of eigenvalues appear
          consecutively with the eigenvalue having the positive imaginary part
  	  first.
  
          double *vl (output), dimension (dim,dim)

  	  If JOBVL = 'V', the left eigenvectors	u(j) are stored	one after
  	  another in the columns of VL,	in the same order as their eigen-
  	  values.  If JOBVL = 'N', VL is not referenced.  If the j-th eigen- 
  	  value	is real, then u(j) = VL(:,j), the j-th column of VL.  If the
  	  j-th and (j+1)-st eigenvalues	form a complex conjugate pair, then
  	  u(j) = VL(:,j) + i*VL(:,j+1) and
  	  u(j+1) = VL(:,j) - i*VL(:,j+1).
    
          double *vr (output), dimension (dim,dim)

  	  If JOBVR = 'V', the right eigenvectors v(j) are stored one after
  	  another in the columns of VR,	in the same order as their eigen-
  	  values.  If JOBVR = 'N', VR is not referenced.  If the j-th eigen-
  	  value	is real, then v(j) = VR(:,j), the j-th column of VR.  If the
  	  j-th and (j+1)-st eigenvalues	form a complex conjugate pair, then
  	  v(j) = VR(:,j) + i*VR(:,j+1) and
  	  v(j+1) = VR(:,j) - i*VR(:,j+1).
  
          double *WORK (workspace/output), dimension (lwork)
 
 	  On exit, if INFO = 0,	WORK(1)	returns	the optimal lwork.
  
          int lwork (input) 

  	  The dimension	of the array WORK.  LWORK >= max(1,3*N), and if	JOBVL
  	  = 'V'	or JOBVR = 'V',	LWORK >= 4*N.  For good	performance, LWORK
  	  must generally be larger.
  
          int INFO (output) 
  
	  = 0:	successful exit
  
	  < 0:	if INFO	= -i, the i-th argument	had an illegal value.
  
	  > 0:	if INFO	= i, the QR algorithm failed to	compute	all the
  	  eigenvalues, and no eigenvectors have	been computed; elements	i+1:N
  	  of WR	and WI contain eigenvalues which have converged.
*/
void CL_DGEEV(char jobvl, char jobvr, int dim, double *A, double *wr, 
              double *wi, double *vl, double *vr, double *work, 
              int lwork, int& info)
{
  F_DGEEV(&jobvl,&jobvr,&dim,A,&dim,wr,wi,vl,&dim,vr,&dim,work,&lwork,&info);
}


/*!
   void CL_DSYSV(char type,int n, int nrhs, double *A, int lda, int *ipiv,
                 double *B, int ldb, double *work, int lwork, int& info);
  
   DSYSV computes the solution to a real system of linear equations  A*X =B
   where A is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.
   The diagonal pivoting method is used to factor A as
   A = U * D * U**T,  if UPLO = 'U', or A = L * D * L**T,  if UPLO = 'L',
   where U (or L) is a product of permutation and unit upper (lower) triangu-
   lar matrices, and D is symmetric and block diagonal with 1-by-1 and 2-by-2
   diagonal blocks.  The factored form of A is then used to solve the system
   of equations A * X = B.

   char uplo:   = 'U':  Upper triangle of A is stored;

                = 'L':  Lower triangle of A is stored.
  
   int N:       The number of linear equations, i.e., the order of the matrix
                A.  N >= 0.

   int nrhs:    The number of right hand sides, i.e., the number of columns of
                the matrix B.  NRHS >= 0.

   double *A:   Array, dimension (LDA,N). On entry, the symmetric matrix A.
                If UPLO = 'U', the leading N-by-N upper triangular part of A
                contains the upper triangular part of the matrix A, and the
                strictly lower triangular part of A is not referenced.
                If UPLO = 'L', the leading N-by-N lower triangular part of A
                contains the lower triangular part of the matrix A, and the
                strictly upper triangular part of A is not referenced.
                On exit, if INFO = 0, the block diagonal matrix D and the
                multipliers used to obtain the factor U or L from the
                factorization A = U*D*U**T or A = L*D*L**T as computed by
                DSYTRF.

   int lda:     The leading dimension of the array A.  LDA >= max(1,N).

   int *ipiv:   Array, dimension (N). Details of the interchanges and the
                block structure of D, as determined by DSYTRF.  If IPIV(k) > 0,
                then rows and columns k and IPIV(k) were interchanged, and
                D(k,k) is a 1-by-1 diagonal block. If UPLO = 'U' and IPIV(k) =
                IPIV(k-1) < 0, then rows and columns k-1 and -IPIV(k) were
                interchanged and D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
                If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, then rows and
                columns k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
                is a 2-by-2 diagonal block.

   double *B:   Array, dimension (LDB,NRHS). On entry, the N-by-NRHS right
                hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS
                solution matrix X.

   int ldb:     The leading dimension of the array B.  LDB >= max(1,N).

   double *work: Array, dimension (LWORK). On exit, if INFO = 0, WORK(1)
                returns the optimal LWORK.

   int lwork:   The length of WORK.  LWORK >= 1, and for best performance
                LWORK >= N*NB, where NB is the optimal blocksize for DSYTRF.
 
  int info:    = 0: successful exit

                < 0: if INFO = -i, the i-th argument had an illegal value

                > 0: if INFO = i, D(i,i) is exactly zero. The factorization has
                been completed, but the block diagonal matrix D is exactly
                singular, so the solution could not be computed.
*/
void CL_DSYSV(char type,int dim, int nrhs, double *A, int lda, int *ipiv,
              double *B, int ldb, double *work, int lwork, int& info)
{
  F_DSYSV(&type,&dim,&nrhs,A,&lda,ipiv,B,&ldb,work,&lwork,&info);
}


/*!
   void CL_DROT(int ntot, double *x, int incx, double *y, int incy,
  	        double costheta, double sintheta);

   This function calculates plane Givens rotation for vectors x,y and angle 
   theta:  x=x*cos+y*sin, y=-x*sin+y*cos; 
 
   int ntot: length of x,y
   int incx,incy: increments for x,y 
*/
void CL_DROT(int ntot, double *x, int incx, double *y, int incy,
             double costheta, double sintheta)
{
  F_DROT(&ntot,x,&incx,y,&incy,&costheta,&sintheta);
}

/*!
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO );
  
        INTEGER       INFO, LDA, LDB, N, NRHS
  
        INTEGER       IPIV( * )
  
        DOUBLE        PRECISION A( LDA, * ), B( LDB, * )
  
  PURPOSE:
    DGESV computes the solution to a real system of linear equations
       A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS
    matrices.
    The LU decomposition with partial pivoting and row interchanges is used to
    factor A as:
  
     A = P * L * U,

    where P is a permutation matrix, L is unit lower triangular, and U is upper
    triangular.  The factored form of A is then used to solve the system of
    equations A * X = B.
  
ARGUMENTS:

    N       (input) INTEGER

            The number of linear equations, i.e., the order of the matrix A.  N
            >= 0.
  
    NRHS    (input) INTEGER

            The number of right hand sides, i.e., the number of columns of the
            matrix B.  NRHS >= 0.
  
    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N).

            On entry, the N-by-N coefficient matrix A.  On exit, the factors L
            and U from the factorization A = P*L*U; the unit diagonal elements
            of L are not stored.
  
    LDA     (input) INTEGER

            The leading dimension of the array A.  LDA >= max(1,N).
  
    IPIV    (output) INTEGER array, dimension (N)

            The pivot indices that define the permutation matrix P; row i of
            the matrix was interchanged with row IPIV(i).
  
    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)

            On entry, the N-by-NRHS matrix of right hand side matrix B.  On
            exit, if INFO = 0, the N-by-NRHS solution matrix X.
  
    LDB     (input) INTEGER

            The leading dimension of the array B.  LDB >= max(1,N).
  
    INFO    (output) INTEGER

            = 0:  successful exit

            < 0:  if INFO = -i, the i-th argument had an illegal value

            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization has
            been completed, but the factor U is exactly singular, so the solu-
            tion could not be computed.
*/
void CL_DGESV(int dim, int nrhs, double *A, int lda, int *ipiv,
              double *B, int ldb, int& info)
{
  F_DGESV(&dim,&nrhs,A,&lda,ipiv,B,&ldb,&info);
}

/*!
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

      .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      ..
      .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
      ..
 
   Purpose
   =======
 
   DGETRF computes an LU factorization of a general M-by-N matrix A
   using partial pivoting with row interchanges.
 
   The factorization has the form
      A = P * L * U
   where P is a permutation matrix, L is lower triangular with unit
   diagonal elements (lower trapezoidal if m > n), and U is upper
   triangular (upper trapezoidal if m < n).
 
   This is the right-looking Level 3 BLAS version of the algorithm.
 
   Arguments
   =========
 
   M       (input) INTEGER
           The number of rows of the matrix A.  M >= 0.
 
   N       (input) INTEGER
           The number of columns of the matrix A.  N >= 0.
 
   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the M-by-N matrix to be factored.
           On exit, the factors L and U from the factorization
           A = P*L*U; the unit diagonal elements of L are not stored.
 
   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,M).
 
   IPIV    (output) INTEGER array, dimension (min(M,N))
           The pivot indices; for 1 <= i <= min(M,N), row i of the
           matrix was interchanged with row IPIV(i).
 
   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value
           > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                 has been completed, but the factor U is exactly
                 singular, and division by zero will occur if it is used
                 to solve a system of equations.
 
*/
void CL_DGETRF(int M, int N, double *A, int lda, int *ipiv,  int& info)
{
  F_DGETRF(&M,&N,A,&lda,ipiv,&info);
}

/*!
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 
   -- LAPACK routine (version 3.1) --
      Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
      November 2006
 
      .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
      ..
      .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
      ..
 
   Purpose
   =======
 
   DGETRI computes the inverse of a matrix using the LU factorization
   computed by DGETRF.
 
   This method inverts U and then computes inv(A) by solving the system
   inv(A)*L = inv(U) for inv(A).
 
   Arguments
   =========
 
   N       (input) INTEGER
           The order of the matrix A.  N >= 0.
 
   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the factors L and U from the factorization
           A = P*L*U as computed by DGETRF.
           On exit, if INFO = 0, the inverse of the original matrix A.
 
   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,N).
 
   IPIV    (input) INTEGER array, dimension (N)
           The pivot indices from DGETRF; for 1<=i<=N, row i of the
           matrix was interchanged with row IPIV(i).
 
   WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
           On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
 
   LWORK   (input) INTEGER
           The dimension of the array WORK.  LWORK >= max(1,N).
           For optimal performance LWORK >= N*NB, where NB is
           the optimal blocksize returned by ILAENV.
 
           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.
 
   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value
           > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                 singular and its inverse could not be computed.
*/

void CL_DGETRI(int N, double *A, int lda, int *ipiv, double *work, int lwork, int& info)
{
  F_DGETRI(&N,A,&lda,ipiv,work,&lwork,&info);
}


#endif

