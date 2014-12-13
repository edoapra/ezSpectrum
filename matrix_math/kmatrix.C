/*! 
\file kmatrix.C
\ingroup (MATRIX_MATH)
AIK, 1997, 2007.
*/

#include "kmatrix.h"
#include "blas_include.h"
#include "tmp_buffer.h"

//!allocates matrix
void KMatrix::Alloc(int d1, int d2)
{
#ifdef DEBUG
  check(IfAlloc(), "KMatrix::Alloc() : already allocated");
  check(!(0!=d1 && 0!=d2), "KMatrix::Alloc() : zero size");
#endif
  dim1=d1;
  dim2=d2;
  size=dim1*dim2;
  matrix=new double[size];
  check(!(matrix>NULL), "KMatrix::Alloc() : alloc failed");
#ifdef DEBUG
  if(!(matrix>NULL))
    {
      printf("KMatrix::Alloc() : alloc failed. Size requested=%d\n",size);
      exit(1);
    }
#endif
}

double KMatrix::Det(const KMatrix& a)
{
  if (a.Dim1()!=a.Dim2())
    {
      std::cout<< "\nDEBUG Error in KMatrix::Det -- determinant of nonsquare matrix is requested\n\n";
      exit(2);
    }
  int n=a.Dim1();

  int i,j,j1,j2;
  double det;
  KMatrix m;
  
  if (n == 1) 
    det = a.Elem2(0,0);
  else if (n == 2) 
    det = a.Elem2(0,0) * a.Elem2(1,1) - a.Elem2(1,0) * a.Elem2(0,1);
  else 
    {
      det = 0;
      for (j1=0;j1<n;j1++) 
	{
	  m.Adjust(n-1,n-1);
	  for (i=1;i<n;i++) 
	    {
	      j2 = 0;
	      for (j=0;j<n;j++) {
		if (j == j1)
		  continue;
		m.Elem2(i-1,j2) = a.Elem2(i,j);
		j2++;
	      }
	    }
	  det += pow(-1.0,1.0+j1+1.0) * a.Elem2(0,j1) * Det(m);
	  m.Free();
	}
    }
  return(det);
}

//For calculations change new-> MyNew, Delete-> MyDelete
//and write class MyAllocMatrix()
KMatrix::KMatrix(int d1, int d2, bool if_set) : dim1(d1), dim2(d2),
size(dim1*dim2),matrix(NULL)
{
  if(size>0)
    {
      matrix=new double[size];
      if(if_set)
	Set();
    }
}

KMatrix::KMatrix(const KMatrix& other,bool if_data) : dim1(other.dim1),
dim2(other.dim2), size(other.size), matrix(NULL)
{
  if(size>0)
    {
      matrix=new double[size];
      if(if_data)
	memcpy(matrix,other.matrix,sizeof(double)*size);
    }
}

KMatrix::~KMatrix()
{
  if(IfAlloc())
    delete [] matrix;
}

void KMatrix::Set() const
{
  memset(matrix,0,sizeof(double)*size);
}

void KMatrix::Set(double to_set) const 
{
  for(int i=0; i<size; i++)
    matrix[i]=to_set;
}

KMatrix& KMatrix::applyThreshold(const double threshold)
{
  for(int i=0; i<size; i++)
    if (fabs(matrix[i]) < threshold)
      matrix[i] = 0.0;
  return *this;
} 


void KMatrix::SetDiagonal(double d, bool if_set) const 
{
  //if if_set=true, matrix is set to be diagonal
  //otherwise diag is set to be equal d, other data unchanged
  if(size)
    { 
       int dim=(int)sqrt((double)size);
#ifdef DEBUG
       check(((size/dim)!=dim),"KMatrix::SetDiagMatrix() : not square matrix");
#endif
       if(if_set)
          Set();
  
       for(int i=0; i<dim; i++)
          matrix[i*dim+i]=d;
    }
}

void KMatrix::SetDiagonal(double d, int from, int to) const
{
  if(size)
    { 
  int dim=(int)sqrt((double)size);
#ifdef DEBUG
  check(((size/dim)!=dim),"KMatrix::SetDiagonal() : not square matrix");
  check((to-1 > dim),"KMatrix::SetDiagonal() : matrix not large enough");
#endif
  Set();

  for(int i=from; i<to; i++)
    matrix[i*dim+i]=d;
    }
}

void KMatrix::AddDiagonal(double d, int from, int to) const 
{
  if(size)
    { 
  int dim=(int)sqrt((double)size);
#ifdef DEBUG
  check(((size/dim)!=dim),"KMatrix::AddDiagonal() : not square matrix");
  check((to-1 > dim),"KMatrix::AddDiagonal() : matrix not large enough");
#endif

  for(int i=from; i<to; i++)
    matrix[i*dim+i]+=d;
    }
}

void KMatrix::GetDiagonal(KMatrix& diag) const //Extracts diagonal {
{ 
 if(size)
    { 
  int dim=(int)sqrt((double)size);
#ifdef DEBUG
  check(((size/dim)!=dim),"KMatrix::GetDiagonal() : not square matrix");
#endif
  diag.Adjust(dim,1);
  for(int i=0; i<dim; i++)
    diag[i]=matrix[i*dim+i];
    }
}

//! Calculates RMS-deviation from diagonal (with 1-s on diag) matrix 
double KMatrix::DeviationFromDiagonal() const
{
#ifdef DEBUG
  check(!(dim1==dim2),"KMatrix::DeviationFromDiagonal() : non-square matrix");
#endif
  KMatrix TMP(*this, false);
  TMP.SetDiagonal();
  return RMS(TMP); //||*this-TMP||
}

//! Check if the matrix is unitary, i.e. UU^+=1
bool KMatrix::CheckIfUnitary(char *ttl, double thresh, FILE *fl) const
{
  KMatrix TMP(*this,true);
  TMP.Transpose();
  TMP*=(*this);
  return TMP.CheckIfUnit(ttl,thresh,fl);
}




//! Check if the matrix is unit, i.e. U=1
bool KMatrix::CheckIfUnit(char *ttl, double thresh, FILE *fl) const
{
  check(!(dim1==dim2),"KMatrix::CheckIfUnit(): non-square matrix");
  bool if_ok=true;
  for(int i=0; i<dim1; i++)
    for(int j=0; j<dim2; j++)
      if((i==j && (fabs(Elem2(i,i)-1.)>thresh) ) ||
	 (i!=j && (fabs(Elem2(i,j))>thresh) ) )
	if_ok=false;

  if(!if_ok)
    Print(ttl,fl,12);
  return if_ok;
}

void KMatrix::Adjust(int d1, int d2)
{
  if( dim1!=d1 || dim2!=d2 )
    if( d1*d2!=size )
      {
	Free();
	if( 0!=d1*d2 )
	  Alloc(d1,d2);
	else//stupid case, makes sure dim-s are set right for 0-length m-s
	  {
	    dim1=d1;
	    dim2=d2;
	  }
      }
    else
      {
	dim1=d1;
	dim2=d2;
      }
}

double KMatrix::Trace() const
{
  if(0==size)
    return 0.;
  
  int dim=(int)sqrt((double)size);
#ifdef DEBUG
  check(((size/dim)!=dim),"KMatrix::Trace() : not square matrix");
#endif
  double r=0;
  for(int i=0; i<dim; i++)
    r+=matrix[i*dim+i];
  return r;
}

double KMatrix::GetAbsMaxElem() const
{
  int indx=CL_IDAMAX(size,matrix,1);
  indx--;
  return fabs(matrix[indx]);
}

double KMatrix::GetAbsMinElem() const
{
  /*
     Cannot find idamin in blas
     int indx=CL_IDAMIN(size,matrix,1);
     indx--;
     return fabs(matrix[indx]);
     */
  double minel=fabs(matrix[0]);
  for(int i=0; i<size; i++)
    minel=(fabs(matrix[i])<minel) ? fabs(matrix[i]) : minel; 
  return minel;
}

int KMatrix::GetIndxAbsMaxElem() const
{
  int indx=CL_IDAMAX(size,matrix,1);
  indx--;
  return indx;
}

int KMatrix::GetIndxAbsMinElem() const
{
  /*
     Cannot find idamin in blas
     indx=CL_IDAMIN(size,matrix,1);
     indx--;
     return fabs(matrix[indx]);
     */
  int indx;
  double minel=fabs(matrix[0]);
  for(int i=0; i<size; i++)
    if(fabs(matrix[i])<minel)
      {
	minel=fabs(matrix[i]);
	indx=i;
      }
  return indx;
}

KMatrix& KMatrix::operator=(const KMatrix& other)
{
  if(this!=&other)
    {
      Adjust(other);
      memcpy(matrix,other.matrix,sizeof(double)*size);
    }
  return *this;
}

void  KMatrix::CopyTo(KMatrix& other) const //copies this to other: other=this
{
  other.Adjust(dim1,dim2);
  memcpy(other.matrix,matrix,sizeof(double)*size);
}

//!other=coeff*this
void KMatrix::CopyScaledTo(KMatrix& other, double coeff) const 
{
  other.Adjust(dim1,dim2);
  memcpy(other.matrix,matrix,sizeof(double)*size);
  CL_DSCAL(size,coeff,other.matrix,1);
}

//!copies other to this
void KMatrix::CopyFrom(const KMatrix& other)
{
  other.CopyTo( *this);
}

void KMatrix::CopyScaledFrom(const KMatrix& other, double coeff)
{
  other.CopyScaledTo( *this, coeff);
}
  
KMatrix& KMatrix::operator+=(const KMatrix& other)
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::operator+=(): invalid operands");
#endif
  CL_DAXPY(size,1.,other.matrix,1,matrix,1);
  return *this;
}

KMatrix& KMatrix::operator-=(const KMatrix& other)
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::operator-=(): invalid operands");
#endif
  CL_DAXPY(size,-1.,other.matrix,1,matrix,1);
  return *this;
}

KMatrix& KMatrix::Add(const KMatrix& other, double coeff)
{
#ifdef DEBUG
  check(!(other.size==size),"KMatrix::Add(): Invalid dimensions");
#endif
  if(0. != coeff)
    CL_DAXPY(size,coeff,other.matrix,1,matrix,1);
  return *this;
}

KMatrix& KMatrix::AddTwo(const KMatrix& A, const KMatrix& B,double coeffa,
			 double coeffb, bool if_add)
{
#ifdef DEBUG
  //Check sizes:
  check(!(A.size==B.size),"KMatrix::AddTwo(): Invalid dimensions");
#endif

  if(!if_add)
    Adjust(A.dim1, A.dim2);
#ifdef DEBUG
  else
    check(!(A.size==size),"KMatrix::AddTwo(): Invalid dimensions");
#endif

  //y=y+coeff*x
  if(0. != coeffa)
    if(!if_add)
      {
	memcpy(matrix,A.matrix,size*sizeof(double));
	CL_DSCAL(size,coeffa,matrix,1);
      }
    else
      CL_DAXPY(size,coeffa,A.matrix,1,matrix,1);
  else
    if(!if_add)
      Set();

  CL_DAXPY(size,coeffb,B.matrix,1,matrix,1);
  return *this;
}

//!This is a right-multiply of this matrix by other matrix, A=A*B
KMatrix& KMatrix::operator*=(const KMatrix& B)
{
  int a_rows = dim1;
  int a_cols = dim2;
  int b_rows = B.dim1;
  int b_cols = B.dim2;
#ifdef DEBUG
  check(!(a_cols==b_rows),"KMatrix::operator*=(): invalid operands");
#endif

  KMatrix tmp(a_rows,b_cols,true);
  
  if(tmp.size) //call dgemm only if result has non-zero size
    /* AIK, 11/02: THIS CALL IS CORRECT! */
    CL_DGEMM('n','n',b_cols,a_rows,a_cols,1.,B.matrix,b_cols,matrix,
	     a_cols,0.,tmp.matrix,b_cols); 
  /* AIK: The call below is wrong!
    CL_DGEMM('n','n',b_cols,a_rows,a_cols,1.,B.matrix,b_rows,
    matrix, a_cols,0.,tmp.matrix,b_cols); 
  */
  *this = tmp;
  return *this;
}


//!This is a left-multiply of this matrix by other matrix.
/*! This function is the companion to *=, which is a right-multiply. A=B*A */
KMatrix& KMatrix::LeftMult(const KMatrix& B, bool if_perm)
{
  int a_rows = dim1;
  int a_cols = dim2;
  int b_rows = (!if_perm) ? B.dim1 : B.dim2; 
  int b_cols = (!if_perm) ? B.dim2 : B.dim1; 

#ifdef DEBUG
  check(!(b_cols==a_rows),"KMatrix::LeftMult: invalid dimensions operands");
#endif

  KMatrix tmp(b_rows,a_cols,true);

  if(tmp.size) //call dgemm only if result has non-zero size
    if(!if_perm) 
      CL_DGEMM('n','n',a_cols,b_rows,a_rows,1.,matrix,a_cols,
	       B.matrix,b_cols,0.,tmp.matrix,a_cols);
    else
      CL_DGEMM('n','t',a_cols,b_rows,a_rows,1.,matrix,a_cols,
	       B.matrix,b_rows,0.,tmp.matrix,a_cols);
  *this = tmp;
  return *this;
}



//! Perform a matrix multiplication C = coeff* A * B
/* Require C to have at least the minimum number of rows/cols to hold the
   result */
KMatrix& KMatrix::Multiply(const KMatrix& A, const KMatrix& B, double coeff,
			 bool if_add)
{
  int a_rows = A.dim1;  int a_cols = A.dim2;
  int b_rows = B.dim1;  int b_cols = B.dim2;
  int c_rows = dim1;    int c_cols = dim2;

#ifdef DEBUG
  check(!(a_cols==b_rows),
	"KMatrix::Multiply(): A and B mismatched dimensions");
#endif

  if(!if_add)
    {
      Adjust(a_rows,b_cols);
      c_rows = dim1;
      c_cols = dim2;
      Set();
    }
#ifdef DEBUG
  else
    {
      check((c_rows < a_rows),
	    "KMatrix::Multiply(): C matrix not enough rows");
      check((c_cols < b_cols),
	    "KMatrix::Multiply(): C matrix not enough cols");
    }
#endif
  if(size) //call dgemm only if result has non-zero size
    CL_DGEMM('n','n',b_cols,a_rows,a_cols,coeff,B.matrix,b_cols,A.matrix,
	     a_cols,(double)if_add,matrix,c_cols);

  return *this;
}


//! Perform a matrix multiplication C = coeff* A * B
/*! Can accept transposed A and B. Require C to have at least the minimum
  number of rows/cols to hold the result */
KMatrix& KMatrix::Multiply(const KMatrix& A, bool trans_a, const KMatrix& B,
			 bool trans_b, double coeff, bool if_add)
{
  int a_rows = (!trans_a) ? A.dim1 : A.dim2;
  int a_cols = (!trans_a) ? A.dim2 : A.dim1;
  int b_rows = (!trans_b) ? B.dim1 : B.dim2;
  int b_cols = (!trans_b) ? B.dim2 : B.dim1;
  int c_rows = dim1;    int c_cols = dim2;

#ifdef DEBUG
  check(!(a_cols==b_rows),
	"KMatrix::Multiply(): A and B mismatched dimensions");
#endif

  if(!if_add)
    {
      Adjust(a_rows,b_cols);
      c_rows = dim1;
      c_cols = dim2;
      Set();
    }
#ifdef DEBUG
  else
    {
      check((c_rows < a_rows),
	    "KMatrix::Multiply(): C matrix not enough rows");
      check((c_cols < b_cols),
	    "KMatrix::Multiply(): C matrix not enough cols");
    }
#endif

  char transa = (trans_a) ? 't' : 'n';
  char transb = (trans_b) ? 't' : 'n';
  
  int ldb=(!trans_b) ? b_cols : b_rows;
  int lda=(!trans_a) ? a_cols : a_rows;

  if(size) //call dgemm only if result has non-zero size  
    CL_DGEMM(transb, transa, b_cols, a_rows, a_cols, coeff, B.matrix, ldb,
	     A.matrix, lda, (double)if_add, matrix, c_cols);

  return *this;
}

KMatrix& KMatrix::Contract(const KMatrix& A, const KMatrix& B, double coeff,
			   bool if_add)   //c=coeff*a*b^+
{
  int a_rows = A.dim1;  int a_cols = A.dim2;
  int b_rows = B.dim1;  int b_cols = B.dim2;

#ifdef DEBUG
  check(a_cols!=b_cols, "KMatrix::Contract(): A and B mismatched dimensions");
#endif

  if(!if_add)
    {
      Adjust(a_rows,b_rows);
      Set();
    }
#ifdef DEBUG
  else
    check(!((dim1==a_rows)&&(dim2==b_rows)),
	  "KMatrix::Contract() : Invalid dimensions of result");
#endif

  int c_cols = dim2;

  if(size) //call dgemm only if result has non-zero size
    CL_DGEMM('t','n',b_rows,a_rows,b_cols,coeff,B.matrix,b_cols,A.matrix,
	     a_cols, (double)if_add,matrix,c_cols);

  return *this;
}

//!c=a(T)*B^+(T)
KMatrix& KMatrix::Contract(const KMatrix& A, bool trans_a, const KMatrix& B,
			 bool trans_b, double coeff, bool if_add)
{
  int a_rows = (!trans_a) ? A.dim1 : A.dim2;
  int a_cols = (!trans_a) ? A.dim2 : A.dim1;
  int b_rows = (!trans_b) ? B.dim1 : B.dim2;
  int b_cols = (!trans_b) ? B.dim2 : B.dim1;

#ifdef DEBUG
  //Check sizes:
  check(a_cols!=b_cols, "KMatrix::Contract():  A and B mismatched dimensions");
#endif

  if(!if_add)
    {
      Adjust(a_rows,b_rows);
      Set();
    }
#ifdef DEBUG
  else
    check(!((dim1==a_rows)&&(dim2==b_rows)),
	  "KMatrix::Contract() : Invalid dimensions of result");
#endif

  int c_cols = dim2;
  
  char transa = (!trans_a) ? 'n' : 't';
  char transb = (!trans_b) ? 't' : 'n';

  int ldb=(!trans_b) ? b_cols : b_rows;
  int lda=(!trans_a) ? a_cols : a_rows;

  if(size) //call dgemm only if result has non-zero size
    CL_DGEMM(transb,transa,b_rows,a_rows,b_cols,coeff,B.matrix,ldb,
	     A.matrix,lda,(double)if_add,matrix,c_cols);
  return *this;
}

//!c_i=coeff*a_i*b_i  
KMatrix& KMatrix::VectorMultiply(const KMatrix& one, const KMatrix& two,
			       double coeff,bool if_add) 
{
#ifdef DEBUG
  check(!(one.size==two.size),"KMatrix::VectorMultiply(): Invalid dimensions");
#endif
  
  if(!if_add)
    {
      Adjust(one.dim1, one.dim2);
      Set();
    }
#ifdef DEBUG
  else
    check(!(one.size==size),"KMatrix::VectorMultiply(): Invalid dimensions");
#endif

  int i;

  if(coeff==1.)
    for(i=0; i<size; i++)
      matrix[i]+=one.matrix[i]*two.matrix[i];
  else
    if(coeff==-1.)
      for(i=0; i<size; i++)
	matrix[i]-=one.matrix[i]*two.matrix[i];
    else
      for(i=0; i<size; i++)
	matrix[i]+=coeff*one.matrix[i]*two.matrix[i];
  return *this;
}

//!c_i*=coeff*other_i  
KMatrix& KMatrix::VectorMultiply(const KMatrix& other, double coeff) 
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::VectorMultiply(): Invalid dimensions");
#endif

  int i;

  if(coeff==1.)
    for(i=0; i<size; i++)
      matrix[i]*=other.matrix[i];
  else
    if(coeff==-1.)
      for(i=0; i<size; i++)
	matrix[i]*=-other.matrix[i];
    else
      for(i=0; i<size; i++)
	matrix[i]*=(coeff*other.matrix[i]);
  return *this;
}

//!c_ij=sqrt(c_ij)
/*! Assumes that all elements of matrix are non-negative */
KMatrix& KMatrix::Sqrt() 
{
  for(int i=0; i<size; i++)
    matrix[i]=sqrt(fabs(matrix[i]));
  return *this;
}


//! Transform
/*!
  Perform a similarity transform on a matrix
  C=B(T)*A*B (forwards transform) or C=B*A*B(T) (backwards transform).
 
  backtr: if true, do a backwards transform (default false)

  AIK, 09/01: I think I have fixed a bug here: the old calls of DGEMM seem
  to be wrong. Just in case I screwed up, here is what old calls looked like:
  CL_DGEMM(transb,'n',b_cols,a_rows,a_cols,1.,B.matrix, b_cols (!!),
           A.matrix,a_cols,0.,tmp,b_cols);  

  CL_DGEMM('n',transb,b_cols,b_cols,b_rows,1.,tmp, b_cols (!!),
           B.matrix,b_rows,0.,matrix,c_cols);

  11/02: GRRRRR!!!!! Spent half a day debuging this thing... Seems to
  work correctly now....
*/
KMatrix& KMatrix::Transform(const KMatrix& A,const KMatrix& B, bool backtr)
{
  int a_rows = A.dim1;  int a_cols = A.dim2;
  int c_rows = dim1;    int c_cols = dim2;
  check((a_cols != a_rows), "KMatrix::Transform(): cols of A != rows of A");
  check((c_cols != c_rows), "KMatrix::Transform(): cols of C != rows of C");
  Set(); //AIX, 06/04
  int b_rows = B.dim1; 
  int b_cols = B.dim2;

  if(true==backtr) // C=B*A*B(T) 
    {
      check((a_cols != b_cols),"KMatrix::Transform(): cols of A != rows of B");
      check((c_rows != b_rows),"KMatrix::Transform(): rows of C != cols of B");
      //tmp = AxB^T
      int tmp_rows=a_rows;
      int tmp_cols=b_rows;
      double *tmp=TmpBufferDouble::theBuffer().GetBuff(tmp_rows*tmp_rows,true);
      CL_DGEMM('t','n',b_rows,a_rows,b_cols,1.,B.matrix, b_cols,
	       A.matrix,a_cols,0.,tmp,tmp_cols); 
      //C = Bxtmp 
      CL_DGEMM('n','n',tmp_cols,b_rows,b_cols,1.,tmp,tmp_cols,
	       B.matrix,b_cols,0.,matrix,c_cols);
      TmpBufferDouble::theBuffer().FreeBuff();
    }
  else //C=B(T)*A*B
    {
      check((a_cols != b_rows),"KMatrix::Transform(): cols of A != rows of B");
      check((c_rows != b_cols),"KMatrix::Transform(): rows of C != cols of B");
      //tmp = AxB
      int tmp_rows=a_rows;
      int tmp_cols=b_cols;
      double *tmp=TmpBufferDouble::theBuffer().GetBuff(tmp_rows*tmp_cols,true);
      CL_DGEMM('n','n',b_cols,a_rows,a_cols,1.,B.matrix, b_cols,
	       A.matrix,a_cols,0.,tmp,tmp_cols);
      //C = B(T)xtmp
      CL_DGEMM('n','t',tmp_cols,b_cols,tmp_rows,1.,tmp,tmp_cols,
	       B.matrix,b_cols,0.,matrix,c_cols);
      TmpBufferDouble::theBuffer().FreeBuff();
    }
  return *this;
}

KMatrix& KMatrix::operator/=(const KMatrix& other)
{
#ifdef DEBUG
  check(!(size==other.size),
	"KMatrix::operator/=(): invalid operands");
#endif
  VectorDivide(other, 0.);
  return *this;
}

KMatrix& KMatrix::VectorDivide(const KMatrix& other, double shift)
{
#ifdef DEBUG
  check(!(size==other.size), "KMatrix::VectorDivide(): invalid operands");
#endif
  int i;
  if (shift == 0.)
    for(i=0; i<size; i++)
      matrix[i]/=other.matrix[i];
  else
    for(i=0; i<size; i++)
      matrix[i] /= (other.matrix[i] + shift);
  return *this;
}

KMatrix& KMatrix::operator*=(double coeff)
{
  Scale(coeff);
  return *this;
}

void KMatrix::Scale(double coeff) const
{
  CL_DSCAL(size,coeff,matrix,1);
}

double KMatrix::DotProduct(const KMatrix& other, short sign) const
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::DotProduct(): invalid sizes");
#endif
  double r=CL_DDOT(size,matrix,1,other.matrix,1);
  return ((sign==1) ? r : -r);
}

double KMatrix::SNorm() const
{
  return CL_DDOT(size,matrix,1,matrix,1);
}

double KMatrix::RMS(const KMatrix& other, short sign) const
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::RMS(): invalid sizes");
#endif
  double r=0., tmp;
  int i;
  if(sign>0)
    for(i=0; i<size; i++)
      {
	tmp=matrix[i]-other.matrix[i];
	r+=tmp*tmp;
      }
  else
    for(i=0; i<size; i++)
      {
	tmp=matrix[i]+other.matrix[i];
	r+=tmp*tmp;
      }
  return r;
}

double KMatrix::GetMaxRMSElem(const KMatrix& other, short sign) const
{
#ifdef DEBUG
  check(!(size==other.size),"KMatrix::RMS(): invalid sizes");
#endif
  double maxelem=0., r=0., tmp;
  int i;
  if(sign>0)
    for(i=0; i<size; i++)
      {
	tmp=matrix[i]-other.matrix[i];
	r=tmp*tmp;
	maxelem=(r > maxelem) ?  r : maxelem;
      }
  else
    for(i=0; i<size; i++)
      {
	tmp=matrix[i]+other.matrix[i];
	r=tmp*tmp;
	maxelem=(r > maxelem) ?  r : maxelem;
      }
  return maxelem;
}

//! transpose the patrix (aka permute)
KMatrix& KMatrix::Transpose()
{
  double *tmatr= TmpBufferDouble::theBuffer().GetBuff(size);
  memcpy(tmatr,matrix,sizeof(double)*size);
  SwapDim();
  int i, j, start;
  start = dim1 % 8;
  // cleanup loop
  for (i = 0; i < start; i++)
    for (j = 0; j < dim2; j++)
	matrix[i*dim2 + j] = tmatr[j*dim1 + i];
  // eight times outer loop unroll
  for (i = start; i < dim1; i = i + 8)
    for (j = 0; j < dim2; j++)
      {
	matrix[i*dim2 + j] = tmatr[j*dim1 + i];
	matrix[(i+1)*dim2 + j] = tmatr[j*dim1 + i+1];
	matrix[(i+2)*dim2 + j] = tmatr[j*dim1 + i+2];
	matrix[(i+3)*dim2 + j] = tmatr[j*dim1 + i+3];
	matrix[(i+4)*dim2 + j] = tmatr[j*dim1 + i+4];
	matrix[(i+5)*dim2 + j] = tmatr[j*dim1 + i+5];
	matrix[(i+6)*dim2 + j] = tmatr[j*dim1 + i+6];
	matrix[(i+7)*dim2 + j] = tmatr[j*dim1 + i+7];
      }
  
  /*for(i=0; i<dim1; i++) //goes by ROWS of permuted matrix
    CL_DCOPY(dim2,&tmatr[i],dim1,&matrix[i*dim2],1); */

  TmpBufferDouble::theBuffer().FreeBuff();
  return *this;
}      

//! inverse the matrix
KMatrix& KMatrix::Inverse()
{
#ifdef DEBUG
  check(!(dim1==dim2),"KMatrix::DeviationFromDiagonal() : non-square matrix");
#endif

  double *t_matr_double;
  t_matr_double=new double[size];
  memcpy(t_matr_double,matrix,sizeof(double)*size);

  int N=dim1;
  int status=0;
  
  int *t_vec_int;
  t_vec_int=new int[dim1]; //IPIV

  double *t_vec_double;
  t_vec_double=new double[dim1]; //WORK

  CL_DGETRF(dim1, dim1, t_matr_double, dim1, t_vec_int, status);
  CL_DGETRI(dim1, t_matr_double, dim1, t_vec_int, t_vec_double, dim1,status);
  check(!(status==0),"KMatrix::MATRIX INVERSION FAILED");

  memcpy(matrix, t_matr_double,sizeof(double)*size);

  delete [] t_vec_int;
  delete [] t_vec_double;
  delete [] t_matr_double;

  return *this;
}

KMatrix& KMatrix::Symmetrize(double coeff, short sgn) 
{
#ifdef DEBUG
  check(dim1!=dim2,"KMatrix::Symmetrize() : non-square matrix");
#endif
  int i,j;
  for(i=0; i<dim1; i++)
    for(j=i; j<dim1; j++)
      {
	matrix[i*dim1+j]+=(sgn>0) ? matrix[j*dim1+i] : -matrix[j*dim1+i];
	matrix[i*dim1+j]=(1.==coeff) ? matrix[i*dim1+j] :
	coeff*matrix[i*dim1+j];
	matrix[j*dim1+i]=matrix[i*dim1+j];
      }
  return *this;
}

//Filter matrix value and set all values below threshold = threshold
//if_filter_neg: true -- filter negative values as well
//Returns if any elements were filtered
bool KMatrix::Filter(double thresh, double to_repl, bool if_filter_negative,
		    bool if_max)
     const
{
  bool filtered = false;
#ifdef DEBUG
  check((thresh<0),"KMatrix::Filter() :  negative threshold");
#endif
  if(!if_max)
    if( if_filter_negative)
      for(int i=0; i<size; i++)
	{
	  if (matrix[i] < thresh)
	    {
	      matrix[i] = to_repl;
	      filtered = true;
	    }
	}
    else
      for(int i=0; i<size; i++)
	{
	  if (fabs(matrix[i]) < thresh)
	    {
	      matrix[i] = SIGN(matrix[i]) * to_repl;
	      filtered = true;
	    }
	}
  else
    if( if_filter_negative)
      for(int i=0; i<size; i++)
	{
	  if (matrix[i] > thresh)
	    {
	      matrix[i] = to_repl;
	      filtered = true;
	    }
	}
    else
      for(int i=0; i<size; i++)
	{
	  if (fabs(matrix[i]) > thresh)
	    {
	      matrix[i] = SIGN(matrix[i]) * to_repl;
	      filtered = true;
	    }
	}
  return filtered;
}

int KMatrix::NNegative() const
{
  int r=0;
  for(int i=0; i<size; i++)
    r+=(matrix[i]<0.) ? 1 : 0 ;
  return r;
}

//!Calculates deviation from hermitian matrix
double KMatrix::GetDevFromHermitian()
{
  KMatrix tmp(*this,true);
  tmp.Transpose();
  return RMS(tmp);
}

//! SVD: MATR=alpha Sigma beta^+;
/*! Overwrites matrix Sigma is dimX1 matrix, column vector.
  Returns sigma in the DESCENDING order in contrast to OSVD.
  Return alpha, beta_transp
*/
void KMatrix::SVD(KMatrix& alpha, KMatrix& beta_transp, KMatrix& sigma)
{
  if(size) //do only for non-zero size
    {
      int tdim=MAX(dim1,dim2);

      int lwork=5*tdim;//min size=5dim-4; for good performance should be larger
      double *tmpwork=TmpBufferDouble::theBuffer().GetBuff(lwork,true);
      int info;
      //AIK, 06/03, avoid UMR's
      alpha.Set();
      beta_transp.Set();
      sigma.Set();
      
      CL_DGESVD('A','A',dim2,dim1,matrix,dim2,sigma.matrix,
		beta_transp.matrix,dim2,alpha.matrix,dim1,tmpwork,lwork,info);
      check(0 != info,"SVD");
      TmpBufferDouble::theBuffer().FreeBuff();
    }
}

//! Diagonolize, only for Hermitian matrices.
/*!
  HC=CE, B=C+, A_new=BAB+.
  Puts egenvectors (B) in the original matrix.
  Ene is dimX1 matrix, column vector
*/
KMatrix& KMatrix::Herm_Diag(KMatrix& ene, bool if_reorder) 
{
  if(size) //do only for non-zero size
    {
#ifdef DEBUG
      check(dim1!=dim2,"KMatrix::Herm_Diag() : non-square matrix");
      check(!(ene.dim1==dim1 && ene.dim2==1),
	    "KMatrix::Herm_Diag() : wrong dimensions");
#endif
      
      int lwork=3*dim1-1, info;
      double *work=TmpBufferDouble::theBuffer().GetBuff(lwork,true);
      //AIK, 06/03, avoid UMR
      ene.Set();
      
      CL_DSYEV('V','U',dim1,matrix,dim1,ene.matrix,work,lwork,info);
      check(0!=info,"KMatrix::Herm_Diag() : diagonalization failed");
      TmpBufferDouble::theBuffer().FreeBuff();
      if( if_reorder ) //Reorder eigenvalues and eigenvectors
	for(int i=0; i<dim1/2; i++)
	  {
	    Swap(ene[i],ene[dim1-1-i]);
	    SwapRows(i,(dim1-1-i));
	  }
    }
  return *this;
}

/*!
  \brief Solves A*x=(lambda)*B*x for Hermitian A and positive definite B.
  
  Overwrites matrix A by eigenvectors, i.e., A(new, old) --
  eigenvectors are stored as ROWS in A. If if_reorder=false,
  eigenenergies are given in the ascending order.
*/
KMatrix& KMatrix::Herm_GenDiag(KMatrix& B, KMatrix& ene, bool if_reorder)
{
  if(size) //do only for non-zero size
    {
      check(dim1!=dim2,"KMatrix::Herm_GenDiag() : non-square matrix");
      check(!(ene.dim1==dim1 && ene.dim2==1),
	    "KMatrix::Herm_GenDiag() : wrong dimensions for ene");
      check(!(B.dim1==dim1 && B.dim2==dim2),
	    "KMatrix::Herm_GenDiag() : wrong dimensions for B");
      int itype=1; //Do A*x = (lambda)*B*x
      char jobz='V'; //Do eigenvectors 
      int lwork=3*dim1-1, info;
      double *work=TmpBufferDouble::theBuffer().GetBuff(lwork,true);
      //AIK, 06/03, avoid UMR's
      ene.Set();
      
      CL_DSYGV(itype,jobz,'U',dim1,matrix,dim1,B.matrix,dim1,ene.matrix,
	       work,lwork,info);
      check(0!=info,"KMatrix::Herm_GenDiag() : diagonalization failed");
      TmpBufferDouble::theBuffer().FreeBuff();
      if( if_reorder ) //Reorder eigenvalues and eigenvectors
	for(int i=0; i<dim1/2; i++)
	  {
	    Swap(ene[i],ene[dim1-1-i]);
	    SwapRows(i,(dim1-1-i));
	  }
    }
  return *this;
}

//!Diagonalize general matrices: solves Left^H A Right = Lambda
/*! Overwrites matrix. real_ene and img_ene
  are dimX1 matrices containing the real and imaginary parts of the
  eigenvalues.  LeftVecs and RightVecs are the left and right eigenvectors.
  Right and left eigenvectors are stored in rows.
  They can be returned biorthonormal by setting if_renorm to true.
*/
void KMatrix::Gen_Diag(KMatrix& real_ene, KMatrix& img_ene, KMatrix& LeftVecs,
		       KMatrix& RightVecs, bool if_reorder, bool if_renorm,
		       double img_zero_tol)
{
  if(size) //do only for non-zero size
    {
#ifdef DEBUG
      check(dim1!=dim2,"KMatrix::Gen_Diag() : non-square matrix");
      check(!(real_ene.dim1==dim1 && real_ene.dim2==1
	      && img_ene.dim1==dim1 && img_ene.dim2==1),
	    "KMatrix::Gen_Diag() : wrong dimensions");
#endif
      int dim=dim1, info;
      int lwork=8*dim; //Minimum is 4*dim for Left and Right vectors
      double *work=TmpBufferDouble::theBuffer().GetBuff(lwork,true);
      //AIK, 06/03, avoid UMR's
      real_ene.Set();
      img_ene.Set();
      RightVecs.Set();
      LeftVecs.Set();
      
      CL_DGEEV('V','V',dim,matrix,real_ene.matrix,img_ene.matrix,
	       RightVecs.matrix,LeftVecs.matrix,work,lwork,info);
      check(0!=info,"KMatrix::Gen_Diag() : diagonalization failed");
      TmpBufferDouble::theBuffer().FreeBuff();
      //Now Right=RightVecs^T and Left=LeftVecs^H
      //Vectors in right are stored in rows, Right(L,k) -- L-th eigenvector
      //Vectors in Left^H are stored in rows, Left(L,k) -- L-th eigenvector, cc
      //DGEEV does not orger eigenvalues. We want to reorder them
      //from small to large
      int i,j;
      if( if_reorder ) //Reorder eigenvalues and eigenvectors
	{
	  for(i=0; i<dim; i++)
	    for(j=i+1; j<dim; j++) //order from small to large
	      if(real_ene[i]>real_ene[j])
		{
		  Swap(real_ene[i],real_ene[j]);
		  Swap(img_ene[i],img_ene[j]);
		  RightVecs.SwapRows(i,j);
		  LeftVecs.SwapRows(i,j);
		}
	  //now take care about complex-conjugate pairs
	  for(i=0; i<dim-1; i++)
	    if((real_ene[i]==real_ene[i+1]) &&
	       (img_ene[i]==-1*img_ene[i+1]))
	      if(img_ene[i]<0.) //pair is swithched
		{
		  Swap(img_ene[i],img_ene[i+1]);
		  RightVecs.SwapRows(i,i+1);
		  LeftVecs.SwapRows(i,i+1);
		}
	}
      
      //now check if there are complex values -- send warning
      for(i=0; i<dim; i++)
	{  
	  //if imaginary part within img_zero_tol of zero, set to zero
	  if (fabs(img_ene[i]) < img_zero_tol)
	    img_ene[i] = 0.;
	  if(img_ene[i] > 0)
	    printf("Matrix::Gen_Diag() : eigenvalue of %.2E+%.2Ei for root %d\n",
		   real_ene[i], img_ene[i], i);
	  else if(img_ene[i] < 0)
	    printf("Matrix::Gen_Diag() : eigenvalue of %.2E%.2Ei for root %d\n",
		   real_ene[i], img_ene[i], i);
	}
      
      if( if_renorm )
	BiOrthogonalizeRows(LeftVecs,RightVecs,real_ene,img_ene);
    }
}

//!solves AX=C, puts X into C, for hermitian undefined matrix.
/*! C is dimX1 matrix, column vector */
bool KMatrix::SolveHermUndef(KMatrix& CCCP) const
{
#ifdef DEBUG
  check(dim1!=dim2,"KMatrix::SolveHermUndef() : non-square matrix");
  check(!(CCCP.dim1==dim1 && CCCP.dim2==1),
	"KMatrix::SolveHermUndef() : wrong dimensions");
#endif
  int dim=dim1;
  int info, nrhs=1,lwork=dim*dim;
  int *ipiv=TmpBufferInt::theBuffer().GetBuff(dim,true);
  double *work=TmpBufferDouble::theBuffer().GetBuff(lwork,true);
  //  CL_DSYSV('U',dim,nrhs,matrix,dim,ipiv,C.matrix,dim,work,lwork,info);

  // Ed's attempt to make this work on DECs with dxml library
  CL_DGESV(dim,nrhs,matrix,dim,ipiv,CCCP.matrix,dim,info);

  TmpBufferDouble::theBuffer().FreeBuff();
  TmpBufferInt::theBuffer().FreeBuff();
  return 0==info;
}

//!Sets data to 0 in a given row 
void KMatrix::SetRow(int i) const
{
  memset(&matrix[i*dim2],0,sizeof(double)*dim2);
}

//!Set data to 0 in a given coloumn
void KMatrix::SetColoumn(int i) const 
{
  MultColoumn(i,0.);
}

void KMatrix::SwapRows(int i, int j) const
{
  if(i!=j)
    {
#ifdef DEBUG
      check(!(i<dim1 && j<dim1),"KMatrix::SwapRows() : invalid i,j");
#endif
      double *tmp=TmpBufferDouble::theBuffer().GetBuff(dim2);
      memcpy(tmp,&matrix[i*dim2],dim2*sizeof(double));
      memcpy(&matrix[i*dim2],&matrix[j*dim2],dim2*sizeof(double));
      memcpy(&matrix[j*dim2],tmp,dim2*sizeof(double));
      TmpBufferDouble::theBuffer().FreeBuff();
    }
}

/*! Does not reallocate matrix, just move data and change dimentions.
  Later can change and de-allocate, if there will be need in it. */
void KMatrix::DeleteRow(int i)
{
#ifdef DEBUG
  check(!(i<dim1),"KMatrix::DeleteRow(int i) : invalid i");
#endif
  if( i<dim1-1 )
    memmove(&matrix[i*dim2],&matrix[(i+1)*dim2],(dim1-i-1)*dim2*
	    sizeof(double));
  dim1--;
  size=dim1*dim2;
}

void KMatrix::SwapColoumns(int i, int j) const
{
  if(i!=j)
    {
#ifdef DEBUG
      check(!(i<dim2 && j<dim2),"KMatrix::SwapColoumns() : invalid i,j");
#endif
      double *tmp=TmpBufferDouble::theBuffer().GetBuff(dim1);
      CL_DCOPY(dim1,&matrix[i],dim2,tmp,1);
      CL_DCOPY(dim1,&matrix[j],dim2,&matrix[i],dim2);
      CL_DCOPY(dim1,tmp,1,&matrix[j],dim2);
      TmpBufferDouble::theBuffer().FreeBuff();
    }
}

//!Givens rotations of rows/coloumns : x=rx*sin + ry*cos; y=rx*cos-ry*sin
void KMatrix::GivensRotRows(int rx, int ry, double costheta, double sintheta)
     const
{
  if(rx!=ry)
    {
#ifdef DEBUG
      check(!(rx<dim1 && ry<dim1),"KMatrix::GivensRotRows() : invalid i,j");
#endif
      int rowsize=dim2;
      int inc=1;
      CL_DROT(rowsize,&matrix[rx*dim2],inc,&matrix[ry*dim2],inc,costheta,
	      sintheta);
    }
}

void KMatrix::GivensRotColoumns(int rx, int ry, double costheta,
			       double sintheta) const
{
  if(rx!=ry)
    {
#ifdef DEBUG
      check(!(rx<dim2 && ry<dim2),"KMatrix::GivensRotColoumns() : invalid i,j");
#endif
      int colsize=dim1;
      int inc=dim2;
      CL_DROT(colsize,&matrix[rx],inc,&matrix[ry],inc,costheta,sintheta);
    }
}

void KMatrix::CopyRowFrom(const KMatrix& other, int i, int j) const
{
  check(!(Dim2()==other.Dim2()), "KMatrix::CopyRowFrom() : invalid Dim2");
  check(!(i<other.Dim1() && j<Dim1() ),
	"KMatrix::CopyRowFrom() : invalid rows");
  memcpy(&matrix[j*dim2],&other.matrix[i*dim2],dim2*sizeof(double));
}

void KMatrix::CopyColoumnFrom(const KMatrix& other, int i, int j) const
{
  check(!(Dim1()==other.Dim1()), "KMatrix::CopyColoumnFrom() : invalid Dim1");
  check(!(i<other.Dim2() && j<Dim2()),
	"KMatrix::CopyColoumnFrom() : invalid coloumns");
  CL_DCOPY(dim1,&other.matrix[i],dim2,&matrix[j],dim2);
}

int KMatrix::MaxValInRow(int i) const
{
  int indx=CL_IDAMAX(dim2,&matrix[i*dim2],1);
  indx--;
  return indx;
}

int KMatrix::MaxValInCol(int i) const
{
  int indx=CL_IDAMAX(dim1,&matrix[i],dim2);
  indx--;
  return indx;
}

//!max val (i,j) for j>j0
int KMatrix::MaxValInRow(int i, int j0) const 
{
  int indx=CL_IDAMAX((dim2-j0),&matrix[i*dim2+j0],1);
  indx--;
  return indx+j0;
}

//!max val (j,i) for j>j0
int KMatrix::MaxValInCol(int i, int j0) const 
{
  int indx=CL_IDAMAX((dim1-j0),&matrix[j0*dim2+i],dim2);
  indx--;
  return indx+j0;
}

//!Multiply i-th row by coeff
void KMatrix::MultRow(int i, double coeff) const
{
  CL_DSCAL(dim2,coeff,&matrix[i*dim2],1);
}

//!Multiply i-th col by coeff
void KMatrix::MultColoumn(int i, double coeff) const
{
  CL_DSCAL(dim1,coeff,&matrix[i],dim2);
}

//!Reorder ROWS of the matrix, to place largest values on diagonal
/*! For reordering rows of transformation matrix (A_new=T A T+),
  i.g., to keep correct ordering of orbitals in C matrix. */
void KMatrix::ReOrderToDiag() const
{
  int i,j;
  for(i=0; i<dim1; i++)
    {
      j=MaxValInCol(i,i);
      SwapRows(i,j);
    }
}

//!Same as ReOrderToDiag except that it also reorders other in the same way
void KMatrix::ReOrderBothToDiag(KMatrix& other)
{
  int i,j;
  for(i=0; i<dim1; i++)
    {
      j=MaxValInCol(i,i);
      SwapRows(i,j);
      other.SwapRows(i,j);
    }
}

//!Simlest way to adjust phase: set max val==positive
/*!
  There might be a problem for degenerate states, were sign is chosen randomly.
  For these cases, when two largest elements of the vector are equal, let us
  make positive one first largest element. */
void KMatrix::AdjustPhaseInRow() const
{
  int jmax, i,j;
  for(i=0; i<dim1; i++)
    {
      jmax=MaxValInRow(i);
      //Check for degeneracy:
      for(j=0; j<dim2 &&
	    (fabs(fabs(matrix[i*dim2+jmax])-fabs(matrix[i*dim2+j]))>0.000001);
	  j++);
      jmax=j; //set first of dejenerate elements, if any, to be positive
      if( matrix[i*dim2+jmax]<0. )
	MultRow(i,-1.);
    }
}

//!Same as above, but also adjusts the phase in other
void KMatrix::AdjustPhasesInBothRows(KMatrix& other)
{
  int jmax, i,j;
  for(i=0; i<dim1; i++)
    {
      jmax=MaxValInRow(i);
      //Check for degeneracy:
      for(j=0; j<dim2 &&
	    (fabs(fabs(matrix[i*dim2+jmax])-fabs(matrix[i*dim2+j]))>0.000001);
	  j++);
      jmax=j; //set first of dejenerate elements, if any, to be positive
      if( matrix[i*dim2+jmax]<0. )
	{
	  MultRow(i,-1.);
	  other.MultRow(i,-1.);
	}
    }
}

//!Add degeneracy test as above
void KMatrix::AdjustPhaseInCol() const
{
  int imax, j;
  for(j=0; j<dim2; j++)
    {
      imax=MaxValInCol(j);
      if( matrix[imax*dim2+j]<0. )
	MultColoumn(j,-1.);
    }
}

//!Scalar products of rows: Row[i]*Row[j]
double KMatrix::GetRowScalarProduct(int i, int j) const
{
#ifdef DEBUG
  check(!(i<dim1 && j<dim1),"KMatrix::GetRowScalarProduct(): invalid sizes");
#endif
  return CL_DDOT(dim2,&matrix[i*dim2],1,&matrix[j*dim2],1);
}

//!Scalar products of coloumns : ol[i]*Col[j]
double KMatrix::GetColoumnScalarProduct(int i, int j) const
{
#ifdef DEBUG
  check(!(i<dim2 && j<dim2),
	"KMatrix::GetColoumnScalarProduct():invalid sizes");
#endif
  return CL_DDOT(dim1,&matrix[i],dim2,&matrix[j],dim2);
}

//!Calculates scalar product between i-th row/coloumn of *this and j-th row/coloumn of other
double KMatrix::GetScalarProduct(const KMatrix& other, int i, bool if_i_row,
				int j, bool if_j_row) const
{
  int dim,tdim,offs1,offs2,inc1,inc2;
  dim=(if_i_row) ? dim2 : dim1;
  inc1=(if_i_row) ? 1 : dim2;
  offs1=(if_i_row) ? i*dim2 : i;
  tdim=(if_i_row) ? other.dim2 : other.dim1;
  inc2=(if_j_row) ? 1 : other.dim2;
  offs2=(if_j_row) ? j*other.dim2 : j;
#ifdef DEBUG
  check(dim!=tdim,"KMatrix::GetScalarProduct() : size does not match");
  check(!( (if_i_row && i<dim1) || (!if_i_row && i<dim2) ),
	"KMatrix::GetScalarProduct() : invalid vector");
  check(!( (if_j_row && j<other.dim1) || (!if_j_row && i<other.dim2) ),
	"KMatrix::GetScalarProduct() : invalid vector");
#endif
  return CL_DDOT(dim,&matrix[offs1],inc1,&other.matrix[offs2],inc2);
}

//! Mixes (adds to *this) i-th row/coloumn of *this and j-th row/coloumn of other
void KMatrix::MixVectors(const KMatrix& other, int i, bool if_i_row, int j,
			  bool if_j_row, double coeff) const
{
  int dim,tdim,offs1,offs2,inc1,inc2;
  dim=(if_i_row) ? dim2 : dim1;
  inc1=(if_i_row) ? 1 : dim2;
  offs1=(if_i_row) ? i*dim2 : i;
  tdim=(if_i_row) ? other.dim2 : other.dim1;
  inc2=(if_j_row) ? 1 : other.dim2;
  offs2=(if_j_row) ? j*other.dim2 : j;
#ifdef DEBUG
  check(dim!=tdim,"KMatrix::MixVectors() : size does not match");
  check(!( (if_i_row && i<dim1) || (!if_i_row && i<dim2) ),
	"KMatrix::MixVectors() : invalid vector");
  check(!( (if_j_row && j<other.dim1) || (!if_j_row && i<other.dim2) ),
	"KMatrix::MixVectors() : invalid vector");
#endif
  CL_DAXPY(dim,coeff,&other.matrix[offs2],inc2,&matrix[offs1],inc1);
}

//!Copy to i-th row/coloumn of this scaled j-th row/coloumn of other
void KMatrix::CopyScaledVectors(const KMatrix& other, int i, bool if_i_row,
			       int j, bool if_j_row, double coeff) const
{
  int dim,tdim,offs1,offs2,inc1,inc2;
  dim=(if_i_row) ? dim2 : dim1;
  inc1=(if_i_row) ? 1 : dim2;
  offs1=(if_i_row) ? i*dim2 : i;
  tdim=(if_i_row) ? other.dim2 : other.dim1;
  inc2=(if_j_row) ? 1 : other.dim2;
  offs2=(if_j_row) ? j*other.dim2 : j;
#ifdef DEBUG
  check(dim!=tdim,"KMatrix::CopyScaledVectors() : size does not match");
  check(!( (if_i_row && i<dim1) || (!if_i_row && i<dim2) ),
	"KMatrix::CopyScaledVectors() : invalid vector");
  check(!( (if_j_row && j<other.dim1) || (!if_j_row && i<other.dim2) ),
	"KMatrix::CopyScaledVectors() : invalid vector");
#endif
  CL_DCOPY(dim,&other.matrix[offs2],inc2,&matrix[offs1],inc1);
  CL_DSCAL(dim,coeff,&matrix[offs1],inc1);
}

//!Row[i]+=coeff*Row[j]
void KMatrix::MixRows(int i, int j, double coeff) const
{
#ifdef DEBUG
  check(!(i<dim1 && j<dim1),"KMatrix::MixRows: invalid sizes");
#endif
  CL_DAXPY(dim2,coeff,&matrix[j*dim2],1,&matrix[i*dim2],1);
}

//!Col[i]+=coeff*Col[j]
void KMatrix::MixColoumns(int i, int j, double coeff) const
{
#ifdef DEBUG
  check(!(i<dim2 && j<dim2),"KMatrix::MixColoumns: invalid sizes");
#endif
  CL_DAXPY(dim1,coeff,&matrix[j],dim2,&matrix[i],dim2);
}

//!Schmidt orthogonalization by rows
void KMatrix::SchmidtOrthogonRows(bool if_go_backw) const
{
  int i,j;
  double norm, sc_prod;

  if(!if_go_backw)
    for(i=0; i<dim1; i++)
      {
	for(j=0; j<i; j++)
	  {
	    sc_prod=GetRowScalarProduct(j,i); //a_j*a_i
	    MixRows(i,j,-1.*sc_prod);
	  }
	norm=GetRowScalarProduct(i,i);  //a_i*a_i
	MultRow(i,(1./sqrt(norm)));
      }
  else //do backward
    for(i=dim1-1; i>=0; i--)
      {
	for(j=i-1; j>=0; j--)
	  {
	    sc_prod=GetRowScalarProduct(j,i); //a_j*a_i
	    MixRows(i,j,-1.*sc_prod);
	  }
	norm=GetRowScalarProduct(i,i);  //a_i*a_i
	MultRow(i,(1./sqrt(norm)));
      }
}

//!Schmidt orthogonalization by coloumns
void KMatrix::SchmidtOrthogonColoumns(bool if_go_backw) const
{
  int i,j;
  double norm, sc_prod;

  if(!if_go_backw)
    for(i=0; i<dim2; i++)
      {
	for(j=0; j<i; j++)
	  {
	    sc_prod=GetColoumnScalarProduct(j,i); //a_j*a_i
	    MixColoumns(i,j,-1.*sc_prod); //a_i-=sc_prod*a_j
	  }
	norm=GetColoumnScalarProduct(i,i);  //a_i*a_i
	MultColoumn(i,(1./sqrt(norm)));
      }
  else //backward orthogonalization
    for(i=dim2-1; i>=0; i--)
      {
	for(j=i-1; j>=0; j--)
	  {
	    sc_prod=GetColoumnScalarProduct(j,i); //a_j*a_i
	    MixColoumns(i,j,-1.*sc_prod); //a_i-=sc_prod*a_j
	  }
	norm=GetColoumnScalarProduct(i,i);  //a_i*a_i
	MultColoumn(i,(1./sqrt(norm)));
      }
}

//!Check for linear dependence of rows
bool KMatrix::IfLinearDependentRows(int i, int j, double thresh) const
{
#ifdef DEBUG
  check(!(i<dim1 && j<dim1),"KMatrix::IfLinearDependentRows() : invalid i,j");
#endif
  double norm=0.,tmp;
  for(int k=0; k<dim2; k++)
    {
      tmp=matrix[i*dim2+k]-matrix[j*dim2+k];
      norm+=tmp*tmp;
    }
  return (norm<thresh);
}

bool KMatrix::IfLinearDependentRows(char *ttl, FILE *out, double thresh) const
{
  AdjustPhaseInRow();  //maximum element positive
  int i,j;
  bool r=false;
  for(i=0; i<dim1; i++)
    for(j=i+1; j<dim1; j++)
      if(IfLinearDependentRows(i,j,thresh) )
	{
	  r=true;
	  fprintf(out,ttl);
	  fprintf(out,"Rows %d,%d are lineary dependent\n",i,j);
	}
  return r;
}

//! Determinant (LAPACK version)
/*
The LU factorization yields three matrices, the product of which is the original
complex matrix. Therefore the determinat is the product of the three determinants
of P, L and U. The determinant of the triangular matrices L and U is the product
of the elements on the diagonal - as for any triangular matrix (for L this is 1
as all elements of the diagonal are one.) The determinant of P is either +1 or -1
depending of whether the number of row permutations is even or odd.
CL_DGETRF returns one matrix (L&U combined); the unit diagonal elements of L are not stored.
 */
double KMatrix::Determinant() const
{
#ifdef DEBUG
  check(!(dim1==dim2),"KMatrix::DeviationFromDiagonal() : non-square matrix");
#endif

  double *t_matr_double;
  t_matr_double=new double[size];
  memcpy(t_matr_double,matrix,sizeof(double)*size);

  int N=dim1;
  int status=0;
  
  int *t_vec_int;
  t_vec_int=new int[dim1]; 

  double *t_vec_double;
  t_vec_double=new double[dim1];

  CL_DGETRF(dim1, dim1, t_matr_double, dim1, t_vec_int, status);
  check(!(status==0),"KMatrix::Determinant() FAILED");

  double det=1;
  for(int i=0; i<N; i++)
    det*=t_matr_double[i*N+i];


  delete [] t_vec_int;
  delete [] t_vec_double;
  delete [] t_matr_double;


  return det;
}


//! Determinant -- non-LAPACK; interface for Det() [recursive function to calculate the determinant of a square function]
double KMatrix::Determinant_nonLAPACK()
{
  return Det(*this);
}

/*
//! Read block with offsets and write with offsets 
void KMatrix::ReadFromMatrix(const KMatrix& read_from, const Index& off_put,
			    const Index& off_take, const Index& nelem,
			    bool if_add) const
{
  Index tmp(off_put);
  tmp+=nelem;
#ifdef DEBUG
  check(!(tmp[0]<=dim1 && tmp[1]<=dim2),
	"KMatrix::ReadFromMatrix() : invalid offsets to put");
#endif
  tmp=off_take;
  tmp+=nelem;
#ifdef DEBUG
  check(!(tmp[0]<=read_from.dim1 && tmp[1]<=read_from.dim2),
	"KMatrix::ReadFromMatrix() : invalid offsets to take");
#endif

  double coeff=1.*if_add;
  int put_size=nelem.Elem(1),i,j;
  for(i=off_put.Elem(0),j=off_take.Elem(0);
      i<off_put.Elem(0)+nelem.Elem(0); i++,j++)
    CL_DAXPY(put_size,coeff,
             &read_from.matrix[j*read_from.dim2+off_take.Elem(1)],1,
             &matrix[i*dim2+off_put.Elem(1)],1);
}
*/
/*
void KMatrix::Print(char *name, FILE *fl, int ndigits) const
{
  static int nchar_inline=80;
  static int nchar_int=5;
 
  int nchar_per_numb=nchar_int+1+ndigits;

  int ncolomns=(nchar_inline-nchar_int)/nchar_per_numb;
  int npages=dim2/ncolomns;
  bool if_extra_page=(0!=dim2%ncolomns);
  int col_in_extra=dim2%ncolomns;

  char format_numb[20];
  char format_row_numb[20];
  char format_el[20];
  sprintf(format_numb, "%%%dd",nchar_per_numb);
  sprintf(format_row_numb, "%%%dd",nchar_int);
  sprintf(format_el, "%%%d.%dlf",nchar_per_numb,ndigits);

  fprintf(fl,"%s\n",name);
  
  int i,j,k,l;

  for(i=0; i<npages; i++)
    {
      for(k=0; k<nchar_int; k++)
	fprintf(fl," "); //space for row numbers
      for(j=0; j<ncolomns; j++)
	fprintf(fl,format_numb,i*ncolomns+j);
      fprintf(fl,"\n");
      for(l=0; l< dim1; l++)
	{
	  fprintf(fl,format_row_numb,l);
	  for(j=0; j<ncolomns; j++)
	    fprintf(fl,format_el,matrix[l*dim2+(i*ncolomns+j)]);
	  fprintf(fl,"\n");
	}
      fprintf(fl,"\n");
    }
  if(if_extra_page)
    {
      for(k=0; k<nchar_int; k++)
	fprintf(fl," "); //space for row numbers
      for(j=0; j<col_in_extra; j++)
	fprintf(fl,format_numb,i*ncolomns+j);
      fprintf(fl,"\n");
      for(l=0; l< dim1; l++)
	{
	  fprintf(fl,format_row_numb,l);
	  for(j=0; j<col_in_extra; j++)
	    fprintf(fl,format_el,matrix[l*dim2+(i*ncolomns+j)]);
	  fprintf(fl,"\n");
	}
      fprintf(fl,"\n");
    }
  fprintf(fl,"\n");
  fflush(fl);
}
*/
void KMatrix::Print(char *name, FILE *fl, int ndigits) const
{
  static int nchar_inline=80;
  static int nchar_int=5;
 
  int nchar_per_numb=nchar_int+1+ndigits;

  int ncolomns=(nchar_inline-nchar_int)/nchar_per_numb;
  int npages=dim2/ncolomns;
  bool if_extra_page=(0!=dim2%ncolomns);
  int col_in_extra=dim2%ncolomns;

  char format_numb[20];
  char format_row_numb[20];
  char format_el[20];
  sprintf(format_numb, "%%%dd",nchar_per_numb);
  sprintf(format_row_numb, "%%%dd",nchar_int);
  sprintf(format_el, "%%%d.%dlf",nchar_per_numb,ndigits);

  fprintf(fl,"%s\n",name);
  
  int i,j,k,l;

  for(i=0; i<npages; i++)
    {
      for(k=0; k<nchar_int; k++)
	fprintf(fl," "); //space for row numbers
      for(j=0; j<ncolomns; j++)
	fprintf(fl,format_numb,i*ncolomns+j);
      fprintf(fl,"\n");
      for(l=0; l< dim1; l++)
	{
	  fprintf(fl,format_row_numb,l);
	  for(j=0; j<ncolomns; j++)
	    {
	      fprintf(fl,format_el,matrix[l*dim2+(i*ncolomns+j)]);
	      //fprintf(fl,":");
	    }
	  fprintf(fl,"\n");
	}
      fprintf(fl,"\n");
    }
  if(if_extra_page)
    {
      for(k=0; k<nchar_int; k++)
	fprintf(fl," "); //space for row numbers
      for(j=0; j<col_in_extra; j++)
	fprintf(fl,format_numb,i*ncolomns+j);
      fprintf(fl,"\n");
      for(l=0; l< dim1; l++)
	{
	  fprintf(fl,format_row_numb,l);
	  for(j=0; j<col_in_extra; j++)
	    {
	      fprintf(fl,format_el,matrix[l*dim2+(i*ncolomns+j)]);
	      //fprintf(fl,":");
	    }
	  fprintf(fl,"\n");
	}
      fprintf(fl,"\n");
    }
  fprintf(fl,"\n");
  fflush(fl);
}

//make it real std:: -vm
void KMatrix::PrintScientific(char *name) const
{
  int ndigits=6;

  static int nchar_inline=80;
  static int nchar_int=5;
 
  int nchar_per_numb=nchar_int+1+ndigits;

  static int nchar_precis=4;

  int ncolomns=(nchar_inline-nchar_int)/nchar_per_numb;
  int npages=dim2/ncolomns;
  bool if_extra_page=(0!=dim2%ncolomns);
  int col_in_extra=dim2%ncolomns;

  std::cout << name << '\n';
  
  int i,j,k,l;

  for(i=0; i<npages; i++)
    {
      for(k=0; k<nchar_int; k++)
	std::cout << ' ';
      for(j=0; j<ncolomns; j++)
	std::cout << std::fixed << std::setw(nchar_per_numb) << i*ncolomns+j;
      std::cout << '\n';
      for(l=0; l< dim1; l++)
	{
	  std::cout << std::fixed << std::setw(nchar_int) << l;
	  for(j=0; j<ncolomns; j++)
	    std::cout << std::scientific << std::setprecision(nchar_precis) << std::setw(nchar_per_numb) << matrix[l*dim2+(i*ncolomns+j)];
	  std::cout << '\n';
	}
      std::cout << '\n';
    }
  if(if_extra_page)
    {
      for(k=0; k<nchar_int; k++)
	std::cout << ' ';
      for(j=0; j<col_in_extra; j++)
	std::cout << std::fixed << std::setw(nchar_per_numb) << i*ncolomns+j;
      std::cout << '\n';
      for(l=0; l< dim1; l++)
	{
	  std::cout << std::fixed << std::setw(nchar_int) << l;
	  for(j=0; j<col_in_extra; j++)
	    std::cout << std::scientific << std::setprecision(nchar_precis) << std::setw(nchar_per_numb) << matrix[l*dim2+(i*ncolomns+j)];
	  std::cout << '\n';
	}
      std::cout << '\n';
    }
  std::cout << '\n';
  std::cout << std::flush;
}



void KMatrix::Read(FILE *fl)
{
  fread(matrix,sizeof(double),size,fl);
  fflush(fl);
}

void KMatrix::Write(FILE *fl) const
{
  fwrite(matrix,sizeof(double),size,fl);
  fflush(fl);
}

double DotProduct(const KMatrix& one, const KMatrix& two, short sign)
{
  return one.DotProduct(two,sign);
}

double RMS(const KMatrix& one, const KMatrix& two, short sign)
{
  return one.RMS(two,sign);
}

//!Calculates C[i,j]=coeff*A[i,k]*B[j,k],
/*! if (if_add==true, calculate C[i,j]+=coeff*A[i,k]*B[j,k] */
void Contract2(KMatrix& res, const KMatrix& one, const KMatrix& two,
	       double coeff, bool if_add)
{
  res.Contract(one,two,coeff,if_add);
}

void Multiply2(KMatrix& res, const KMatrix& one, const KMatrix& two,
		      double coeff, bool if_add)
{
  res.Multiply(one,two,coeff,if_add);
}

//! Calculates res[i,j](+)=(coeff1*)one[i,j]*two[i,j]
void VectorMultiply2(KMatrix& res, const KMatrix& one, const KMatrix& two,
		     double coeff, bool if_add)
{
  res.VectorMultiply(one,two,coeff,if_add);
}

//! res+=coeff*one
void Add2(KMatrix& res, KMatrix& one, double coeff)
{
  res.Add(one,coeff);
}

//! Calculates res[i,j]=coeff1*one[i,j]+coeff2*two[i,j]
/* if if_add==true, calculate res[i,j]+=coeff1*one[i,j]+coeff2*two[i,j] */
void AddTwo2(KMatrix& res, KMatrix& one, KMatrix& two, double coeff1,
	    double coeff2, bool if_add)
{
  res.AddTwo(one,two,coeff1,coeff2,if_add);
}

// These operators can be efficient if compiler supports
//"Return Value Optimization"
KMatrix operator+(const KMatrix& A, const KMatrix& B)
{
  return KMatrix(A,true)+=B;
}

KMatrix operator-(const KMatrix& A, const KMatrix& B)
{
  return KMatrix(A,true)-=B;
}

KMatrix operator*(const KMatrix& A, const KMatrix& B)
{
  return KMatrix(A,true)*=B;
}

KMatrix operator*(double coeff, const KMatrix& A)
{
  return KMatrix(A,true)*=coeff;
}

KMatrix operator*(const KMatrix& A, double coeff)
{
  return KMatrix(A,true)*=coeff;
}

/*! DGEEV gives Left and Riht vectors satisfying HR_i=E_i R_i, and
   L^+_iH=L^+ i E_i. For nondegenerate states this means that
   L_j=L_i=coeff*\delta_ij. However, for degenerate states
   situation is more complicated and requires extra steps. */
void BiOrthogonalizeRows(KMatrix& LeftVecs, KMatrix& RightVecs,
			 KMatrix& re_ene, KMatrix& img_ene)
{
#ifdef DEBUG
  check( !(LeftVecs.dim1==RightVecs.dim1 && RightVecs.dim2==LeftVecs.dim2 &&
	   LeftVecs.dim1==LeftVecs.dim2),
	 "BiOrthogonalizeRows() : invalid vectors");
#endif
  int i,j,k,l, dim=LeftVecs.dim1;

  //Before transformation, adjust phase in row:
  //This potentially may cause problems for complex roots
  //LeftVecs.AdjustPhaseInRow();
  //RightVecs.AdjustPhaseInRow();
  //Check for linear dependence
  /*
  LeftVecs.IfLinearDependentRows("BiOrthogonalizeRows(), LeftVecs: \n",stdout,
				 0.00001);
  RightVecs.IfLinearDependentRows("BiOrthogonalizeRows(), RightVecs: \n",
				  stdout,0.00001);
  */
  //Now bi-orthogonolize left vectors to right vectors
  //Left^H*Right=1
  //For non-degenerate vectors we have just to normalize L_i^+R_i=1,
  //since for nondegenerate states these vectors are already orthogonal
  //For degenerate blocks we have to perform transformation in degenerate block
  // S_ij=L^+_j R_i=alpha Sigma beta^+
  //L=L*alpha, R=beta^+R, R_i,L_i/=sigma_i

  double overlap;
  /*
    //Normalize first right vectors to unity,
    //then normalize left vectors to right
    for(i=0; i<dim; i++)
    if(img_ene[i]==0.)
    {
    overlap=RightVecs.GetScalarProduct(RightVecs,i,true,i,true);
    RightVecs.MultRow(i,1./sqrt(overlap));
    }
  */
  bool *if_done = new bool[dim];
  for(i = 0; i < dim; i++)
    if_done[i] = false;
  int *degen = new int[dim];

  for(i=0; i<dim; i++)
    {
      if (if_done[i])
	continue;
      int ndeg_states=1;
      /* LVS changed how to manage complex roots. Now hey are treated as
	 degenerate roots. This corresponds to ... */
      //if( img_ene[i]==0. ) for real vectors only
      //{
      //determine if there is degeneracy for this root:
      for(j=i+1; j<dim; j++)
	if (fabs(re_ene[i] - re_ene[j]) < 0.000001)
	  {
	    if_done[j] = true;
	    degen[0] = i;
	    degen[ndeg_states++] = j;
	  }
      //printf("i=%d, ndeg_states=%d\n",i, ndeg_states);
      if(1==ndeg_states) //non-degenerate
	{
	  overlap=LeftVecs.GetScalarProduct(RightVecs,i,true,i,true);
	  RightVecs.MultRow(i,1./sqrt(fabs(overlap)));
	  LeftVecs.MultRow(i,SIGN(overlap)/sqrt(fabs(overlap)));
	}
      else //there is a block of degenerate states
	{
	  KMatrix DOV(ndeg_states,ndeg_states,true);
	  //Calculate overlapp:
	  for(k=0; k<ndeg_states; k++)
	    for(l=0; l<ndeg_states; l++)
	      DOV.Elem2(k,l)=  //L^+_k*R_l
		LeftVecs.GetScalarProduct(RightVecs,degen[k],true,
					  degen[l],true);
	  
	  KMatrix alpha(ndeg_states,ndeg_states,true);
	  KMatrix beta_plus(ndeg_states,ndeg_states,true);
	  KMatrix sigma(ndeg_states,1,true);
	  
	  DOV.SVD(alpha,beta_plus,sigma);
	  KMatrix tmp_vecs(ndeg_states,dim,true);
	  //Transform Left vectors: L^+= alpha^+ L^+,
	  //Lp_ij=\sum_k alpha_ki Lp_kj
	  for(k=0; k<ndeg_states; k++)
	    for(l=0; l<dim; l++)
	      for(j=0; j<ndeg_states; j++)
		tmp_vecs.Elem2(k,l)+=alpha.Elem2(j,k)
		  *LeftVecs.Elem2(degen[j],l);
	  //Now copy left vectors back
	  for(k=0; k<ndeg_states; k++)
	    for(l=0; l<dim; l++)
	      LeftVecs.Elem2(degen[k],l)=tmp_vecs.Elem2(k,l);
	  //Do the same for right vectors:
	  tmp_vecs.Set();
	  //Transform Right vectors:R=R beta,R^+=beta_p R^+,R_p=beta_p R_p
	  //Rp_ij=\sum_k beta_p_ik Rp_kj
	  for(k=0; k<ndeg_states; k++)
	    for(l=0; l<dim; l++)
	      for(j=0; j<ndeg_states; j++)
		tmp_vecs.Elem2(k,l)+=beta_plus.Elem2(k,j)*
		  RightVecs.Elem2(degen[j],l);
	  //Now copy left vectors back
	  for(k=0; k<ndeg_states; k++)
	    for(l=0; l<dim; l++)
	      RightVecs.Elem2(degen[k],l)=tmp_vecs.Elem2(k,l);
	  //Now normalize Left and Right vectors:
	  for(k=0; k<ndeg_states; k++)
	    {
	      RightVecs.MultRow(degen[k],1./sqrt(sigma[k]));
	      LeftVecs.MultRow(degen[k],1./sqrt(sigma[k]));
	    }
	}
      
      /* Old rocedure for handling complex roots
      //}
      //else 
      //{
      //determine if there is degeneracy for this root:
      double real = LeftVecs.GetScalarProduct(RightVecs,i,true,i,true);
      double imaginary =
      LeftVecs.GetScalarProduct(RightVecs,i+1,true,i+1,true);
      overlap =sqrt(fabs(real*real + imaginary*imaginary));
      
      LeftVecs.MultRow(i,1./sqrt(fabs(overlap)));
      RightVecs.MultRow(i,SIGN(overlap)/sqrt(fabs(overlap)));
      LeftVecs.MultRow(i+1,1./sqrt(fabs(overlap)));
      RightVecs.MultRow(i+1,SIGN(overlap)/sqrt(fabs(overlap)));
      */
    }

  /*
    //Debug check
    KMatrix OV(LeftVecs,false);
    OV.Contract(LeftVecs,false,RightVecs,false);
    
    bool if_ok=true;
    for(i=0; i<dim && if_ok; i++)
    for(j=0; j<dim && if_ok; j++)
    if(((i!=j) && (fabs(OV.Elem2(i,j))> 0.0000001)) ||
    ((i==j) && (fabs(OV.Elem2(i,j)-1.) > 0.0000001)))
    if_ok=false;
    if(!if_ok)
    {
    OV.Print("BiOrthogonalize(): Overlap matrix");
    LeftVecs.Print("BiOrthogonalize(): LeftVecs",stdout,3);
    RightVecs.Print("BiOrthogonalize(): RightVecs",stdout,3);
    }
  */

  delete [] if_done;
  delete [] degen;
}

/*!
  This function calculates N-nvecs new vectors biorthoganal to nvecs vectors
  provided, e.i., biorthogonalize N unit vectors to the nvecs vectors given.
  Matrices Left(Right)Vectors are N-row matrices and contain nvecs vectors
  given, and all the rest (N-nvecs) vectors will be overwritten.
*/
void BiOrthogonalizeToRows(int nvecs, KMatrix& LeftVecs, KMatrix& RightVecs)
{
#ifdef DEBUG
  check( !(LeftVecs.dim1==RightVecs.dim1 && RightVecs.dim2==LeftVecs.dim2 &&
	   LeftVecs.dim1==LeftVecs.dim2),
	 "BiOrthogonalizeToRows() : invalid vectors");
#endif
  int dim1=LeftVecs.dim1;
  int dim2=LeftVecs.dim2;
  int vect_tofind=dim1-nvecs;
#ifdef DEBUG
  check(!(nvecs<dim1),"BiOrthogonalizeToRows() : invalid nvecs/size");
#endif

  int nfound=0; //how much new vectors are already found
  int current=0; //what vector is generated now
  double sc_prod, sc_prod_rl, sc_prod_lr;
  int i,j;
  for(i=0,current=nvecs; i<dim2 && nfound<vect_tofind; i++)
    {
      LeftVecs.SetRow(current);
      RightVecs.SetRow(current);
      LeftVecs.Elem2(current,i)=1.;  //Start from (1,0,...0)
      RightVecs.Elem2(current,i)=1.; //Start from (1,0,...0)
      for(j=0; j<current; j++) //orthogonolize to all previous vectors
	{
	  sc_prod_lr=LeftVecs.GetScalarProduct(RightVecs,current,true,j,true);
	  sc_prod_rl=RightVecs.GetScalarProduct(LeftVecs,current,true,j,true);
	  LeftVecs.MixRows(current,j,-1.*sc_prod_lr);
	  RightVecs.MixRows(current,j,-1.*sc_prod_rl);
	}
      sc_prod=LeftVecs.GetScalarProduct(RightVecs,current,true,current,true);
      if(fabs(sc_prod)> 0.0000001)
	{
	  RightVecs.MultRow(current,1./sqrt(fabs(sc_prod)));
	  LeftVecs.MultRow(current,SIGN(sc_prod)/sqrt(fabs(sc_prod)));
	  current++;
	  nfound++;
	}
    }
  
  check(!(nfound==vect_tofind),
	"BiOrthogonalizeToRows: BiOrthogonalization failed");
  
  /*
  //Debug check
  KMatrix OV(LeftVecs,false);
  OV.Contract(LeftVecs,false,RightVecs,false);
    
  bool if_ok=true;
  int dim=dim1;
  for(i=0; i<dim && if_ok; i++)
    for(j=0; j<dim && if_ok; j++)
      if(((i!=j) && (fabs(OV.Elem2(i,j))> 0.0000001)) ||
	 ((i==j) && (fabs(OV.Elem2(i,j)-1.) > 0.0000001)))
	if_ok=false;
  if(!if_ok)
    {
      OV.Print("BiOrthogonalize(): Overlapp matrix");
      LeftVecs.Print("BiOrthogonalize(): LeftVecs",stdout,3);
      RightVecs.Print("BiOrthogonalize(): RightVecs",stdout,3);
    }
  */
}




