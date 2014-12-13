/*! \file kmatrix.h 
\ingroup MATRIX_MATH
\brief This file contains our KMatrix class declaration.

AIK, 1997.
AIK, 2007: cleaned up tensor-related and obsolete stuff

TODO: 1. repalce bool by bool, TRUE by true
      2. replace c i/o by c++ 
      3. May be add sorting functions 
*/

#ifndef _kmatrix_h
#define _kmatrix_h

#include "genincludes.h"
#include "mathutil.h"
#include "genutil.h"
#include "genincludes.h"

//#include "index.h"


//! KMatrix is a rectangular dim1Xdim2 matrices.
class KMatrix 
{
  //!# of rows
  int dim1;            
  //!# of columns
  int dim2;
  //!# total size of data: dim1xdim2
  int size; 
  //! data
  double *matrix;      
  //! Swap dimensions: dim1=dim2, dim2=dim1
  void SwapDim() { Swap(dim1,dim2); }
  //! Adjust itself to be the same as other
  void Adjust(const KMatrix& other)  
    {
      Adjust(other.dim1,other.dim2);
    }
  //! Allocates matrix
  void Alloc(int dim1, int dim2); 
  bool IfAlloc() const { return size!=0; }

  //! recursive function to calculate the determinant of a square function
  double Det(const KMatrix& a);
  
public: 
  //! Adjust itself to new dim1,dim2
  void Adjust(int dim1,int dim2);    
  
  void Free()
  {
    if(IfAlloc())
      {
	delete [] matrix;
	dim1=dim2=size=0;
	matrix=NULL;
      }
  }
  
 public:
  //! Empty constractor, no allocation
  /*! Need it for tensors and tensors lists */
  KMatrix(): dim1(0), dim2(0), size(0), matrix(NULL) {}
  KMatrix(int dim1, int dim2,bool if_set=FALSE);
  KMatrix(const KMatrix& other,bool if_copy_data=TRUE);
  ~KMatrix();
  
  int Size() const { return size;}
  int Dim1() const { return dim1;}
  int Dim2() const { return dim2;}
  
  double& Elem2(int i, int j) const
    {
#ifdef DEBUG
      check( (i>=dim1 || j>=dim2),"KMatrix::Elem2() : invalid i,j");
#endif
      return matrix[i*dim2+j];
    }

  double& Elem1(int i) const
    {
#ifdef DEBUG
      check((i>=size),"KMatrix::Elem1() : invalid i");
#endif
      return matrix[i];
    }
  
  double& operator[](int i) const
    {
#ifdef DEBUG
      check((i>=size),"KMatrix::operator[] : invalid index");
#endif
      return matrix[i];
    }
  //! Returns pointer to the data array, use with caution.
  double *GetMatrixPtr() const {return matrix;}
  //! Calculates trace of the matrix (if it is square m-x).
  double Trace() const;
  //! Finds a value of the element whose absolute value is the largest.
  double GetAbsMaxElem() const;
  //! Finds a value of the element whose absolute value is the smallest.
  double GetAbsMinElem() const;
  //! Finds an address of  the element whose absolute value is the largest.
  int GetIndxAbsMaxElem() const;
  //! Finds an address of  the element whose absolute value is the smallest.
  int GetIndxAbsMinElem() const;

  KMatrix& operator=(const KMatrix& other);
  //! Copies this to other: other=this
  void CopyTo(KMatrix& other) const;
  //! other=coeff*this
  void CopyScaledTo(KMatrix& other, double coeff) const; 
  //! copies other to this
  void CopyFrom(const KMatrix& other); 
  void CopyScaledFrom(const KMatrix& other, double coeff); 
  KMatrix& operator+=(const KMatrix& other);
  KMatrix& operator-=(const KMatrix& other);
  //! *this+=coeff*other
  KMatrix& Add(const KMatrix& other, double coeff=1.);
  KMatrix& AddTwo(const KMatrix& A, const KMatrix& B,double coeffa=1.,
		 double coeffb=1., bool if_add=FALSE);
  //! Regular matrix mult (on rt)
  KMatrix& operator*=(const KMatrix& other);   
  //! Regular matrix mult (on lt): A=other*A (if_perm==TRUE) use other^T
  KMatrix& LeftMult(const KMatrix& other, bool if_perm=FALSE);
  //! c=coeff*a*b
  KMatrix& Multiply(const KMatrix& A, const KMatrix& B, double coeff=1.,
		   bool if_add=FALSE);
  //! c=coeff*a(T)*b(T)
  KMatrix& Multiply(const KMatrix& A, bool trans_a, const KMatrix& B,
                   bool trans_b, double coeff=1., bool if_add=FALSE);  
  //! c=coeff*a*b^+
  KMatrix& Contract(const KMatrix& A, const KMatrix& B, double coeff=1.,
		   bool if_add=FALSE);
  //! c=a(T)*B^+(T)
  KMatrix& Contract(const KMatrix& A, bool trans_a, const KMatrix& B,
                   bool trans_b, double coeff=1., bool if_add=FALSE);
  //! c_ij=coeff*a_ij*b_ij
  KMatrix& VectorMultiply(const KMatrix& one, const KMatrix& two,
                   double coeff=1., bool if_add=FALSE); 
  //! c_i*=coeff*other_i
  KMatrix& VectorMultiply(const KMatrix& other, double coeff=1.);
  //! c_ij=sqrt(c_ij)
  KMatrix& Sqrt(); 
  
  //! c=b^T*a*b
  KMatrix& Transform(const KMatrix& A,const KMatrix& B, bool backtr=FALSE); 
  KMatrix& operator/=(const KMatrix& other);
  KMatrix& VectorDivide(const KMatrix& other, double shift=0.);
  KMatrix& operator*=(double);
  void Scale(double coeff) const; 
  double DotProduct(const KMatrix& other, short sign=1) const;
  double SNorm() const;

  //! Calculates \sum_i (this[i]-sign*other[i])^2.
  double RMS(const KMatrix& other, short sign=1) const;
  //! Finds MAX value of the (this[i]-sign*other[i])^2 vector.
  double GetMaxRMSElem(const KMatrix& other, short sign=1) const;
  //! transpose the patrix 
  KMatrix& Transpose();
  //! Inverse matrix
  KMatrix& Inverse();
  //! res[i,j]=coeff*(matrix[i,j]+sgn*matrix[j,i])
  KMatrix& Symmetrize(double coeff=1.,short sgn=1);
  //! matrix[i]=0
  void Set() const; 
  //! if element<threshold, then set it to 0
  KMatrix& applyThreshold(const double threshold); 
  //! matrix[i]=to_set
  void Set(double to_set) const; 
  void SetDiagonal(double d=1., bool if_set=TRUE) const;
  void SetDiagonal(double d, int from, int to) const;
  void AddDiagonal(double d, int from, int to) const;
  void SetDiagonal(double d, int n) const ;
  //! Extracts diagonal 
  void GetDiagonal(KMatrix& diag) const; 

  /*! 
   \brief Filter matrix value and set all values to a trheshold value.
  
   Sets it below threshold (if_filter_max=FALSE)
   or above threshold (if_filter_max=TRUE) = to_repl
   if_filter_neg: TRUE -- filter negative values as well.
   Returns if any elements were filtered. */
  bool Filter(double thresh, double to_repl, bool if_filter_neg=TRUE,
	      bool if_filter_max=FALSE) const;
  //! Calculate number of negative elements.
  int NNegative() const;
  //! Calculates RMS(A,A.Permute());
  double GetDevFromHermitian();
  //! Calculates RMS-deviation from unit diagonal matrix 
  double DeviationFromDiagonal() const;

  //! Check if the matrix is unitary, i.e. UU^+=1
  bool CheckIfUnitary(char *ttl, double thresh=0.000001, FILE *fl=stdout)
    const;
  //! Check if the matrix is unit, i.e. U=1
  bool CheckIfUnit(char *ttl, double thresh=0.000001, FILE *fl=stdout) const;
  //! Single Value Decomposition u=alpha D beta^+, overwrites matrix.
  void SVD(KMatrix& alpha, KMatrix& beta_c, KMatrix& d);
  //! Diagonolize, only for Hermitian matrices. Returns eigenvectors.
  KMatrix& Herm_Diag(KMatrix& ene, bool if_reorder); 
  /*!
    \brief Solves A*x=(lambda)*B*x for Hermitian A and positive definite B.

    Overwrites matrix A by eigenvectors, i.e., A(new, old) --
    eigenvectors are stored as ROWS in A. If if_reorder=FALSE,
    eigenenergies are given in the ascending order. */
  KMatrix& Herm_GenDiag(KMatrix& B, KMatrix& ene, bool if_reorder); 
  /*!
   \brief General matrix diagonalization Solves Left^H A Right = Lambda.
 
   Overwrites matrix. */
  void Gen_Diag(KMatrix& real_ene, KMatrix& img_ene, KMatrix& LeftVecs,
		KMatrix& RightVecs, bool if_reorder, bool if_renorm,
		double img_zero_tol=0.);
  //! Solves AX=C, puts X into C, for hermitian undefined matrix. 
  bool SolveHermUndef(KMatrix& CCCP) const;
  
  //! Sets data to 0 in a row i. 
  void SetRow(int i) const; 
  //! Sets data to 0 in a coloumn i. 
  void SetColoumn(int i) const; 

  void SwapRows(int i, int j) const;
  void SwapColoumns(int i, int j) const;
  void DeleteRow(int i);
  //void DeleteColoumn(int i);
  //! Givens rotations of rows/coloumns : x=rx*sin + ry*cos; y=rx*cos-ry*sin.
  void GivensRotRows(int rx, int ry, double costheta, double sintheta) const;
  void GivensRotColoumns(int rx, int ry, double costh, double sinth) const;
  //! Copies ith row from other to jth row of this
  void CopyRowFrom(const KMatrix& other, int i, int j) const;
  //! Copies ith coloumn from other to jth coloumn of this
  void CopyColoumnFrom(const KMatrix& other, int i, int j) const;
  
  int MaxValInRow(int i) const; 
  int MaxValInCol(int i) const;
  //! Finds the maximum value in i-th row for j>j0.
  int MaxValInRow(int i, int j0) const; 
  //! Finds the maximum value in i-th coloum for j>j0.
  int MaxValInCol(int i, int j0) const;
  //! Multiplies i-th row by coeff.
  void MultRow(int i, double coeff) const;     
  //! Multiplies i-th col by coeff.
  void MultColoumn(int i, double coeff) const; 
  //! Reorders ROWS of the matrix, to place largest values on diagonal.
  void ReOrderToDiag() const;
  void ReOrderBothToDiag(KMatrix& other);
  //! Simlest way to adjust phase: set max val==positive.
  void AdjustPhaseInRow() const;
  void AdjustPhasesInBothRows(KMatrix& other);
  void AdjustPhaseInCol() const;
  //! Scalar products of rows/columns row[i]*row[j].
  double GetRowScalarProduct(int i, int j) const;
  double GetColoumnScalarProduct(int i, int j) const;
  //! Calculates scalar product between i-th row/coloumn of *this and j-th row/coloumn of other.
  double GetScalarProduct(const KMatrix& other, int i, bool if_i_row, int j,
			  bool if_j_row) const;
  //! Mixes (adds to *this with coeff) i-th row/coloumn of *this and j-th row/coloumn of other.
  void MixVectors(const KMatrix& other, int i, bool if_i_row, int j,
		  bool if_j_row, double coeff) const;
  //! Copy to i-th row/coloumn of this scaled j-th row/coloumn of other.
  void CopyScaledVectors(const KMatrix& other, int i, bool if_i_row, int j,
			 bool if_j_row, double coeff) const;
  //! row[i]+=coeff*row[j].
  void MixRows(int i, int j, double coeff) const;
  void MixColoumns(int i, int j, double coeff) const;
  //! Schmidt orthogonalization by rows.
  void SchmidtOrthogonRows(bool if_go_backw=FALSE) const;
  //! Schmidt orthogonalization by columns.
  void SchmidtOrthogonColoumns(bool if_go_backw=FALSE) const;
  //! Check for linear dependence of rows.
  bool IfLinearDependentRows(int i, int j, double thresh=0.000001) const;
  //! Check for linear dependence of rows for the whole matrix.
  bool IfLinearDependentRows(char *ttl, FILE *out=stdout,
			     double thresh=0.000001) const;
  
  //! Read block with offsets and write with offsets
  //  void ReadFromMatrix(const KMatrix& read_from, const Index& off_put,
  //  		      const Index& offs_take, const Index& nelem,
  //  		      bool if_add) const;

  //! Determinant (LAPACK version)
  double Determinant() const;
  //! Determinant -- non-LAPACK; interface for Det() [recursive function to calculate the determinant of a square function]
  double Determinant_nonLAPACK();
  
  void Print(char *name, FILE *fl=stdout, int ndigits=6) const;
  void PrintScientific(char *name) const;
  void Read(FILE *fl);
  void Write(FILE *fl) const;
  
friend double DotProduct(const KMatrix& one, const KMatrix& two,short sign=1);
friend double RMS(const KMatrix& one, const KMatrix& two,short sign=1);
friend void Contract2(KMatrix& res, const KMatrix& one, const KMatrix& two,
		      double coeff=1., bool if_add=FALSE);
friend void Multiply2(KMatrix& res, const KMatrix& one, const KMatrix& two,
		      double coeff=1., bool if_add=FALSE);
friend void VectorMultiply2(KMatrix& res, const KMatrix& one, const KMatrix& two,
			    double coeff=1.,bool if_add=FALSE);  
friend void AddTwo2(KMatrix& res, KMatrix& one, KMatrix& two,
		    double coeff1=1.,double coeff2=1.,bool if_add=FALSE);
friend void Add2(KMatrix& res, KMatrix& one, double coeff);

  //These functions do not work, because we have wrong copy-constructor
  //Have to change default arguments for copy-constructor   
friend KMatrix operator+(const KMatrix&, const KMatrix&);
friend KMatrix operator-(const KMatrix&, const KMatrix&);
friend KMatrix operator*(const KMatrix&, const KMatrix&);
friend KMatrix operator*(double, const KMatrix&);  
friend KMatrix operator*(const KMatrix&, double);

  //! Bi-orthogonalization.
friend void BiOrthogonalizeRows(KMatrix& LeftVecs, KMatrix& RightVecs,
				KMatrix& ReEne, KMatrix& ImEne);

  //! Bi-orthoganalization to nvecs given vectors.
friend void BiOrthogonalizeToRows(int nvecs, KMatrix& LeftVecs, 
     KMatrix& RightVecs);
};
#endif




