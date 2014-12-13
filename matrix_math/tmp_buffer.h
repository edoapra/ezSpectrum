/*! 
\file tmp_buffer.h
\brief This file contains tools to manage tmp-buffers for low level routines.
\ingroup MATRIX_MATH

AIK, 1998.

   Only one buffer at a time can exist. Use in low-level routines only,
   as a local tmp-variables.
   Now use it in Matrix amd Tensor for tmp-arrays (blas, permutaions, etc).

   AIK, 09/2000: I am changing this class from template to two explicit
   classes (for ints and doubles), since on some platforms we have problems
   with multiple definitions for templates, which can cause serious problems
   for functions with static variables, as we have here.
*/

#ifndef _tmp_buffer_h
#define _tmp_buffer_h

#include "genincludes.h"
#include "genutil.h"

/*! 
\brief TMP-buffer of integers.
\ingroup (TENSOR)
*/
class TmpBufferInt
{
  //!chunk for resize
  int minsize;
  //!max memory can be allocated
  int maxsize; 
  bool in_use;
  //!total size of allocated memory
  int buff_size; 
  int *tmp_buffer;
  
  TmpBufferInt();
  void Resize(int buffs);
  bool IfInUse() const { return in_use; } 
  TmpBufferInt& operator=(const TmpBufferInt& );
  //!For debug purposes
  /*! To track the use of uninitialized memory: after
    buffer is freed, fill it with large negative numbers */
  void PutDeadBeef() const;
public:
  ~TmpBufferInt();
  //! returns TMP-buffer 
  int *GetBuff(int buffs, bool if_set=false); 
  void FreeBuff() { in_use=false; }
  //! deallocate buffer 
  void DeleteBuff(); 
  int Size() const { return buff_size; }
  void SetMaxSize(int maxsz) { maxsize=maxsz; }
  void SetMinSize(int minsz) { minsize=minsz; }
  int GetMaxSize() const { return maxsize; } 
  static TmpBufferInt& theBuffer(); 
};

/*! 
\brief TMP-buffer of doubles.
\ingroup (TENSOR)
*/
class TmpBufferDouble
{
  //!chunk for resize
  int minsize;
  //! max memory can be allocated
  int maxsize; 
  bool in_use;
  //! total size of allocated memory
  int buff_size; 
  double *tmp_buffer;
  
  TmpBufferDouble();
  void Resize(int buffs);
  bool IfInUse() const { return in_use; } 
  TmpBufferDouble& operator=(const TmpBufferDouble& );
  //!For debug purposes
  /*! To track the use of uninitialized memory: after
    buffer is freed, fill it with large negative numbers */
  void PutDeadBeef() const;
public:
  ~TmpBufferDouble();
  //! returns TMP-buffer 
  double *GetBuff(int buffs, bool if_set=false); 
  void FreeBuff() { in_use=false; }
  //! Deallocate buffer 
  void DeleteBuff(); 
  int Size() const { return buff_size; }
  void SetMaxSize(int maxsz) { maxsize=maxsz; }
  void SetMinSize(int minsz) { minsize=minsz; }
  int GetMaxSize() const { return maxsize; } 
  static TmpBufferDouble& theBuffer(); 
};
#endif


