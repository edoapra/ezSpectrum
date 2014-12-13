/*! 
\file tmp_buffer.C
\ingroup (MATRIX_MATH)
AIK, 1998
*/

#include "tmp_buffer.h"

//! Default minsize for buffer (512 K for doubles).
#define MINSIZE (64*1024)      
//! Default maxsize for buffer (256 MB for doubles).
#define MAXSIZE (32*1024*1024)

///// Buffer for integers  /////////

TmpBufferInt::TmpBufferInt() : minsize(MINSIZE),
maxsize(MAXSIZE), in_use(false),buff_size(0),tmp_buffer(NULL)
{
  //tmp_buffer=new T[buff_size];
}

//For debug purposes (to track the use of uninitialized memory: after
//buffer is freed, fill it with large negative numbers
void TmpBufferInt::PutDeadBeef() const
{
  int dead_beef=666666;
  for(int i=0; i<buff_size; i++)
    tmp_buffer[i]=dead_beef;
}

void TmpBufferInt::Resize(int buffs)
{
#ifdef DEBUG
  check(!(buffs>=0),"TmpBufferInt::Resize() : invalid size");
#endif
  if(NULL!=tmp_buffer)
    delete [] tmp_buffer;
  buff_size=(buffs/minsize+((buffs%minsize)!=0))*minsize;
  check(buffs>buff_size,"TmpBufferInt::Resize() : insufficient buff_size");
  tmp_buffer=new int[buff_size];
}

void TmpBufferInt::DeleteBuff() //deallocate buffer 
{
  check(in_use,"TmpBufferInt::DeleteBuff() : Buffer is in use");
  buff_size=0;
  if(NULL!=tmp_buffer)
    {
      delete [] tmp_buffer;
      tmp_buffer=NULL;
    }
}

TmpBufferInt::~TmpBufferInt()
{
  DeleteBuff();
}

int *TmpBufferInt::GetBuff(int buffs, bool if_set) 
{
  check(in_use,"TmpBufferInt::GetBuff() : Buffer is already in use");
  in_use=true;
  if(buffs>buff_size)
    Resize(buffs);
  if(if_set)
    memset(tmp_buffer,0,sizeof(int)*buffs);
#ifdef DEBUG
  else
    PutDeadBeef();
#endif
  return tmp_buffer;
}

TmpBufferInt& TmpBufferInt::theBuffer()
{
  static TmpBufferInt TheTmpBuffer;
  return TheTmpBuffer;
}


///// Buffer for doubles  /////////

TmpBufferDouble::TmpBufferDouble() : minsize(MINSIZE),
maxsize(MAXSIZE), in_use(false),buff_size(0),tmp_buffer(NULL)
{
  //tmp_buffer=new T[buff_size];
}

//For debug purposes (to track the use of uninitialized memory: after
//buffer is freed, fill it with large negative numbers
void TmpBufferDouble::PutDeadBeef() const
{
  double dead_beef=666666.;
  for(int i=0; i<buff_size; i++)
    tmp_buffer[i]=dead_beef;
}

void TmpBufferDouble::Resize(int buffs)
{
#ifdef DEBUG
  check(!(buffs>=0),"TmpBufferDouble::Resize() : invalid size");
#endif
  if(NULL!=tmp_buffer)
    delete [] tmp_buffer;
  buff_size=(buffs/minsize+((buffs%minsize)!=0))*minsize;
  check(buffs>buff_size,"TmpBufferDouble::Resize() : buff_size is too small");
  tmp_buffer=new double[buff_size];
}

void TmpBufferDouble::DeleteBuff() //deallocate buffer 
{
  check(in_use,"TmpBufferDouble::DeleteBuff() : Buffer is in use");
  buff_size=0;
  if(NULL!=tmp_buffer)
    {
      delete [] tmp_buffer;
      tmp_buffer=NULL;
    }
}

TmpBufferDouble::~TmpBufferDouble()
{
  DeleteBuff();
}

double *TmpBufferDouble::GetBuff(int buffs, bool if_set) 
{
  check(in_use,"TmpBufferDouble::GetBuff() : Buffer is already in use");
  in_use=true;
  if(buffs>buff_size)
    Resize(buffs);
  if(if_set)
    memset(tmp_buffer,0,sizeof(double)*buffs);
#ifdef DEBUG
  else
    PutDeadBeef();
#endif
  return tmp_buffer;
}

TmpBufferDouble& TmpBufferDouble::theBuffer()
{
  static TmpBufferDouble TheTmpBuffer;
  return TheTmpBuffer;
}
