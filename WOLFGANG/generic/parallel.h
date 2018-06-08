/*TEX
%
% IPACK - GENERIC
% Quantum Chemistry Project: GENERIC TYPES
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the
% license, see the LICENSE section in the file "main.c" distributed
% with this package or write to:
%      Wolfgang Wenzel, Theoretical Physics I,
%      University of Dortmund,
%      D-44221 Dortmund, Germany
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: parallel.h,v 1.2 2008/05/03 22:09:55 igor Exp $

\section{ParMachine}
*/

#ifndef PARALLEL_H
#define PARALLEL_H
#include <mpi.h>
#include <stringc.h>
#include <tarray.h>
#include <errorip.h>
#include <numer.h>
#include "vtype.h"

const int default_tag = 1;

class ParMachine
{
public:
     ParMachine(int *argc,char*** argv);
    ~ParMachine();
void abort()
     {
       error("ParMachine::abort()");       
     }
virtual int  master() const { return 1; }
virtual int  id()     const { return 0; }  
virtual int  nodes()  const { return 1; }

virtual void barrier  (const String& string)
             { cout << " +++ Synchronization Msg: " << string << endl; }
virtual int  broadcast(char* buffer,int& size,int sender) { return 0; }
virtual int  broadcast(void* buffer,int  size,int sender)  { return 0; }
virtual int  broadcast(int&  n,int sender) { return 0; }
virtual int  broadcast(Mat&  m,int sender) { return 0; }
virtual int  broadcast(NUMARRAY2D<float>&  m,int sender) { return 0; }
virtual int  nb_probe (int& flag) 
    { error("ParMachine::cant probe a simulated ParMachine");
      return 0;
    } 
virtual int  nb_send  (void* buf,int count,int node,int tag = default_tag) 
  { return 0; }
virtual int  receive  (void* buf,int& count) { return 0; }
virtual int  all_sum  (int src   ,int&    tar) { tar = src; return 0; }
virtual int  all_sum  (double src,double& tar) { tar = src; return 0; }
virtual int  all_sum  (IVec& src,IVec& tar);
virtual int  all_sum  (Mat&  src,Mat&  tar);
virtual int  all_sum  (NUMARRAY2D<float>&  src,NUMARRAY2D<float>&  tar);
virtual int  source()  const  { return 0; }
virtual int  all_to_all(IVec& in,IVec& out,int items);
virtual int  all_to_all_displ(IVec&  in,IVec& sendcount,IVec& senddispl,
			      IVec& out,IVec& recvcount,IVec& recvdispl);
virtual int  all_to_all_displ(void*  in,IVec& sendcount,IVec& senddispl,
			      void* out,IVec& recvcount,IVec& recvdispl)
    {
      error("not implemented");
      return 0;
    }
virtual int  tag()     const  
    { 
      error("ParMachine::get_count not implemented."); 
      return 0; 
    }
virtual int    error_code() const  { return 0;  }

virtual int    max_loc(double val){return id();}
virtual int    min_loc(double val){return id();}
virtual int    max_loc(int val){return id();}
virtual int    min_loc(int val){return id();}
virtual double max(double& val){return val;}
virtual double min(double& val){return val;}
virtual int    max(int& val){return val;}
virtual int    min(int& val){return val;}    
};

const int NB_CANT_SEND =  0;
const int NB_SENT      =  1;

class MPIParMachine : public ParMachine 
{
  int source_id;
  int msg_tag;
  int errcode;
  int count;
  
  int my_id;
  int numprocs;
  int world_id;
  int world_numprocs;

  ARRAY<MPI_Request*> rq;
  
public:
     MPIParMachine(int *argc,char*** argv);
    ~MPIParMachine();
  
void abort(int code = 0);
int  master() const { return my_id == 0; }
int  id()     const { return my_id; }  
int  nodes()  const { return numprocs; }
void barrier  (const String& string);
int  broadcast(char* buffer,int& size,int sender);
int  broadcast(void* buffer,int size,int sender);
int  broadcast(int&  n,int sender);
int  broadcast(Mat&  m,int sender);
int  broadcast(NUMARRAY2D<float>&  m,int sender);
int  all_to_all(IVec& in,IVec& out,int items);
int  all_to_all_displ(void*  in,IVec& sendcount,IVec& senddispl,
		      void* out,IVec& recvcount,IVec& recvdispl);
int  all_to_all_displ(IVec&  in,IVec& sendcount,IVec& senddispl,
		      IVec& out,IVec& recvcount,IVec& recvdispl);
int  nb_probe (int& flag);  
int  nb_send  (void* buf,int count,int node,int tag = default_tag);
int  receive  (void* buf,int& count);
int  all_sum  (int src,int& tar);
int  all_sum  (double src,double& tar);
int  all_sum  (IVec& src,IVec& tar);
int  all_sum  (Mat&  src,Mat&  tar);
int  all_sum  (NUMARRAY2D<float>& src,NUMARRAY2D<float>& tar);
int  source()  const     { return source_id; }
int  tag()     const     { return msg_tag;    }
int  error_code() const  { return errcode;  }

int    max_loc(double val);
int    min_loc(double val);
int    max_loc(int val);
int    min_loc(int val);
double max(double& val);
double min(double& val);
int    max(int& val);
int    min(int& val);
};

extern ParMachine* pm;

class TaskCount 
{
  int n;
public:
    TaskCount()    { n = 0; } 
int operator++(int)   { return n++; }
int local() { return (n % pm -> nodes()) == pm -> id(); }

};
  
#endif






