/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
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
% $Id: four_stor.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\section{Management of Four-Index Arrays of \name{Storage}}

\subsection{Description}  
This section implements a class to manage the storage of four-index
arrays of \name{Storage} array. These classes are used for symmetry
adaption and spherical symmetrization. We implement one polymorphic
class where zero to four of the indices may be of type \name{Momentum}
to label cartesian angular momentum indices or alternately of type
\name{int} to label spherical angular momentum indices. The types of
the indices are given by an flag of type \name{FOURSTORTYPE}  , which
takes the values \name{T0},...,\name{T4} , which indicate the presence
of zero to four indices of type \name{Momentum} respectively.  

The class maintains an array of \name{pointers} to objects of type
\name{Storage}. The destructor will destroy all objects with non-zero
pointers at exit. The relevant member functions are: 
\bi 
\item Constructor (default) Initialize an empty array of type
      \name{T0}. Used to allow arrays of type \name{FourStorage}.
\item \name{find(m1,m2,m3,m4)}, returns the pointer  indexed by the
      arguments. This function is polymorphic, type checking is implemented.
\item \name{insert(m1,m2,m3,m4,p)}, Inserts the pointer p in the element
      indexed by the
      arguments. This function is polymorphic, type checking is implemented.
\item \name{reset(l1,l2,l3,l4,type)} clear the present array and resets the
dimensions according to the values of the angular momenta and the give
type for the array.
\ei
\index{FourStorage} \index{Angular Momentum Arrays, Storage}  
\subsection{Header} 
*/
#ifndef FOURSTOR_H
#define FOURSTOR_H


enum FOURSTORTYPE { T0,T1,T2,T3,T4 }; 

class FourStorage
{
  friend class TwoElectronSphere;
  
  FOURSTORTYPE fst;
  Allocator         alloc;
  D_Array<Storage*> fields;
  Quadruple<int>    max;
  Quadruple<int>    no;
  Quadruple<int>    sphere;
  
int      compute_index(Momentum& m1,Momentum& m2,Momentum& m3,Momentum& m4)
         {
	   assert(fst == T0);
	   assert(m1.l() == max[0] && m2.l() == max[1] && m3.l() == max[2]
		  && m4.l() == max[3]);
	   
	   return ((m1.local_index()  * no[1] + m2.local_index())*no[2] +
		       m3.local_index()) * no[3] + m4.local_index();
	 }
int      compute_index(int m1,Momentum& m2,Momentum& m3,Momentum& m4)
         {
	   assert(fst == T1);
	   assert(m1 < max[0] && m2.l() == max[1] && m3.l() == max[2]
                  && m4.l() == max[3]);
	   return ((m1  * no[1] + m2.local_index())*no[2] +
		    m3.local_index()) * no[3] + m4.local_index();
	 }
int      compute_index(int m1,int m2,Momentum& m3,Momentum& m4)
         {
	   assert(fst == T2);
	   assert(m1 < max[0] && m2 < max[1] && m3.l() == max[2]
                  && m4.l() == max[3]);
	   return ((m1  * no[1] + m2)*no[2] +
		    m3.local_index()) * no[3] + m4.local_index();
	 }
int      compute_index(int m1,int m2,int m3,Momentum& m4)
         {
	   assert(fst == T3);
	   assert(m1 < max[0] && m2 < max[1] && m3 < max[2] && 
		  m4.l() == max[3]);
	   return ((m1  * no[1] + m2)*no[2] + m3) * no[3] + m4.local_index();
	 }
int      compute_index(int m1,int m2,int m3,int m4)
         {
	   assert(fst == T4);
	   assert(m1 < max[0] && m2 < max[1] && m3 < max[2] && m4 < max[3]);
	   return ((m1  * no[1] + m2)*no[2] + m3) * no[3] + m4;
	 }

public:
void     clear();  
         FourStorage(int s1,int s2,int s3,int s4,int nb) :alloc(s1,s2,s3,s4,nb)
	   { fst = T0; max[0] = max[1] = max[2] = max[3] = 0; } 
        ~FourStorage() { clear(); }  

void     reset(int m1,int m2,int m3,int m4,FOURSTORTYPE ft);

Storage* find  (Momentum& m1,Momentum& m2,Momentum& m3,Momentum& m4)
         { return fields[compute_index(m1,m2,m3,m4)]; }    
void     insert(Momentum& m1,Momentum& m2,Momentum& m3,Momentum& m4,
		Storage* ptr)
         { fields[compute_index(m1,m2,m3,m4)] = ptr; }    

Storage* find  (int m1,Momentum& m2,Momentum& m3,Momentum& m4)
         { return fields[compute_index(m1,m2,m3,m4)]; }    
void     insert(int m1,Momentum& m2,Momentum& m3,Momentum& m4,Storage* ptr)
         { fields[compute_index(m1,m2,m3,m4)] = ptr; }    

Storage* find  (int m1,int m2,Momentum& m3,Momentum& m4)
         { return fields[compute_index(m1,m2,m3,m4)]; }    
void     insert(int m1,int m2,Momentum& m3,Momentum& m4,Storage* ptr)
         { fields[compute_index(m1,m2,m3,m4)] = ptr; }    

Storage* find  (int m1,int m2,int m3,Momentum& m4)
         { return fields[compute_index(m1,m2,m3,m4)]; }    
void     insert(int m1,int m2,int m3,Momentum& m4,Storage* ptr)
         { fields[compute_index(m1,m2,m3,m4)] = ptr; }    

Storage* find  (int m1,int m2,int m3,int m4)
         { return fields[compute_index(m1,m2,m3,m4)]; }    
void     insert(int m1,int m2,int m3,int m4,Storage* ptr)
         { fields[compute_index(m1,m2,m3,m4)] = ptr; }    
Allocator& allocator() { return alloc; }

};


ostream& operator<<(ostream& os, const FourStorage& fs);

#endif
