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
% $Id: f.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Computation of Auxialiry Error Functions} 

This set of routines deals with the computation of:
\be
   F(T,m) = \int_0^1 dt \ t^{2m} \exp(-T t^2)
\ee
for all values of T and $0 \le m \le {\rm MMAX}.$ 
This functiom must be evaluated for all quadruplets of primitve basis
functions several times (sum of all agular momenta), hence a significant 
fraction of the computational effort is associated with the evaluation of 
this function. 

For a given $m$ all arguments for $t < T_m$ are evaluated as a sixth
order Taylor-expansion. For $t > T_m$ we use an approximate formula. 
This strategy is described by OS. The evaluation of the function for a
single argument, as well as for an array of values is accomplished by 
the function \name{F} \index{F}. 
*/

extern ARRAY<Vec*> ftable;

void   read_ftable(String& filename);  
void   init_ftable();  
void   clean_ftable();  
double fastcomp(IntegType T,int m);
double F(IntegType  T,int m);
void   F(IntegType* T,IntegType* result,int no,int m);

void   compute_ftable(String& filename);


extern double* binbin;
