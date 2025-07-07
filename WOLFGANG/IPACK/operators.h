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
% $Id: operators.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Evaluation of Single Particle Operators}

Given two \name{BasisSets} $b_1$ and $b_2$, 
we can evaluate the matrix elements for a number of operators. 
Angular momentum recursion relations are used to evaluate the
operators, for a given sets of primitive basis functions with fixed angular 
total momenta on either side of the operator. Such a set is called an
angular momentum shell. The shells of primitive basis functions of
fixed angular momentum are constructed from each basis set by the
function \name{prepare}. 
\index{operators,single-particle}
\index{integrals, overlap}
\index{integrals, kinetic energy}
\index{integrals, nuclear repulsion}
\index{integrals, normalization of}
*/
#ifndef OPERATORS_H
#define OPERATORS_H
void compute_overlap(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
		     int same); 
void overlap(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
                     int same); 
void overlap(InternalBasis& b1,ARRAY<IntegMat>& matrix); 
void overlap(InternalBasisList& blist,ARRAY<IntegMat>& matrix); 

void kinetic(InternalBasisList& blist,ARRAY<IntegMat>& matrix); 
void kinetic(InternalBasis& b1,ARRAY<IntegMat>& matrix);     
void kinetic(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
                     int same);     
void compute_kinetic(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
		     int same);     

void moment(InternalBasisList& bl,MomentData& data,ARRAY<IntegMat>& matrix);
void moment(InternalBasis& b1,MomentData& data,ARRAY<IntegMat>& matrix);
void moment(InternalBasis& b1,InternalBasis& b2,
	    MomentData& data,ARRAY<IntegMat>& matrix,int same);
void compute_moment(InternalBasis& b1,InternalBasis& b2,
	    MomentData& data, ARRAY<IntegMat>& matrix,int same);


void nuclear(InternalBasisList& bl,Molecule& molec,ARRAY<IntegMat>& matrix);
void nuclear(InternalBasis& b1,Molecule& molec,ARRAY<IntegMat>& matrix);
void nuclear(InternalBasis& b1,InternalBasis& b2,Molecule& molec,
	     ARRAY<IntegMat>& matrix,int same);
void compute_nuclear(InternalBasis& b1,InternalBasis& b2,Molecule& molec,
                     ARRAY<IntegMat>& matrix,int same);

void e_field(InternalBasisList& bl,Location& loc, ARRAY<IntegMat>& matrix);
void e_field(InternalBasis& b1,Location& loc,ARRAY<IntegMat>& matrix);
void e_field(InternalBasis& b1,InternalBasis& b2,
             Location& loc, ARRAY<IntegMat>& matrix,int same);
void compute_e_field(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same);

void angular(InternalBasisList& bl,Location& loc, ARRAY<IntegMat>& matrix);
void angular(InternalBasis& b1,Location& loc,ARRAY<IntegMat>& matrix);
void angular(InternalBasis& b1,InternalBasis& b2,
             Location& loc, ARRAY<IntegMat>& matrix,int same);
void compute_angular(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same);

void spinorb(InternalBasisList& bl,Location& loc, ARRAY<IntegMat>& matrix);
void spinorb(InternalBasis& b1,Location& loc,ARRAY<IntegMat>& matrix);
void spinorb(InternalBasis& b1,InternalBasis& b2,
             Location& loc, ARRAY<IntegMat>& matrix,int same);
void compute_spinorb(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same);

void normalize(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
               int max_mu = 0);
void normalize(InternalBasisList& blist,ARRAY<IntegMat>& matrix,int max_mu = 0);

class IntegFile;

void twoel(InternalBasisList& blist,IntegFile& out);

#endif







