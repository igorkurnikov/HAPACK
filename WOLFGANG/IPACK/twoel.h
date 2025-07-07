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
% $Id: twoel.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\section{Evaluation of Symmetry Adapted Integrals}

*/ 
class Coefficient_Set;
class HRR_Set;
class OputputBuffer;

#include "tquadruple.h"
#include "symop.h"
#include "symdesig.h"
#include "four_stor.h"

class HRR2;
class VRR;
class Allocator;

class TwoElectronIntegrals
{
  InternalBasis& b1;
  InternalBasis& b2;
  InternalBasis& b3;
  InternalBasis& b4;
  int same12,same34,same_pair;
  IntegFile& obuf;
  int debug;
  
  ARRAY2D<ARRAY2D<Coefficient_Set*> > cset12;
  ARRAY2D<ARRAY2D<HRR_Set*> >         hset12;
  ARRAY2D<ARRAY2D<Coefficient_Set*> > cset34;
  ARRAY2D<ARRAY2D<HRR_Set*> >         hset34;

  SStack<Quadruple<SymOp>,max_sym >   primitive;
  SStack<Quadruple<SymDesignator>,max_sym>  target;

  SStack <Quadruple<SymOp>,max_sym >  sources;
  SStack <SymOp,max_sym >             operation;

  ARRAY<FourStorage*>                 tconf;                      
  
public:
  TwoElectronIntegrals(InternalBasis& bb1,InternalBasis& bb2,
		       InternalBasis& bb3,InternalBasis& bb4,
		       int s12,int ss34,int ssp,IntegFile& oobuf);
 ~TwoElectronIntegrals();  
void logic();  
void evaluate(int l1,int l2,int l3,int l4);
void generate_targets(int l1,int l2,int l3,int l4);
void generate_targets();
void generate_primitives();
void generate_sources(Quadruple<SymOp>& q); 
void clear_coefs();  
void generate_coefs();  
void distribute(int l1,int l2,int l3,int l4,Quadruple<SymOp>& primitive,
		HRR2& hrr);
void distribute(int l1,int l2,int l3,int l4,Quadruple<SymOp>& primitive,
		VRR& vrr);
void distribute(Storage& ss,Quadruple<Momentum>& m,
		Quadruple<SymOp>& source,SymOp& operation);
void distribute(Storage& ss,Quadruple<int>& m,Quadruple<Momentum>& meff,
		Quadruple<SymOp>& source,SymOp& oper);
double compute_factor(Quadruple<SymDesignator>&  tar,Quadruple<Momentum>&  m,
		      Quadruple<SymOp>        &  src,SymOp& oper);
};


  

