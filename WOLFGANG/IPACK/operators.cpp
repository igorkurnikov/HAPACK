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
% $Id: operators.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Overlap Integrals} 

\index{Overlap Integrals} The following
section implements the top-level function to evaluate the
single-particle overlap matrix elements between two sets of basis
sets. We loop over all unique pairs of angular momenta.  IntegMatrix
elements for one shell of angular momenta are evaluated using the
class \name{Overlap}, which is based on \name{Recursion}. In case that
the two basis sets are identical, the function stored the newly
computed normalizing factors. 

We briefly summarize the strategy: The \name{prepare} function of the
\name{BasisSet} rearranges the basis function of the basis sets
according to their total angular momentum. The set of all functions
with the same angular momentum is called a {\em shell}. Given two sets
of such functions we can use angular momentum recursion relations to
evaluate single particle operators between them. This is done in the
recursion relation classes, which are called from the functions in
this file.

\index{integrals, overlap}
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "basis.h"
#include "function.h"
#include "integ_array.h"
#include "coef_array.h"
#include "coef_set.h"
#include "rec.h"
#include "overlap.h"
#include "kinetic.h"
#include "moment.h"
#include "nuclear.h"
#include "e_field.h"
#include "operators.h" 
#include "one_el.h"
#include "symop.h"

void overlap(InternalBasis& b1,ARRAY<IntegMat>& matrix)
{
  matrix.reset(1);      // mu index is not used ( mu = 0 )

  int nof1 = b1.no_function();
  if (matrix[0].size1() != nof1 || matrix[0].size2() != nof1)
  {
    matrix[0].reset(nof1,nof1);
    matrix[0].set(0);
  }
  compute_overlap(b1,b1,matrix,true);
}

void compute_overlap(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
		     int same)
{
  OverlapOp ovlp(b1,b2,same,matrix);

  if (same)
  {
    int nof1 = b1.no_function();
    //cout << matrix;
    b1.norm.reset(nof1); 
    for(int i=0; i< nof1;i++)
      b1.norm[i] = 1/sqrt(matrix[0](b1.min_index() + i,b1.min_index() + i));    
    // cout << "NORM: " << endl;
    //b1.norm.print(cout,10,"%8.3f");
  }
}


void overlap(InternalBasisList& blist,ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  for(int i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  matrix.reset(1);      // mu index not used ( mu = 0 )
  
  if (matrix[0].size1() != nof || matrix[0].size2() != nof)
  {
    matrix[0].reset(nof,nof);
    matrix[0].set(0);
  }

  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      compute_overlap(*blist[b1],*blist[b2],matrix,b1 == b2);
}


/*TEX
\subsection{Kinetic Energy Integrals}
\index{integrals,kinetic energy} 
The following section implements the
computation of kinetic energy matrix elements for one or two basis
sets.
*/

void kinetic(InternalBasis& b1,ARRAY<IntegMat>& matrix)
{
  int nof1 = b1.no_function();

  matrix.reset(1);        // no mu index used ( mu = 0 )
  
  if (matrix[0].size1() != nof1 || matrix[0].size1() != nof1)
  {
    matrix[0].reset(nof1,nof1);
    matrix[0].set(0);
  }
  compute_kinetic(b1,b1,matrix,true);
}

void kinetic(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
             int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  matrix.reset(0);     // no mu index used ( mu = 0 )
  
  if (matrix[0].size1() != nof1 || matrix[0].size2() != nof2)
  {
    matrix[0].reset(nof1,nof2);
    matrix[0].set(0);
  }
  compute_kinetic(b1,b2,matrix,same);
}


void kinetic(InternalBasisList& blist,ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  for(int i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  matrix.reset(1);      // no mu index is used ( mu = 0 )
  

  if (matrix[0].size1() != nof || matrix[0].size2() != nof)
  {
    matrix[0].reset(nof,nof);
    matrix[0].set(0);
  }
  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      compute_kinetic(*blist[b1],*blist[b2],matrix,b1 == b2);
}


void compute_kinetic(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
		     int same)
{
  KineticOp kin(b1,b2,same,matrix);
}
/*TEX
\subsection{Moment Integrals}
The following section implements the
computation of ALL moment integrals around a given set of \name{centers}
for one or two basis sets up to the total moment \name{mu}. 
*/

void moment(InternalBasis& b1,MomentData& data,ARRAY<IntegMat>& matrix)
{
  Momentum m;
  int nof1 = b1.no_function();

  matrix.reset(m.min(data.max_mu + 1));
  for(int i = 0; i < m.min(data.max_mu + 1); i++)
  { 
    if (matrix[i].size1() != nof1 || matrix[i].size1() != nof1)
    {
      matrix[i].reset(nof1,nof1);
      matrix[i].set(0);
    }
  }
  
  compute_moment(b1,b1,data,matrix,true);
}


void moment(InternalBasisList& blist,MomentData& data,ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  int i;
  for(i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();

  Momentum m;
  matrix.reset(m.min(data.max_mu + 1));
  for(i = 0; i < m.min(data.max_mu + 1); i++)
  {
    if (matrix[i].size1() != nof || matrix[i].size2() != nof)
    {
      matrix[i].reset(nof,nof);
      matrix[i].set(0);
    }
  }
  
  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      compute_moment(*blist[b1],*blist[b2],data,matrix,b1 == b2);
}


void moment(InternalBasis& b1,InternalBasis& b2,
	    MomentData& data,ARRAY<IntegMat>& matrix,int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  Momentum m;
  matrix.reset(m.min(data.max_mu + 1));
  for(int i = 0; i < m.min(data.max_mu + 1); i++ )
  {
    if (matrix[i].size1() != nof1 || matrix[i].size2() != nof2)
    {
      matrix[i].reset(nof1,nof2);
      matrix[i].set(0);
    }
  }
  compute_moment(b1,b2,data,matrix,same);
}


void compute_moment(InternalBasis& b1,InternalBasis& b2,
		     MomentData& data,ARRAY<IntegMat>& matrix,int same)
{
//  cout << " +++ Moment at Location: " << data.loc << " max total mu " 
//       << data.max_mu <<endl;
  MomentOp (b1,b2,same,matrix,data);      
}

/*TEX
\subsection{Nuclear Repulsion Integrals}
The following section implements the
computation of nuclear repulsion matrix elements for one or two basis
sets and a given set of nuclear \name{centers} with charges as
specified in \name{charges}. The sum of all terms is added to the matrix.
*/

void nuclear(InternalBasis& b1,Molecule& molec,ARRAY<IntegMat>& matrix)
{
  int nof1 = b1.no_function();

  matrix.reset(1);       // no mu index used ( mu = 0 )
  
  if (matrix[0].size1() != nof1 || matrix[0].size1() != nof1)
  {
    matrix[0].reset(nof1,nof1);
    matrix[0].set(0);
  }
  compute_nuclear(b1,b1,molec,matrix,true);
}


void nuclear(InternalBasisList& blist,Molecule& molec,ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  for(int i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  matrix.reset(1);
  
  if (matrix[0].size1() != nof || matrix[0].size2() != nof)
  {
    matrix[0].reset(nof,nof);
    matrix[0].set(0);
  }

  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      compute_nuclear(*blist[b1],*blist[b2],molec,matrix,b1 == b2);
}


void nuclear(InternalBasis& b1,InternalBasis& b2,
	     Molecule& molec,ARRAY<IntegMat>& matrix,int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  matrix.reset(1);
  
  if (matrix[0].size1() != nof1 || matrix[0].size2() != nof2)
  {
    matrix[0].reset(nof1,nof2);
    matrix[0].set(0);
  }
  compute_nuclear(b1,b2,molec,matrix,same);
}


void compute_nuclear(InternalBasis& b1,InternalBasis& b2,
		     Molecule& molec,ARRAY<IntegMat>& matrix,int same)
{
  int chg;
  for(int loc_indx = 0; loc_indx < molec.size(); loc_indx++)
    if (fabs((double)(chg = molec.element(loc_indx).charge())) > 0)
    {
      Location c(molec.center(loc_indx));
      for(SymOp s1(c.sym_equiv(X),c.sym_equiv(Y),c.sym_equiv(Z));
	  s1.in_range(); ++s1)
      {
	Location cc(c);
	cc.apply(s1);
	
//	cout << " +++ Nuclear Interaction for Center: " << cc << " Charge: " 
//	  << chg << endl;
	NuclearOp (b1,b2,same,matrix,cc,chg);      
      }
    }
  
}

/*TEX
\subsection{Electric Field and Electric Field Gradient Integrals}
*/

void e_field(InternalBasis& b1,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof1 = b1.no_function();

  Momentum m;
  matrix.reset(m.min(3));
  for(int i = 0; i < m.min(3); i++)
  { 
    if (matrix[i].size1() != nof1 || matrix[i].size1() != nof1)
    {
      matrix[i].reset(nof1,nof1);
      matrix[i].set(0);
    }
  }
  
  compute_e_field(b1,b1,loc,matrix,true);
}


void e_field(InternalBasisList& blist,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  int i;
  for(i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  Momentum m;
  matrix.reset(m.min(3));
  for(i = 0; i < m.min(3); i++)
  {
    if (matrix[i].size1() != nof || matrix[i].size2() != nof)
    {
      matrix[i].reset(nof,nof);
      matrix[i].set(0);
    }
  }

  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      {
//      cout << " +++ InternalBasis Combination : " << b1 << " " << b2 << endl;
	compute_e_field(*blist[b1],*blist[b2],loc,matrix,b1==b2);
      }
}

void e_field(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  Momentum m;
  matrix.reset(m.min(3));
  for(int i = 0; i < m.min(3); i++ )
  {
    if (matrix[i].size1() != nof1 || matrix[i].size2() != nof2)
    {
      matrix[i].reset(nof1,nof2);
      matrix[i].set(0);
    }
  }
  
  compute_e_field(b1,b2,loc,matrix,same);
}


void compute_e_field(InternalBasis& b1,InternalBasis& b2,
                     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
//  cout << " +++ Electric Field for Center: " << loc << endl;
  E_FieldOp (b1,b2,same,matrix,loc);      
}

/*TEX
\subsection{Angular Momentum Integrals}
*/

void angular(InternalBasis& b1,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof1 = b1.no_function();

  Momentum m;
  matrix.reset(m.min(2));
  for(int i = 0; i < m.min(2); i++)
  { 
    if (matrix[i].size1() != nof1 || matrix[i].size1() != nof1)
    {
      matrix[i].reset(nof1,nof1);
      matrix[i].set(0);
    }
  }
  
  compute_angular(b1,b1,loc,matrix,true);
}


void angular(InternalBasisList& blist,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  int i;
  for(i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  Momentum m;
  matrix.reset(m.min(2));
  for(i = 0; i < m.min(2); i++)
  {
    if (matrix[i].size1() != nof || matrix[i].size2() != nof)
    {
      matrix[i].reset(nof,nof);
      matrix[i].set(0);
    }
  }

  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      compute_angular(*blist[b1],*blist[b2],loc,matrix,b1==b2);
}


void angular(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  Momentum m;
  matrix.reset(m.min(2));
  for(int i = 0; i < m.min(2); i++ )
  {
    if (matrix[i].size1() != nof1 || matrix[i].size2() != nof2)
    {
      matrix[i].reset(nof1,nof2);
      matrix[i].set(0);
    }
  }
  
  compute_angular(b1,b2,loc,matrix,same);
}


void compute_angular(InternalBasis& b1,InternalBasis& b2,
                     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
//  cout << " +++ Angular Momentum for Center: " << loc << endl;
  AngularOp (b1,b2,same,matrix,loc);      
}

/*TEX
\subsection{Spin-Orbital Interaction Integrals}
*/

void spinorb(InternalBasis& b1,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof1 = b1.no_function();

  Momentum m;
  matrix.reset(m.min(2));
  for(int i = 0; i < m.min(2); i++)
  { 
    if (matrix[i].size1() != nof1 || matrix[i].size1() != nof1)
    {
      matrix[i].reset(nof1,nof1);
      matrix[i].set(0);
    }
  }
  
  compute_spinorb(b1,b1,loc,matrix,true);
}


void spinorb(InternalBasisList& blist,Location& loc, ARRAY<IntegMat>& matrix)
{
  int nof = 0;
  int i;
  for(i=0; i< blist.size(); i++)
    nof += blist[i] -> no_function();
  
  Momentum m;
  matrix.reset(m.min(2));
  for(i = 0; i < m.min(2); i++)
  {
    if (matrix[i].size1() != nof || matrix[i].size2() != nof)
    {
      matrix[i].reset(nof,nof);
      matrix[i].set(0);
    }
  }

  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 <= b1; b2++)
      {
//	cout << " +++ InternalBasis Combination : " << b1 << " " << b2 << endl;
	compute_spinorb(*blist[b1],*blist[b2],loc,matrix,b1==b2);
      }
}


void spinorb(InternalBasis& b1,InternalBasis& b2,
	     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();

  Momentum m;
  matrix.reset(m.min(2));
  for(int i = 0; i < m.min(2); i++ )
  {
    if (matrix[i].size1() != nof1 || matrix[i].size2() != nof2)
    {
      matrix[i].reset(nof1,nof2);
      matrix[i].set(0);
    }
  }
  
  compute_spinorb(b1,b2,loc,matrix,same);
}


void compute_spinorb(InternalBasis& b1,InternalBasis& b2,
                     Location& loc,ARRAY<IntegMat>& matrix,int same)
{
//  cout << " +++  Spin-Orbital Interaction at Center: " << loc << endl;
  SpinOrbOp (b1,b2,same,matrix,loc);      
}

/*TEX
\subsection{Normalize\label{section-normalize}}
The following function provides the normalization
\index{normalization} of the matrices for single-particle operators.
It is assumed that the norm-field of both bais sets has been
initalized by a call to overlap.
*/
void normalize(InternalBasis& b1,InternalBasis& b2,ARRAY<IntegMat>& matrix,
               int max_mu)
{
  int nof1 = b1.no_function();
  int nof2 = b2.no_function();
  int i,j;
  
  int min1 = b1.min_index();
  int min2 = b2.min_index();
  
  //cout << "NORM B1: " << min1 << " " << b1.norm << endl;
  //cout << "NORM B2: " << min2 << " " << b2.norm << endl;
  
  Momentum m;

  for(int mu = 0; mu < m.min(max_mu+1);mu++)
  for(i=0; i< nof1;i++)
  for(j=0; j< nof2;j++)
    matrix[mu](i + min1,j + min2) *= b1.norm[i]*b2.norm[j];  
}

void normalize(InternalBasisList& blist,ARRAY<IntegMat>& matrix,int max_mu)
{
  for(int b1=0; b1 < blist.size(); b1++)
    for(int b2=0; b2 < blist.size(); b2++)
      normalize(*blist[b1],*blist[b2],matrix,max_mu);
}







