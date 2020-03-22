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
% $Id: vault.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\section{Class DataVault}

This class keeps track of bits and pieces of information as a series
of calculations progresses.  The information is characterized by its
name and type.  Below we tabulate the various data which can be found
in a DataVault -- the horizontal lines seperate the data into groups
which are generated for the first ime together.
%
\begin{table}
\begin{tabular}{lcl}
\hline
   Name               &  Type     &  Meaning and possible values   \\
\hline\hline
  TOTAL\_CORE         &  int      &  total number of electrons in core   \\
  E\_NUCLEAR          &  double   &  nuclear--nuclear repulsion energy \\
  STATUS              &  String   &  every program adds a character when successful,
                      &           &  so the mrci can know whther the 'init' worked etc.\\
  MOLTYPE             &  String   &  "M" or "A"   \\
~~IPACK.SYMTABLE~~    &  IVector  &  the good old orb_info   \\
  IPACK.ORBMAP        &  IVector  &     \\
\hline
  INTERNAL\_INFO      &  IVector  &  Number of internal orbitals for each symm \\
  VIRTUAL\_INFO       &  IVector  &  Number of vitual orbitals for each symm   \\
\hline
  CORE\_INFO          &  int      &  Number of core orbitals for each symm      \
\
  NELEC\_UP           &  int      &  Number of up electrons     \\
  NELEC\_DN           &  int      &  Number of down electrons  \\
  SCF\_ENERGY         &  double   &  SCF energy   \\
  CORE\_ENERGY        &  double   &  Core energy   \\
\hline
  MB\_SYMMETRY        &  int      &  symmetry of the state in the many body calculation   \\
  EXCLUDE\_INFO       &  IVector  &  presently a name for the bsr orbitals \\
  MRCI.OPS            &  double   &  number of mega ops in the mrci   \\
  MRCI.NORM           &  double   &  norm of mrci state  \\
  MRCI.ENERGY         &  double   &  energy of mrci state (excluding core and nur-rep)   \\
  XITER.OPS           &  double   &  megaops in xiter   \\
  XITER.NORM          &  double   &  norm of state after xiter   \\
  XITER.ENERGY        &  double   &  energy of mrci state after xiter   \\
\hline
  BSR.BARE\_CORR      &  double   &  perturbative correction for bsr   \\
  BSR.NORM\_CORR      &  double   &  bsr norm contribution   \\
  BSR.ENXX\_CORR      &  double   &  <psi|h_xx|psi> in BSR   \\
  BSR.OPS             &  double   &  mega-ops for BSR   \\
\hline
\end{tabular}
\ec
\end{table}
%
The constructor creates a new DataVault or reads an existing DataVault
from the disk-file "Data".

 and write the existing
*/
#ifndef VAULT_H
#define VAULT_H

#include <numer.h>
#include <string.h>
#include <tsstack.h>
#include "vtype.h"

enum  DataVaultCode  { vault_ignore,vault_warning,vault_severe };   
enum  DataVaultState { vault_old,vault_new };
enum  DataVaultType  { v_string,v_boolean,v_int,v_double,v_vector,
		       v_matrix,v_ivector,v_imatrix,v_svector };       

ostream& operator<<(ostream& os,DataVaultType& t);
istream& operator>>(istream& os,DataVaultType& t);

class DataVault
{
  SStack<       String,100> nname;
  SStack<DataVaultType,100> ntype;
  SStack<       String,100> ncont;
  String                    filename;

int     find(const String& name,DataVaultType = v_string,
	     DataVaultCode=vault_ignore);
void    unpack(char** ptr,double& d);
void    unpack(char** ptr,long& d);
void    pack(String& string,double d);
void    pack(String& string,long d);
  
public:
        DataVault(String name = String("Data"),
		  DataVaultState newcode = vault_old);
        DataVault(istream& is) { read(is); }
       ~DataVault();
  
void    read (istream& is);
void    write(ostream& os); 
void    read ();
void    write();

int     check    (const String& name);
int     check    (const String& name,DataVaultType type); 
int     insert   (const String& name,const String& cont,
		  DataVaultType type = v_string, 
		  DataVaultCode code = vault_severe);
String  retrieve (const String& name,DataVaultCode code = vault_severe);
//
// type specific items
//
void    insert   (const String& name,double  d,
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,double& d);  

void    insert   (const String& name,long   d, 
		  DataVaultCode code = vault_ignore);  
void    insert   (const String& name,int   d, 
		  DataVaultCode code = vault_ignore)
                 { insert(name,(long) d,code); } 
void    retrieve (const String& name,long&  d);
void    retrieve (const String& name,int &  d);

void    insert   (const String& name,const Vec&  d,
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,Vec& d);  

void    insert   (const String& name,const NUMARRAY<int>&  d,
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,NUMARRAY<int>& d);  

void    insert   (const String& name,const NUMARRAY2D<int>&  d,
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,NUMARRAY2D<int>& d);  

void    insert   (const String& name,const Mat&  d,
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,Mat& d);  

void    insert   (const String& name,const ARRAY<String>& s, 
		  DataVaultCode code = vault_ignore);  
void    retrieve (const String& name,ARRAY<String>& s);
void    print    (ostream& os);  
};

#endif
