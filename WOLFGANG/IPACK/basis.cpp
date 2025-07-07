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
% $Id: basis.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Preparing Primitive Information for BasisSets}

This function prepares the various array which are required for the
internal of the basis set as it is used in the evaluation of the
recursion relations.  The individual array have already been discussed
in the previous sections, here the values are merely assembled. This
is accomplished by simply running though the \name{PureBasisSets} 
sets as they have been added to the \name{BasisSet} and to accumulate the 
information on all primitives gaussian orbitals.
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "const.h"
#include "basis.h"
#include "storage.h"
#include "typesip.h"
#include "function.h"
#include "symdesig.h"
#include "sphere.h"
#include "parallel.h"

int PureBasisSetList::print_flag = 0;

InternalBasis::InternalBasis(BasisSet& basis,const Location& target,
			     BasisSelectionMode mode,
			     OrbInfo& oinf) : oinfo(oinf)
     
{
  debug = basis.debug;
  global_start_index = oinfo.no_orbs();
  
  /*
     determine the maximum angular momentum
     for each momentum, gather all the radial basis functions
  */


  /*
     count the number of functions and primitive functions 
     of a given momentum
  */

  IVec nofun_tmp(::max_momentum);
  IVec noprim_tmp(::max_momentum);
  
  nofun_tmp.set(0);  
  noprim_tmp.set(0);
  
  int set,func,prim,l,l1,i;
  
  for(set = 0; set < basis.loc.size(); set++) // run through all sets
  {
    //
    // check whether the location satisfies the selection criteria
    //
    switch(mode)
    {
    case ALL:
      break;
    case ALL_CENTER:
      if (basis.loc[set] != target)
	continue;
      break;
    default:
      error(" InternalBasis::Select Mode  not implemented. ");
    }
    
    PureBasisSet& pset = *(basis.pbasis(set));
    for(func=0; func < pset.size(); func++)
    {
      PureBasisFunction& pfun = *pset(func);
      nofun_tmp [pfun.momentum()]++;
      noprim_tmp[pfun.momentum()]+= pfun.contraction();
    }
  }

  maxl = 0;
  for(i=0; i < ::max_momentum;i++)
    if (nofun_tmp[i] > 0)
      maxl = i;

  maxl++;  
  if (debug)
    cout << "Maximum Momentum: " << maxl-1 << endl;
     
  nofun .reset(maxl);
  noprim.reset(maxl);
  
  for(i= 0; i < maxl; i++)
  {
    nofun[i]   = nofun_tmp[i];
    noprim[i]  = noprim_tmp[i];
  }
  
  if (debug)
  {
    cout << "Number of Functions and Primitives: "  << endl
         << nofun << noprim << endl;
  }
  
  radials.reset(maxl);
  fcenter.reset(maxl);
  s_index.reset(maxl);
  name   .reset(maxl);

  // contraction.reset(maxl);
 
  for(l=0; l<maxl;l++)
  {
    s_index[l].reset(nofun[l]);    
    fcenter[l].reset(nofun[l]);
    name   [l].reset(nofun[l]);
  }
  
  nofun.set(0);
  int function_count = 0;
  IVec pbasiscount(::max_momentum);
  char buf[200];
  
  Momentum mcount(0);
  for(set = 0; set < basis.loc.size(); set++) // run through all sets
  {
    //
    // check whether the location satisfies the selection criteria
    //
    switch(mode)
    {
    case ALL:
      break;
    case ALL_CENTER:
      if (basis.loc[set] != target)
	continue;
      break;
    default:
      error(" InternalBasis::Select Mode  not implemented. ");
    }

    String newname = basis.name[set] + "*";
        
    PureBasisSet& pset = *(basis.pbasis(set));
    pbasiscount.set(0);
    
    for(func=0; func < pset.size(); func++)
    {
      PureBasisFunction& pfun = *pset(func);
      int ll = pfun.momentum();
      
      sprintf(buf,"%i",++pbasiscount[ll]);
      String tmpname(newname + String(buf));

      switch(ll)
      {
      case 0: tmpname = tmpname + String("-S"); break;
      case 1: tmpname = tmpname + String("-P"); break;
      case 2: tmpname = tmpname + String("-D"); break;
      case 3: tmpname = tmpname + String("-F"); break;
      case 4: tmpname = tmpname + String("-G"); break;
      case 5: tmpname = tmpname + String("-H"); break;
      case 6: tmpname = tmpname + String("-I"); break;
      default:
	sprintf(buf,"M%i",ll);
	tmpname = tmpname +  String(buf);
      }	
		
      s_index[ll][nofun[ll]] = function_count;
      fcenter[ll][nofun[ll]] = basis.loc[set];
      name   [ll][nofun[ll]] = tmpname; 
      if(flaglist(Spherical))
	function_count += 2*ll+1;
      else
	function_count += mcount.number(ll);
      for(prim=0;prim < pfun.contraction(); prim++)
        radials[ll].add(RadialFunction(basis.loc[set],pfun.exponent(prim),
				       nofun[ll],pfun.coefficient(prim)));
      nofun[ll]++;
    }    
  }
  
  for(l1=0; l1 < maxl; l1++)
    noprim[l1] = radials[l1].size();
 
#if !defined(__DECCXX)
  if (debug)
    for(l1=0; l1 < maxl; l1++)
    {
      cout << "RadialFunctions of Momentum: " << l1 << endl;
      cout << radials[l1]; 
      cout << s_index[l1]; 
    }  
#endif

  
  int count = 0;
  
  if (pm -> master() && debug)
  {
    cout << " +++ Basis Function Indices: " << endl;
    cout << "  +-----+-------------------+---------------------------+------------+-----+\n";
  
    cout << "  |   # |        Designator |                    Center |   Momentum | SYM | \n";
    cout << "  +-----+-------------------+---------------------------+------------+-----+\n";
  }
  
  for(SymDesignator s1(sym_equiv(X),sym_equiv(Y),sym_equiv(Z));
      s1.in_range(); ++s1)
  {
    for(int ll=0; ll < maxl; ll++)
    for(int fno=0; fno < nofun[ll]; ++fno)
    {
      if (!flaglist(Spherical))
      {
	for(Momentum m(ll); m(); ++m)
	{
	  if (pm -> master() && debug)
	    (ostream&) cout << "  | " 
               << setw(3) << global_start_index + count << " |" 
	       << setw(18) << name[ll][fno] << " |"
	       << fcenter[ll][fno] << " | " << m << " | ";
	  count++;	

	  SymDesignator seff(s1);

#ifdef GNU
	  for(Direction d = X; d <= Z; op_plus(d))
#else
	  for(Direction d = X; d <= Z; ++d)
#endif
	  if (sym_equiv(d))
	  {
	    if (m(d) % 2 == 0)
	    {
	      seff[d] = even;
	      if (pm -> master() && debug)
		cout << "E";
	    }
	    else 
	    {
	      seff[d] = odd;
	      if (pm -> master() && debug)
		cout << "O";
	    }
	  }
	  else
	    if (pm -> master() & debug)
	      (ostream&) cout << s1[d];

	  if (pm -> master() & debug)
	    cout << " |"<< endl;
	  
	  oinfo.record(global_start_index + count - 1,ll,seff);

	}
      }
      else
      {
	for(int m=0 ; m<2*ll+1; m++)
	{
	  if (pm -> master() && debug)
	    (ostream&) cout << "  | " 
               << setw(3) << global_start_index + count << " |" 
	       << setw(18) << name[ll][fno] << " |"
	       << fcenter[ll][fno] << " | " << setw(10) << m << " | "; 
	  count++;
	  
	  Momentum& meff = Sphere::syminfo(ll,m);
	  SymDesignator seff(s1);
	  
#ifdef GNU
	  for(Direction d = X; d <= Z; op_plus(d))
#else
	  for(Direction d = X; d <= Z; ++d)
#endif
	    if (sym_equiv(d))
	    {
	      if (meff(d) % 2 == 0)
	      {
		if (pm -> master() && debug)
		  cout << "E";

		seff[d] = even;
	      }
	      else 
	      {
		if (pm -> master() && debug)
		  cout << "O";
		seff[d] = odd;
	      }
	    }
	    else
	      if (pm -> master() && debug)
		(ostream&) cout << s1[d];

	  if (pm -> master() && debug)
	    cout << " |"<< endl;	  
	  oinfo.record(global_start_index + count - 1,ll,seff);

	}
      }
    }
  }

  if (pm -> master() && debug)
    cout << "  +-----+-------------------+---------------------------+------------+-----+\n\n";  
}

InternalBasis::~InternalBasis()
{
  if(debug)
    cout << "kill internal basis. " << endl; 
}


/*TEX
\subsection{Utilities for Basis Sets} 
  
The following functions are extremely simple, each counting the number
of functions in a basis set of returning the maximum angular momentum
encountered.
*/

int InternalBasis::no_function() const 
{
  int count = 0;
  for(SymDesignator s1(sym_equiv(X),sym_equiv(Y),sym_equiv(Z));
      s1.in_range(); ++s1)
    count++;  
  return count*no_sym_function();
}

 
int InternalBasis::no_sym_function() const 
{
  int count = 0;
  Momentum m;
  for(int i=0; i< nofun.size(); i++)
    if (flaglist(Spherical)) 
      count+= (2*i+1)*nofun(i);
    else
      count+= m.number(i)*nofun(i);
  return count;
}

int InternalBasis::sym_start_index(SymDesignator& s) const 
{
  int count = 0;
  for(SymDesignator s1(sym_equiv(X),sym_equiv(Y),sym_equiv(Z));
      s1.in_range(); ++s1)
  {
    if (s == s1)
      return count*no_sym_function();
    count++;
  }
  
  cout << " +++ Error: Illegal Symmetry Designator for Basis: " << s << endl;
  abort();
  return 0;
} 

/*
\subsection{Symmetry Equivalence}

The function \name{sym_equiv} returns \name{true} if the basis
transforms onto itself under the application of the reflection of the
coordinate in the argument.  This is relevant only for those
directions for which the symmetry flag is true. The function signals
an error if the basis set contains some functions which are invariant
and others which are not.
*/

int InternalBasis::sym_equiv(Direction d) const 
{
  if (!SymOp::symmetry(d))
    return false;

  int sym = -1;
  
  for(int l = 0; l < fcenter.size(); l++)
    for(int f = 0; f < fcenter(l).size(); f++)
      if (fcenter(l)(f).sym_equiv(d))
      {
	if (sym == -1) sym = true;
	else
	  if (sym == false)
	    error(" InternalBasis::sym_equiv :: Basis not symmetry adapted."); 
      }
      else
      {
	if (sym == -1) sym = false;
	if (sym == true)
	  error(" InternalBasis::sym_equiv :: Basis not symmetry adapted.");
      } 
  if (sym == -1)
      error(" InternalBasis::sym_equiv :: Basis is empty.");
  return sym;
}
/*TEX
\subsection{Destructors for Pointer Lists}

There are three classes in which we have stored pointers to
objects. In these containers we must delete the non-zero pointers on
exit. This implies that each pointer may be present in only one
instantiation at a time. It is presently the responsibiltiy of the
programmer to insure this.  
*/

InternalBasisList::~InternalBasisList ()
{
  for(int i=0; i< size(); i++)
    if (base[i]) delete base[i];
}

PureBasisSet::~PureBasisSet()
{ 
  for(int i=0; i< size(); i++)
    if (base[i]) delete base[i]; 
}

PureBasisSetList::~PureBasisSetList()
{ 
  for(int i=0; i< bs.size(); i++)
    if (bs[i]) delete bs[i]; 
}



/*TEX
\subsection{Printing of PureBasisFunctions}
*/

ostream& operator<<(ostream& os,const PureBasisFunction& bf)
{
  switch(bf.momentum())
  {
  case 0: os << "S   "; break;
  case 1: os << "P   "; break;
  case 2: os << "D   "; break;
  case 3: os << "F   "; break;
  case 4: os << "G   "; break;
  case 5: os << "H   "; break;
  case 6: os << "I   "; break;
  default:    
    os << "M" << setw(2) << bf.momentum();
    
  }

  char buf[100];
  os << endl 
     << "         EXPONENT         COEFFICIENT\n";

  for(int i=0; i < bf.contraction(); i++)
  {
    sprintf(buf,"  %15.6f        %10.7f \n",bf.exponent(i),
	    bf.coefficient(i));
    os << buf;
  }
  return os;
} 

void PureBasisFunction::print_gaussian_format(ostream& os,int& no_bas,
					      int& no_prim) const 
{
  os << " ";
  no_bas  += Momentum::number(momentum());
  no_prim += Momentum::number(momentum()) * contraction();
  
  switch(momentum())
  {
  case 0: os << "S   "; break;
  case 1: os << "P   "; break;
  case 2: os << "D   "; break;
  case 3: os << "F   "; break;
  case 4: os << "G   "; break;
  case 5: os << "H   "; break;
  case 6: os << "I   "; break;
  default:    
    os << "M" << setw(2) <<momentum();
  }

  char buf[100];
  os << setw(2) << contraction() << " 1.00\n";
  for(int i=0; i < contraction(); i++)
  {
    sprintf(buf,"  %16.8e  %16.8e\n",exponent(i),coefficient(i));
    os << buf;
  }
} 

void PureBasisSet::print(String& name,ostream& os)
{
  if (printed) return;
 
  int count[::max_momentum];
  int i;
  for(i=0; i< ::max_momentum; i++)
    count[i] = 0;

  os << " +++ Contracted Function Names: " << endl; 
  
  for(i=0; i < size(); i++)
  {
    int m = base[i] -> momentum();
    os << "     *** Designator: (NAME*NUMBER-MOMEMTUM): " 
       << name << "*" << ++count[m] << "-"  << *base[i] << endl;
  }
  printed = true;
  
}

void PureBasisSet::print_gaussian_format(ostream& os,int& no_bas,int& no_prim)
{
  for(int i=0; i < size(); i++)
    base[i] -> print_gaussian_format(os,no_bas,no_prim);
  os << " ****" << endl;
}
