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
% $Id: input.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\section{Input of Basis Sets}

\subsection{Error Reporting}

The first function deals only with parse errors. It always terminates
the program. Better safe than sorry.
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "basis.h"
#include "function.h"
#include "symop.h"


void ioerror(istream& is,char* msg)
{
  cerr << " +++ IO Error: " << msg << endl; 
  if (is.eof())
     cerr << "     Unexpected End of File. " << endl;
  else
    if (is.fail())
      cerr << "     Stream IO Error. General Failure. " << endl;
    else
      {
	cerr << "     Parse Error near: \"";
	char buf[20];
	is.get(buf,20);
	cerr << buf << "\"";
      }
  abort();
}

/*TEX
\subsection{Input of \name{PureBasisFunctions}} 
  
A \name{PureBasisFunctions} with exponents $\zeta_1,\zeta_2,\zeta_3,\dots$
and contraction coeffcicients $c_1,c_2,\dots$ is entered as:
   $$ \{\ \zeta_1\  c_1\ \zeta_2\ c_2\ \dots\ \} $$ 
All whitespace is ignored.
\index{basis function, input of pure}
*/

istream& operator>>(istream& is,PureBasisFunction& bf)
{
  char c;
  double xx,cc;
  is >> c;
  if (c != '{')
    ioerror(is,"Input PureBasisFunction: Opening Brace Expected. ");
  
  is >> c;
  if (!is) 
    ioerror(is,"Input PureBasisFunction: General Read Error. ");
  
  while(c != '}')
  {
    is.putback(c);
    if (!is) 
      ioerror(is,"Input PureBasisFunction: General Read Error. ");
    is >> xx;
    if (!is) 
      ioerror(is,"Input PureBasisFunction: Exponent Expected. ");
    is >> cc;
    if (!is) 
      ioerror(is,"Input PureBasisFunction: Coefficient Expected. ");
    bf.add(xx,cc);
    is >> c;
  }
  if (bf.contraction() == 1 && fabs(bf.coefficient(0) - 1.0) > delta)
    error("BasisFunction Error: The contraction coeffcient of a function\n\
                        with a single exponent must be 1.0");

  return is;
}
/*TEX
\subsection{Input of \name{PureBasisSets}} 
  
A \name{PureBasisSet} is entered as a series of 
\name{PureBasisFunctions} $f_i$, with angular momentum $l_i$ and exponents 
exponents $\zeta_{i1},\zeta_{i2},\zeta_{i3},\dots$
and contraction coefficients $c_{i1},c_{i2},\dots$ is entered as:

\ba
\{ &  l_1 \ \{\ \zeta_{11}\  c_{11}\ \zeta_{12}\ c_{12}\ \dots\ \}&\nonumber \\
   &  l_2 \ \{\ \zeta_{21}\  c_{21}\ \zeta_{22}\ c_{22}\ \dots\ \}&\nonumber \\
   &   \dots &\} \nonumber
\ea
All whitespace is ignored. Presently allowed are values for $0 \le l_i < 10$,
which may be given as $s,p,d,f,g,h,i,j,k,l$ or their upper-case variants, respectively. 
Futhermore it is possible to include an entire existing \name{PureBasisSet}
with the symbolic name NAME into the present set by entering  the  sequence 
$$  \# {\rm NAME} $$
anywhere in the list of \name{PureBasisFunctions}.
\index{basis sets, input of pure}
\index{basis sets, combining}

*/

int PureBasisSet::parse(istream& is,PureBasisSetList* list)
{
  char c;
  is >> c;
  if (c != '{')
    ioerror(is,"Input PureBasisSet: Opening Brace Expected. ");

  is >> c;  
  if (!is) 
    ioerror(is,"Input PureBasisSet: General Read Error. ");

  while(c != '}')
  { 
    int l;
    switch(c)
    {
    case 's':
    case 'S':
    case '0':
      l = 0;
      break;
    case 'p':
    case 'P':
    case '1':
      l = 1;
      break;
    case 'd':
    case 'D':
    case '2':
      l = 2;
      break;
    case 'f':
    case 'F':
    case '3':
      l = 3;
      break;
    case 'g':
    case 'G':
    case '4':
      l = 4;
      break;
    case 'H':
    case 'h':
    case '5':
      l = 5;
      break;
    case 'I':
    case 'i':
    case '6':
      l = 6;
      break;
    case 'j':
    case 'J':
    case '7':
      l = 7;
      break;
    case 'k':
    case 'K':
    case '8':
      l = 8;
      break;
    case 'l':
    case 'L':
    case '9':
      l = 9;
      break;
    case '#':
    {
      String name;
      is >> name;
      if (!is) 
	ioerror(is,"Input PureBasisSet: Symbolic Name expected.");
      PureBasisSet* ptr = 0;
      if (list)
	ptr = list -> find(name);
      if (!ptr)
      {
	cerr << "Symbolic Name: " << name << " not found." << endl;
	abort();
      }  
      add(*ptr);
      continue;
    }
    default:  
      is.putback(c);
      ioerror(is,"Input PureBasisSet: Expected Angular Momentum Designator.");
    }

    PureBasisFunction* bf = new PureBasisFunction(l);
    is >> *bf;
    add(bf);
    is >> c;
  } 
  return true;
}

/*TEX
\subsection{Input of PureBasisSetLists}  
  
A \name{PureBasisSetList} is entered as a series of
\name{PureBasisSet} $s_i$ in the following fashion. The input field is
preceded by a keyword, which designates what comes next.  

\begin{tabular}{ll} 
END & terminates the input of \name{PureBasisSets}. \\
BASIS {\em name} $s_i$ & enters the \name{PureBasisSet} $s_i$. \\
 & $s_i$ designates the series of tokens described in the section  \\
 & on input of such basis sets. \\
INCLUDE path & reads an entire file with name \name{path}  \\
             & of \name{PureBasisSets} which contains records just as  \\
             & those defined here. This file must terminate in an END line.
\end{tabular}
\index{PureBasisSets, input}
*/


istream& operator>>(istream& is,PureBasisSetList& list)
{
  String token;
  
  String token_comment      ("/*");
  String token_end_comment  ("*/");
  String token_end      ("END");
  String token_include  ("INCLUDE");
  String token_basis    ("BASIS");
   
  if (pm->master())
    cout << " +++ Reading Pure Basis Sets: " << endl; 

  while(1)
  {
    is >> token;
    if (!is)
      ioerror(is," PureBasisSetList::General Read Error. ");
    
    if (token == token_comment)
    {
      while(token != token_end_comment)
      {
	is >> token;
	if (!is)
	  ioerror(is," BasisSet::General Read Error parsing comment. ");
      }
      continue;
    }

    if (token == token_end) return is;
    if (token == token_include)
    {
      is >> token;
      if (!is)
	ioerror(is," PureBasisSetList::General Read Error. ");
      cout << " +++ PureBasisSetList::read including file: " << token << endl;
      char *fname = token.contents();
      ifstream infile(fname,ios::in);
      delete fname;
      
      if (!infile)
	error("IO Error -- Could not open file. ");
      infile >> list;
      continue;
    }
    if (token == token_basis)
    {
      is >> token;
      if (!is)
	ioerror(is," PureBasisSetList::General Read Error. ");
      if (list.print_flag && pm -> master())
	cout << " +++ Reading PureBasisSet named: " << token << endl;

      PureBasisSet* newbasis = new PureBasisSet;
      
      newbasis -> parse(is,&list);
      if (list.print_flag && pm -> master())
	(ostream&) cout << *newbasis;
      list.add(token,newbasis);
      continue;
    }
    cerr << " PureBasisSetList::read token: " << token << endl;
    ioerror(is," PureBasisSetList::Illegal Token Encountered. ");
  }
  
  return is;
} 

/*TEX
\subsection{Reading BasisSets and Integral Information}


A \name{BasisSet} is specified as a list of records of the format:
\bc
 CENTER \hspace{1em} \ x\ y\ z\hspace{1em} E \hspace{1em} BNAME
\ec
%
where CENTER is a keyword, $x,y,z$ are the locations of the center,
E is the name of the element  (X indicates ghost centers) and \name{BNAME} is a
symbolic name for a \name{PureBasisSet}, which must be present in the
\name{PureBasisSetList} \name{list}. The keywords KINETIC, NUCLEAR
TWOEL, MOMENT, EFIELD, ANGMOM and
SPINORB followed by the word ON or OFF enable to select the
calculated integrals. Overlap integrals must always be calculated for 
normalization.
In case of MOMENT ON the location and the maximal
total moment maxM must be defined, too:
\bc
 MOMENT\hspace{1em}ON\hspace{1em} \ x\ y\ z\hspace{1em} maxM.
\ec
Two more keywords are ANGSTROM and SPHERICAL followed by ON or OFF, which
allow to decide between angstrom or atomic units and spherical or cartesian
functions, respectively.
The list is terminated by the keyword END.  The function returns the list of
centers, which have nonzero charge and the respective charges collected
in molecule and the data for moment integrals.
\index{BasisSet,input}
\index{Flags, setting}
\index{flaglist, declaration}


*/

void BasisSet::parse(istream& is,Molecule& molecule,MomentData& mu_data,
		     Location& e_field_loc,Location& anglr_loc,
                     Location& spin_orb_loc,PureBasisSetList& list)
{
  if (pm->master())
    cout << " +++ Reading Information on Integrals, Basis Sets, Centers "
	 << "and Symmetry: " << endl; 

  String token;

  String token_comment      ("/*");
  String token_end_comment  ("*/");
  String token_end          ("END");
  String token_kinetic      ("KINETIC");
  String token_nuclear      ("NUCLEAR");
  String token_twoelec      ("TWOEL");
  String token_moment       ("MOMENT");
  String token_efield       ("EFIELD");
  String token_angular      ("ANGMOM");
  String token_spinorb      ("SPINORB");
  String token_spherical    ("SPHERICAL");
  String token_angstrom     ("ANGSTROM");
  String token_symmetry     ("SYMMETRY");
  String token_center       ("CENTER");

  int center_count = 0;
  int tot_charge   = 0;
  int no_prim      = 0;
  int no_bas       = 0;
  
  ofstream  basis_list("BASIS_LIST",ios::out);
  ofstream center_list("CENTER_LIST",ios::out);
  if (pm -> master())
  {
    center_list 
      << " This is part of the Gaussian 94(TM) system of programs.\n\n"  
      << " *********************************************\n"
      << " Gaussian 94:  x86-Linux-G94RevD.3  1-May-1996\n"
      << "                   Standard orientation:\n" 
      << " ----------------------------------------------------------\n"
      << " Center     Atomic              Coordinates (Angstroms)\n"
      << " Number     Number             X           Y           Z\n"
      << " ----------------------------------------------------------\n";
    basis_list <<  " Basis set in the form of general basis input:\n";
  }
  
  while(1)
  {
    is >> token;
    if (!is)
      ioerror(is," BasisSet::General Read Error. ");
    if (token == token_comment)
    {
      while(token != token_end_comment)
      {
	is >> token;
	if (!is)
	  ioerror(is," BasisSet::General Read Error parsing comment. ");
      }
      continue;
    }
    
    if (token == token_end) break;

    if (token == token_kinetic)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Kinetic Integrals are turned off. " << endl;	
        flaglist[Kinetic] = false;  
        continue;
      }
      if(token == "ON")
      {
	if (pm->master())
	  cout << " +++ Kinetic Integrals are turned on. " << endl;	
        flaglist[Kinetic] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_nuclear)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Nuclear Integrals are turned off. " << endl;	
        flaglist[Nuclear] = false;  
        continue;
      }
      if(token == "ON")
      {
	if (pm->master())
	  cout << " +++ Nuclear Integrals are turned on. " << endl;	
        flaglist[Nuclear] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_twoelec)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Two Electron Integrals are turned off. " << endl;
        flaglist[TwoElec] = false;  
        continue;
      }
      if(token == "ON")
      {
	if (pm->master())
	  cout << " +++ Two Electron Integrals are turned on. " << endl;
        flaglist[TwoElec] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_moment)
    {
      Location loc;
      int      max_mu;

      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ MOMENT Integrals are turned off. " << endl;	
        flaglist[Moment] = false;  
        continue;
      }
      if(token == "ON")
      {
	is >> loc;
	if (!is)
	  ioerror(is," BasisSet::Location expected. ");

	is >> max_mu;
	if (!is)
	  ioerror(is," BasisSet::MAX_MU Integer expected. ");
      
	if (pm->master())
	  cout << " +++ Moment Integrals at " << loc << " up to mu = "
	       << max_mu <<endl;
	mu_data.loc = loc;
	mu_data.max_mu = max_mu;
	flaglist[Moment] = true;  
	continue;
      }
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }

    if (token == token_efield)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Electric Field Integrals are turned off. " << endl; 
        flaglist[eField] = false;  
        continue;
      }
      if(token == "ON")
      {
	is >> e_field_loc;
	if (!is)
	  ioerror(is," BasisSet::Location expected. ");

	if (pm->master())
	  cout << " +++ Electric Field Integrals at " << e_field_loc << endl;   

        flaglist[eField] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_angular)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Angular Momentum Integrals are turned off. " << endl;	
        flaglist[Angular] = false;  
        continue;
      }
      if(token == "ON")
      {
	is >> anglr_loc;
	if (!is)
	  ioerror(is," BasisSet::Location expected. ");

	if (pm->master())
	  cout << " +++ Angular Momentum Integrals at " << anglr_loc << endl;   

        flaglist[Angular] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_spinorb)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Spin Orbital Interaction Integrals are turned off. " 
	       << endl;	
        flaglist[SpinOrb] = false;  
        continue;
      }
      if(token == "ON")
      {
	is >> spin_orb_loc;
	if (!is)
	  ioerror(is," BasisSet::Location expected. ");

	if (pm->master())
	  cout << " +++ Spin Orbital Interaction Integrals at " 
	       << spin_orb_loc   << endl;
        flaglist[SpinOrb] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    
    if (token == token_spherical)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Using cartesian functions. " << endl;	
        flaglist[Spherical] = false;  
        continue;
      }
      if(token == "ON")
      {
	if (pm->master())
	  cout << " +++ Using spherical functions. " << endl;	
        flaglist[Spherical] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }
    if (token == token_angstrom)
    {
      is >> token;
      if (!is)
        ioerror(is," BasisSet::General Read Error. ");
      if(token == "OFF")
      {
	if (pm->master())
	  cout << " +++ Distance Units in Atomic Units. " << endl;	
        flaglist[Angstrom] = false;  
        continue;
      }
      if(token == "ON")
      {
	if (pm->master())
	  cout << " +++ Distance Units in Angstrom. " << endl;	
        flaglist[Angstrom] = true;  
        continue;
      }      
      ioerror(is," BasisSet:: Read Error: ON or OFF expected. ");
      continue;
    }

    if (token == token_symmetry) 
    {
      Direction d;
      is >> d;
      if (!is)
	ioerror(is," Basis Set::Symmetry Direction expected. ");
      SymOp::set_symmetry(d);
      if (pm->master())
	cout << " +++ Switch on Reflection Symmetry in Direction: " << d 
	     << endl;
      continue;
    }
    
    if (token == token_center) 
    {
      Location loc;
      is >> loc;
      if (!is)
	ioerror(is," BasisSet::Location expected.");
      Element elem;
      is >> elem;
      if (!is)
	ioerror(is," BasisSet::Element expected.");
      is >> token;
      if (!is)
	ioerror(is," BasisSet::Basis Set Name expected.");
      String bname(elem.name());
      bname = bname + ":";
      bname = bname + token;
      if (pm->master())
	cout << " +++ Center: " << loc << " Element: " << elem 
	     << " BasisSet: " << token << endl;
      PureBasisSet* pb = list.find(bname);
      if (!pb)
      {
	cerr << "  BasisSet::read PureBasisSet name: " << bname << endl;
	ioerror(is,"  Basis Set Not Found. ");
      }
      add(pb,bname,loc);
      molecule.add(elem,loc,bname);

      char buf[100];
      for(SymOp s1(loc.sym_equiv(X),loc.sym_equiv(Y),loc.sym_equiv(Z));
	  s1.in_range(); ++s1)
      {
	tot_charge += elem.charge();
	Location tmp(loc);
	tmp.apply(s1);
	sprintf(buf,"           %8.4f    %8.4f    %8.4f",
		bohr*tmp.x(),bohr*tmp.y(),bohr*tmp.z());
	center_list << setw(5) << ++center_count << setw(11) 
		    << elem.charge() << buf << endl;
	basis_list << setw(3) << center_count << " 0" << endl;
	pb -> print_gaussian_format(basis_list,no_bas,no_prim);
      }      
      continue;
    } 
    cerr << " BasisSet::read token: " << token << endl;
    ioerror(is," BasisSet::Illegal Token Encountered. ");
  }

  center_list << " ----------------------------------------------------------\n";
  int nu = tot_charge / 2;
  int nd = tot_charge - nu;
  
  basis_list 
    << "   " << setw(3) << no_bas << " basis functions      " 
    << setw(3) << no_prim << " primitive gaussians\n"
    << "    " << setw(2) << nu << " alpha electrons       " 
    << setw(2) << nd << " beta electrons\n\n";
}

/*TEX
\subsection{Output of PureBasisSets}
*/  		  
ostream& operator<<(ostream& os,PureBasisSet& set)
{
  for(int i = 0;  i < set.size(); i++)
    os  << *(set[i]) << endl;
  return os << endl;
  
}
/*TEX
\subsection{Output of PureBasisFunctions}
*/  		  
/* OLD CODE -- DEAD

ostream& operator<<(ostream& os,const PureBasisFunction& bf)
{
  os << "     --> Contracted Function with Momentum: " << bf.momentum()
     << " Degree of Contraction: " << bf.contraction() << endl
     << "         Exponent       Coefficient " << endl;

  char buf[80];
  for(int i=0; i < bf.contraction(); i++)
  {
    sprintf(buf,"    %15.8f %15.8f ",bf.exponent(i),bf.coefficient(i));
    os << buf << endl;
  }  
  os << " +++ ------------------------------------------------------------------ +++\n";
  return os;  
}
*/
