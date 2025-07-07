/*TEX
%
% QCHEM - IPACK
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: ipack.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\chapter{IPACK}
\section{Introduction}

\subsection{Goals}

The generation of hamiltonian matrix elements is one of the most basic
prerequisites of any quantum chemistry calculation. The goal of this
package is to provide a {\bf simple, documented and adaptable} tools
to perform basic tasks in the computationally intensive field of
quantum chemistry. 

The most important goal of this program was to implement a coherent
simple hierarchical structure to perform this calculation. My hope was
that such a structure would provide the highest likelihood that the
program will find wide use in the development of new quantum chemistry
methods. Given this goal, efficiency considerations were secondary,
when the clarity of the underlying concepts and methods was at
stake. The formal developments used in this program were developed in:

\bc 
\begin{tabular}{ll} OS & S. Obara \& A.Saika :
Computation of Molecular Integrals, \\ & JCP {\bf 84} 3963,(1986) \\
HGP & Head-Gordon \& Pople: 
\end{tabular} 
\ec 
which should be consulted by those who wish to use and modify this
package and which offer a good introduction into the subject.

This program should be viewed as a toolkit, which supplies a set of
building blocks used to evaluate the matrix elements needed in quantum
chemistry. Using these blocks I have generated one program, which uses
them in the most conventional sense to compute the ``standard'' matrix
elements calculations for a given basis set. This however, I consider
the minimal application of the program. I hope that the package is
written and structured in such a way that the user can {\bf modify} is
to suit the individual needs. This means first and foremost to alter
the conventions for input and output, but also more fundamental
aspects of the integral generation.  Of course it is impossible to
anticipate all possible derivations one contemplate from a
``standard'' matrix element implementation. However, I hope that this
program was written (and documented) in a sufficiently modular and
transparent fashion that ``high-level'' modifications can be
accomplished by reconfiguring the high-level objects in this package,
while low-level modifications should derivable locally in the
(hopefully) sufficiently well documented routines. This particularly
pertains to the adoption of the low-level routines which the user
might wish to optimize for specific architectures.

As the reader has probably already guessed from my use of codewords
like ``object'' and ``derivation'' this package was written using C++
and not in good old FORTRAN. Those who whish to modify it will need
some command of the language, while those who merely want to run it
need an ATT cfront 3.0 compatible C++ compiler.

\subsection{Contents}

The package presently allows the evaluation of
\begin{itemize}
\item one-electron overlap integrals
\item kinetic energy integrals
\item nuclear repulsion integrals
\item electric field integrals
\item electric field gradient integrals and
\item two-electron coulomb repulsion integrals
\end{itemize}
It can be adopted to do compute a number of other
single-particle integrals:
\begin{itemize}
\item magnetic moment integrals,
\item angular momentum integrals and
\item spin-orbit integrals
\end{itemize}
I think the computation of forces is equally feasable, but presently
not incorporated. 

Presently the biggest drawback is the lack of implementaion of
symmetries, a problem which we hope to be able to fix soon. However,
since symmetries can (at a high cost) always be introduced through a
transformation of unsymmetrized integrals, this is more an issue of
effciency than of principle.

\subsection{LICENSE} 

I am not really an expert in dealing with this part, however I would
like to excercise some kind of control regarding the distribution and
spread of this package. I hope that this program proves to be useful
tool and will be widely used, it is, however, not in the public
domain, since it it involved a substantial effort to produce. My main
intent of the following LICENSE is to keep track of who is using this
program and for what purpose. In general, at least for academic users,
you will be able obtain a license to use this program free of charge
and install it on as many computers as you like as long as you insure
that only people which are covered under your licensing agreement can
use it.

Upon request you may further obtain a license to attach this package
to whatever application you intend to develop and obtain the right
to grant licenses for this package to the users of your application as
provided in the licensing agreement. Should you intend to, you may
further be granted a license to modify this program, which is afterall
one goal for its development, and distribute the altered version under
an appropriate licensing agreement. If you inquire for a licensing
agreement, you should therefore state clearly your intent as to how you
use this package in order to keep the beaurocratic hassle at a minimum.

For these reasons the IPACK integral package is distributed under
the following

\bc
{ \bf \Large LICENSE }
\ec
\begin{itemize} 
\item[0.]  The following is the SOFTWARE LICENSE for the IPACK Integral Package, 
henceforth referred to as the "program" written by
\bc
 Wolfgang Wenzel \\
 Theoretical Physics I \\
 Dortmund University  \\
 D-44221 Dortmund \\
 Germany \\
 e-mail: wenzel@cip.physik.uni-dortmund.de 
\ec

This program and its attached documentation is propriatary software,
it is therefore illegal to use, modify or distribute this program
without the WRITTEN PERMISSION in form of an LICENSING AGREEMENT of
the author or a person expressely authorized to do so. To obtain such
permission contact the author or authorized person under the address
above with a request for a license, stating the purpose and
circumstances for which you wish to use it. The author or authorized
person will then send you a LICENSING AGREEMENT stating the terms and
conditions under which you may use/modify and distribute this program,
which becomes part of this LICENSE.

\item[1.]  You may not copy, modify, sublicense, link with, or
distribute the program except as expressly provided under this
License.  Any attempt otherwise to copy, modify, sublicense, link
with, or distribute the program is void, and will automatically
terminate your rights under this License. 

\item[2.]  You are not required to accept this License, since you have
not signed it. However, nothing else grants you permission to use or
modify the program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
using or modifying this program (or any work based on this program),
you indicate your acceptance of this License to do so, and all its
terms and conditions for copying, distributing or modifying this
program or works based on it.

\bc
			    NO WARRANTY
\ec
\item[3.] THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT
PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING
THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS"
WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY
AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE
DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
CORRECTION.

\item[4.]  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN
WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY
AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU
FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.
\end{itemize}

\section{HOWTO}

\subsection{How to obtain IPACK}
Read the LICENSE Section first. Then e-mail:
\bc
     wenzel@cip.physik.uni-dortmund.de
\ec
with a request for a license. 

I think there are three ``typical users'' of this package and have
supplied a brief HOWTO guide for each.

\begin{itemize}
\item Everyone interested should need a copy  of the program and the
first section, ``How to Obtain IPACK'' is therefore required reading
for all. 
\item If you wish to use the program as is, the reader needs to read 
hardly any of the documentation to this program and 
none of the code. The section ``How to Use IPACK'' is for you and
should contain all you want to know about this package. 
\item If you wish to optimize IPACK for a parcticular computer
architecture, consult the section ``How To Optimize IPACK''. Again you
will be able to skip most of this document and focus only on the
elementary operation section.
\item If you wish to modify the input and output conventions, work
with multiply basis sets at the same time and/or wish to compute only
subsets of the integrals the section ``How to Modify Top-Level Objects
in IPACK'' is for you. 
\item If you wish to implement some deep revolutionary change the way
integrals are computed or go beyond the previous three points you will
probalby have to work through a good part of this manual. The section
``How to Read the IPACK Document'' might help.
\end{itemize}

\subsection{How to use IPACK}
Although running this program does not require any familiarity with
either the code or the language it is written in, some concepts
developed in the documentation are nevertheless required. Since one
goal of this package was to supply self-documenting code (which is
supposedly easy to write for the programmer) there is no auxialiary
manual and the reader will have to work through some of the sections
of the documentation to understand the concepts involved in the
formulation of the problem. In order to run IPACK, you should:
\begin{itemize}
\item Read the paper by Obara \& Saika and the paper by Head-Gordon
\& Pople.
\item Read the following sections of documentation, while skipping the
attached code:
\begin{itemize}
\item the chapter~\ref{chapter-basis} entitled: Basis Sets, skipping the subsections on
``Internal Variables'', ``Class Header'' and ``Preparing Primitive 
Information''.
\item the section~\ref{section-output}: ``Output of Two-Electron Integrals''
\item the section~\ref{twoel-concepts} and the
section~\ref{twoel-impl} on ``Computation of Two-Electron Repulsion
Integrals'' 
\item the section~\ref{ipack-main}: ``IPACK - Main Program''
\end{itemize} 
\item Look at the sample basis files `sample1' supplied with this
program and try to run it. 
\end{itemize}

\subsection{How to optimize IPACK}

We have tried to express the computationally intensive steps of
this package in a small set of routines, which are contained in the
described in the section\ref{chapter-storage} ``Storage'', section\ref{section-f-eval} 
``Evaluation of F'' and~\ref{section-contract} ``Contractions''.
 
\subsection{How to modify Top-Level Objects in IPACK\label{howto-upper-level}}

It is quite straightforward to split the computation of the matrix 
elements into two conceptually distinct, hirarchical levels. The
upper level deals with the admistration of such things as orbitals, 
basis sets, symmetries and input and output. Given two basis sets
the lower level deals with the computation of all the unique 
matrix elements where any function of the first basis set may appear
on the left side of the operator and any basis set of the second set
may appear on the right side of the operator. There are a (neccesarily
limited) number of modifications, which can be performed to adopt this
program for your purpose, which require knowledge and operation only
of the higher and vastly simpler level of this program. Such
adoptations may pertain to, but are hopefully not limited to
\begin{itemize} 
\item Costumization of Input and Ouput
\item Partial Evaluation of Subsets of Integrals
\end{itemize}
If this is what you need to do, you need to read and manipulte 
only a relatively small section of documentation and code of this
package, which is listed here:
\begin{itemize}
\item Learn how to run the program, when reading the information on
basis sets, study the code. 
\item Look at the section~\ref{ipack-main} below. 
\item To deal with high-level manipulations, you will probably only to
modify the code in the main routine of IPACK and the IO classes.
\end{itemize} 

\subsection{How to read the IPACK Document}
Well if none of the above was  sufficient, you have to deal with the
full package. What follows is a suggested order in which to read the
documentation in order to subdivide the code into smaller subsets. It
may be beast to first just read the documentation and the class
headers, leaving implementation details for later.

\begin{itemize}
\item Work through the section~\ref{howto-upper-level} ``How To Modify
Upper-Level Objects in IPACK'' if you have not already done so. This
should supply you with a top-down view of IPACK. At its bottom, you
will find calls to the routines and classes which do the computation,
such as \name{overlap},\name{twoel},\name{kinetic} etc. 
\item These routines, however, use a completely different
representation of the \name{BasisSets} and it may be best to learn
this representation as well as other concepts required in the
implementation of the recursion relations. All of these are introduced
in the chapter~\ref{auxiliary-classes} ``Auxiliary Classes'', which
you should work through completely. 
\item The next step is to understand the concepts involved in the 
implementation of a {\bf single-particle recursion relation}, 
as discussed in the chapter~\ref{chapter-recursion}.
If you have understood the class \name{Recursion} and the associated 
data-types class \name{IntegArray} and \name{Coefficient\_Set} the concepts 
of the implementation of the single-particle recursion relations. The
rest is filling in the blanks, specific operator recursion relations
are discussed in sections~\ref{section-overlap},\ref{section-kinetic}
and \ref{section-nuclear} for instance. 
\item It may surprise you that all the {\bf two-particle recursion
relation}, for the electron repulsion integrals can be forumulated
entirely in terms of the \name{Recursion} class for single particle
operators. The evaluation of these recursion relations is discussed in
chapter~\ref{chapter-twoel}.The required classes, \name{VRR},
\name{HRR1} and \name{HRR2} are derived from class \name{Recursion}.
\end{itemize}

After completing the above program you have been introduced to all
important concepts and classes of this package. If I have succeeded in
my goal to write a documented and transparent program, you should now
have the tools at your disposal to augment this package to suit your
needs. 

\section{IPACK --- Main Program\label{ipack-main}} */
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include "typesip.h"
#include "basis.h"
#include "operators.h"
#include "f.h"
#include "integ_array.h"
#include "integ_file.h"
#include "nuc_rep.h"
#include "nuc_mom.h"

#include "nuc_e_field.h"
#include "vault.h"
#include "sphere.h"
#include "ipack.h"

// Timer* timer;

int main(int argc,char** argv)
{
  pm = new ParMachine(&argc,&argv);

  timer = new Timer;
  init_timer();
#ifdef MEMCHK_IPACK
  memory_manager = new Memory_Manager(100000);
  memory_manager -> on();

  if (argc > 1)
  {
    int allocno;
    sscanf(argv[1],"%i",&allocno);
    cout << "intercepting: " << allocno << endl;
    memory_manager -> intercept(allocno);
  }
#endif
  
  DataVault   vault("Data",vault_new);
  IntegFile** ifile = new IntegFile*;

  ipack_main(cin,ifile,vault,false);

  return 0;
}

/* 
  int debug_flag  = 0;
  timer -> start(t_main);
  timer -> start(t_sp);

  if (!input_skip(cin,"IPACK")){
    if (pm->master())
      cout << "IPACK-Header not found" << endl;
    pm -> abort();
  }

  init_ftable();

  //String ftablename(FTBL_NAME);
  //read_ftable(ftablename);

  
  // superceded by HEIKO's flags spherical    = false;   
  // cartesian gaussians by default

  BasisSet basis;  
  PureBasisSetList purebasislist;
  
  ifstream bstream("BASIS",ios::in);
  bstream >> purebasislist;  

  Molecule molecule;
  MomentData mu_data;
  Location   e_field_loc;
  Location   anglr_loc;
  Location   spin_orb_loc;
  
  SStack<Location,max_cntr>   centers;
  SStack<  double,max_cntr>   charges;
 
  basis.parse(cin,molecule,mu_data,e_field_loc,anglr_loc,spin_orb_loc,
	      purebasislist);

  int core_size = 0;
  int i;
  
  ARRAY<IntegMat> ovlp;
  ARRAY<IntegMat> kin;
  ARRAY<IntegMat> mom;
  ARRAY<IntegMat> nuc;
  ARRAY<IntegMat> e_fld;
  ARRAY<IntegMat> ang;
  ARRAY<IntegMat> sp_orb;

  InternalBasisList iblist;
  OrbInfo           orbinfo; 
  int no_orbs = 0;
  
  for(i=0; i < molecule.size(); i++)
  {
    InternalBasis* ib1 = new InternalBasis(basis,molecule.center(i),ALL_CENTER,
					   orbinfo);
    iblist.add(ib1);
  }

  ARRAY<unsigned int> orbital_map;
  no_orbs  = orbinfo.make_orbital_map(orbital_map);
  
  IntegDiskFile ifile("INTEG",no_orbs);
  ifile.set_orbital_map(orbital_map);
 
  orbinfo.save(vault);
  
  const int buf_size = 2560;
  char  molec_buffer[buf_size];
  // ostrstream  molec_geom(molec_buffer,buf_size);
  molecule.write(vault);

  double nuc_rep = nuclear_repulsion(molecule,core_size);
  vault.insert("E_NUCLEAR",nuc_rep);
  vault.insert("TOTAL_CORE",core_size);

  int is_atom = centers.size() == 1 && 
      (fabs(centers[0](X)) < delta) &&  
      (fabs(centers[0](Y)) < delta) &&  
      (fabs(centers[0](Z)) < delta);

  vault.insert("MOLTYPE","M"); // always assume molecular calculation

  if(flaglist(Overlap))
  {
    overlap(iblist  ,ovlp);
    normalize(iblist,ovlp);
    ifile.write_mat(ovlp[0],IntegOvlp);
  }
  
  if(flaglist(Kinetic))
  {
    kinetic(iblist,kin);
    normalize(iblist,kin);
    ifile.write_mat(kin[0],IntegKinetic);
  }
 
  if(flaglist(Nuclear))
  {
    nuclear(iblist,molecule,nuc);
    normalize(iblist,nuc);
    ifile.write_mat(nuc[0],IntegNuclear);
  }
       
  if(flaglist(Moment))
  {
    nuclei_moment(molecule,mu_data,vault);
    moment(iblist,mu_data,mom);
    normalize(iblist,mom,mu_data.max_mu);

    vault.insert("MAX_MOMENT",mom.size());
    String name="MOMENT";
    for(i = 1; i < mom.size(); i++)
    {
      ifile.write_mat(mom[i],IntegMoments + i - 1);
    }
  }

  if(flaglist(eField))
  {
    nuclei_e_field(molecule,e_field_loc,vault);
    e_field(iblist,e_field_loc,e_fld);
    normalize(iblist,e_fld,2);

    for(i = 1; i < e_fld.size(); i++)
      ifile.write_mat(e_fld[i],IntegEField+i-1);
  }
  
  if(flaglist(Angular))
  {
    angular(iblist,anglr_loc,ang);
    normalize(iblist,ang,1);

    for(i = 1; i < ang.size(); i++)
      ifile.write_mat(ang[i],IntegAngular + i - 1);
  }

  if(flaglist(SpinOrb))
  {
    spinorb(iblist,spin_orb_loc,sp_orb);
    normalize(iblist,sp_orb,1);

    String name="SPINORB";
    for(i = 1; i < sp_orb.size(); i++)
      ifile.write_mat(sp_orb[i],IntegSpinOrb + i - 1);
  }

  timer -> stop(t_sp);
 
  twoel(iblist,ifile);

  ifile.write_all();
  
  timer -> stop(t_main);
  if (debug_flag)
    timer -> print();
  delete timer;
  timer = 0;

  
  vault.insert("STATUS","I");  
  clean_ftable();
  Element::clean();
  Sphere::transform_clear();
  Storage::clear_work();
  
  if (debug_flag)
    Storage::print_stats();    
}
#ifdef MEMCHK_IPACK
  memory_manager -> print(cout);
#endif
  return 0;

}

 
*/
