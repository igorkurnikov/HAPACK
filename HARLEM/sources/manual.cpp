//  manual.cpp
//
//  Main Pages of HARLEM Manual
//
//  Igor Kurnikov, Northwestern University 
//
//  Created:   April 9 2002
//

/*! \mainpage HARLEM Programmers Manual
  
  \section sec_abstract Abstract

  HARLEM (HAmiltonians to Research LargE Molecules) was initially developed 
  as a program to compute electronic donor-acceptor interactions using <i>ab initio</i>
  Divide-and-Conquer approach but now is gradually evolving to become a multipurpose 
  molecular simulation package with some unique capabilities not to be found 
  in other programs.
  
  Currently the program implements algorithms to compute parameters of 
  long-distance Electron Transfer reactions (donor/acceptor interaction, reorganization energies)
  in biological and non-biological systems. It also allows (alone and in a combinationation with other programs) 
  to perform Molecular Mechanics, Quantum Chemical and continuum electrostatics calculations of 
  macromolecular systems.

  The program has an object oriented design and is written mostly in C++ programming language.
  I also uses some functions and libraries written in C and FORTRAN.
  Different computational tasks are performed by separate program modules, organized as
  C++ classes. This permits easy expandability of the program functionality. The program
  has a Graphical User Interface (GUI) based on MFC(Microsoft Foundation Classes)  and 
  Scripting language interface based on PYTHON. 
  It runs on WINDOWS, LINUX and Digital UNIX ( unix versions run without GUI or 
  using commercial WIN32 porting libraries).

  \ref page_motivations

  \ref page_install

  \ref page_prog_struct

  \ref page_gui

  \ref page_et_calc

  \ref page_qchem

*/
  
/*! \page page_motivations Motivations

  In the design of computational algorithms for quantum chemical problems a very important
and probably most time-consuming part is the design of a computer code, which implement
this algorithm. A wide use of ab initio and semiempirical methods for electronic struc-
ture calculations in a great part due to the availability of su.ciently convinient Quantum
Chemical programs. "Black Box" program packages GAUSSIAN, SPARTAN, GAMESS,
MOPAC, HYPERCHEM and others allowed increasingly wide community of scientists to
get useful theoretical insights without bothering about fine details of computational algo-
rithms and without spending too much time for learning the program operation.
A creation of a quantum chemical package require very substantial programming efforts.For
example a popular GAUSSIAN package has been developed for more than 30 years with
more than 50 developers contributed to the code. While designing new computational
algorithms it is reasonable to try to maximally reuse the code already existed and to
implement a new method as a module of a some large package to join to the efforts of
other researchers and developers. However for several reasons we decided for our purposes
to develop our own program package while still trying maximally integrate our code with
existing molecular simulation software. To indicate our main research target
properties of macromolecules we used the acronym HARLEM (HAmiltonians to Research
LargE Molecules) for the name of our program.
Our goals in creating HARLEM package were the following:

<ol>
<li> We aim to overcome shortcomings of the architecture of the existing quantum chemical
software. The divide-and-conquer methods for electronic structure calculations we develop
differ substantially from the standard quantum chemical methods which make the popular
quantum chemical packages like GAUSSIAN to be inadequate for their implementation.
For example the concepts of the bond, chemical group are actively used in the fragmen-
tation algorithms while GAUSSIAN consider molecules just as a set of atoms placed in
the space. Thus we wanted our program to be "more smart" about the chemical concepts,
which require a more sophisticated organization of the structural data and the methods to
handle them then is typically done in existing quantum chemisry software.

<li> We want to achieve a high degree of the program usability, so it can be a valuable to a
wide range of scientists, not neccessary experts in the field of theoretical chemistry. For this
purpose the program was designed with a sophisticated Graphical User Interface (GUI).
GUI permits fast learning of the program functionality, reduces the time to setup calcula-
tions and helps to avoid errors in setup of input parameters. Graphical presentation
of the results of the calculations in many cases make them much easier to understand and
often give valuable insights that is usually a goal of a scientific study.

<li> The program should have a modular structure which permit easy modifications and
extensibility of the program functionality. Many developers with different expertize and
goals may add and modify the code. The program architecture should be extensible enough to
allow gradual evolution of the program, so the different changes to the code will typically
not effect each other. Also chemical and biological problems are rarely can be classified as
pure electronic. Usually one would need to perform many different types of calculations
such as quantum chemical calculations and molecular dynamics simulations for a compre-
hesive theoretical analysis of the system. Therefore we target in the program design to put
"under the same roof" different types of molecular modeling. All these properties we tried
to achieve by an object-oriented design of the program.

<li> The methods for theoretical analysis of chemical and biological problems proliferate
very rapidly and the functionality of the program we develop will always cover only a
small fraction of the possible theoretical approaches. We tried therefore to insure that our
program will interact smoothly with other modeling software, permiting easy exchange of
the data and perform a particular task on the particular software package which does it
best

</ol>

  \subsection subsection1 The first subsection
  Text.
  \subsection subsection2 The second subsection

*/

/*! \page page_gui Graphical User Interface

<img border="0" src="../images/ScreenHunter_001.gif" width="600" height="500">

<p> The Graphical User Interface(GUI) is important part of the program which very much
simplify the work with a complicated program and permit intuitive dislpay of results of
simulations. GUI designing require substantial programming eorts and until only com-
ercial simulation packages were based on the GUI. Our nal goal is to achieve the same
type GUI functionality and 
exibilty as in the best comercial programs as InsightII or
HYPERCHEM, while maintain the code for developers so new functionality can be added
to the program.
GUI of HARLEM is build on the basis of MFC(Microsoft Foundation Classes) which is a
standard library for the development of applications for WINDOWS operating system, X-
WINDOWS MOTIF library (a standard for UNIX environment) and RASMOL a powerful
and portable freeware program for molecular visualization. HARLEM adopts a so called
Document/View architecture, which is supported by MFC classes. MolSet class which
as we described in the previous section the core class of the program specing the geometry
molecular system under study is derived from the CDocument class of the MFC library,
while 3D representation of the molecular system is given by HaMolView class which has
a parent CView class of MFC library. We employ so called MDI (Multiple Document
Interface) organization of the program, which immediately permits several instances of
CDocument class (in this case several molecular sets) to coexist in the same instance of
the program.
HaMolView class which shows a 3D image of the molecule is built on the basis of RASMOL
program. RASMOL allow several methods of 3D representation of the molecules. Available
display modes include: (1)wireframe model (bonds are represented as lines) (2) spacell
(atoms are represented as spheres) (3) sticks (bonds represented as cylinders) (4) ball and
stick (a combination of sticks model and spacell model with decreased atom radii) (5)
ribbons and cartoon(solid ribbon) represenation of the biomolecular secondary structure.
RASMOL permits 
exible coloring schemes for atoms and bonds, allows the display of
119
the atomic labels, boundary boxes, dot surfaces, hydrogen bonds (as dotted lines). It has
simple menu interface for the most frequently used functions and a quite developed text
command language. RASMOL can read and write molecular geometry in several popular
formats and save displayed molecular image is several graphical formats.
While incorporating RASMOL to our program we tried to preserve RASMOL functionality
and compatibility of the text command interface. Thus a person familiar with RASMOL
can immediately perform basic molecular display operations on HARLEM.
RASMOL is written in C programming language and because of the way the data are
organized inside the program allows only one molecule and one molecular image to be
displayed. We translated RASMOL code into C++, incapsulating into separate classes
data structures and code of the program performing specic tasks and eliminating global
functions and variables. We achieved a substantial localization of the code of the program,
which among other allowed to present simultaneusly several molecules in the same window,
each of each can be indepedently created and manipulated. This functionality will be very
useful for future simulations of intermolecular interactions with the program, and the
addition of molecular editing functionality to the program.
The method of allocation of memory and storage of data in HARLEM substantially changed
compare to RASMOL. In RASMOL the collections of item such as atoms are organized
as linked lists, with each item containing a pointer of the next item in the collection.
Allocation of the memory for items were performed with specialized utility functions. We
straighforward memory allocation and storage of the data using contaner template classed
of the STL library.
RASMOL is running on UNIX, WINDOWS and MAC-OS operationg systems with a min-
imal amount of a code specific for a given platform. Such a transferability is achieved
thanks to the development of the specialized library of functions realizing 3D primitives,
such as line, cylinder of sphere. These functions ll a buer which represent 3D image of
the molecular system. This buer can be displayed on the screen, with just a few system
dependent commands, which makes most of the code of the program system independent.
We incapsulated graphical primitives and the associated data structures into a class Can-
vas3D which is contained in the main molecular display class HaMolView. Substitution of
the Canvas3D class by the class which maintain the same member function call interface
but using different libraries of 3D primitives for example OpenGL library, provide a path
120
to modify molecular 3D display in HARLEM without great changes in the program text.
HaMolView class maintain a pointer to an instance MolSet class which contains a list
of pointers to 3D objects being displayed. An C++ abstruction for the 3D object is
Object3D class. Every 3D object can be rotated and translated independently of each
other in screen coordinates (only tranformation matricies of the objects are chanaged) or
in real space (reference coordinates of the object are changed). There are currently only
two types of 3D objects : HaMolecule class (abstraction for a molecule) and HaSurface class
(abstruction for a contour surface). However new graphical elements can be easily if they
can be displayed with graphical primitives of Canvas3D class. The possibility to display
contour surfaces is a new functionality of HARLEM compare to RASMOL. Introduction
of these new graphical elements was very much simplied by the object-oriented design of
the program.
Interaction of a user with the program functional modules is mostly accomplished through
the dialogs. In the Windows version of HARLEM every dialog is a class, deirved from
MFC CDialog class. Dialogs of the Unix-Motif version of the program are organized as
functions because Motif library written in C not in C++. There about 20 dialogs currently
in the program. Every dialog maintains a pointer to the class which is associated with, for
example a pointer of some computational module (an instance HaCompMod class). When
the user click mouse close to a particular atom in the 3D display window associated with
HaMolView class, open dialogs which can use this information are notied. This way the
user can interactively choose atoms or other graphical elements in the 3D display.
Unix/Motif version of the program was not updated for several mothes due to the neccessety
to rewrite dialogs in Motif, which were relatively easily generated under Windows using
MFC classes. We plan to reorganize Unix Menu and Dialog interface so it will be more
compatible with Windows interface, eliminating or reducing the need to write the different
code of the same functionality for different platforms. One path can be the use of comercial
or freeware WIN32 clone libraries, which allow to compile MFC library on UNIX and
directly transfer the code from Windows to UNIX. Another alternative can be to write our
own clone of the main MFC classes with a minimum neccessary functionality to encapsulate
MOTIF functions into classes with the same member function calling specications as MFC
classes.
RASMOL mechanism of text command processing was modified in HARLEM to better
fit the object-oriented design of the program. Classes of the program which need handle
text commands are derived from HaTextCmdTarget class. A command entered through
dialog interface a by other means (for example form the processing of the script or remote
command) First handled by currently active instance of MolSet class which then redirect
it if neccessary to be processed in other class derived from HaTextCmdTarget for example
some computational module. Such code organization permits developers working on a
particular module introduce text commands specific for their modules without changing
basic classes used by other modules.

*/

/*! \page page_install Program Installation
 
   Subsections: \ref win_install and \ref unix_install
  
    \section win_install Windows Installation

  <ol>
       <li> Download harlem.zip file from <a href=http://kurnikov.org/harlem_download/>HARLEM download page </a>
	   <li> Unpack content of harlem.zip file into c:\harlem directory, 
	        this should create 
			<ul>
                 <li> HARLEM executable harlem.exe and several dll files in the directory c:\harlem
				 <li> Subdirectory c:\harlem\residues_db with residues and force-field databases 
				      (currently files amber_94_ff.dat, aminoacids.hlm,cofactors.hlm.nucleotides.hlm,water.hlm)
				 <li> Subdirectory c:\harlem\scripts with Python scripts
				 <li> Subdirectory c:\harlem\basis with basis sets in DALTON format
				 <li> Subdirectoty c:\harlem\harlem_manual containing HARLEM documentation.
				 <li> Subdirectory c:\harlem\examples with a few input files and scripts to run HARLEM jobs
            </ul>
	   <li> If installation has been performed into the directory other than c:\harlem (for example d:\harlem_prog\)
	        please set environment variable HARLEM_HOME to this directory (SET HARLEM_HOME=d:\harlem_prog\)
	   <li> Include the directory c:\harlem into your PATH enivronment variable, also add to path the directory 
	        containing Internet Explorer binary iexpore.exe, to be able to invoke manual from inside the program 
   </ol>
    \section unix_install Unix and Linux Installation
   
*/

/*! \page page_prog_struct Program Structure

HARLEM has an object-oriented design and is mostly written in C++ language. We
target the programs for the computers running Windows and UNIX operating systems.
During the last ten years the C++ standard stabilized and good compilers appear for all
popular computer platforms. A lot of useful mathematical software exist which is written
116
in FORTRAN, also a FORTRAN function being compiled will typically run faster than an
analagous C or C++ function. Typically basic mathematical data transformations such
as matrix multiplication or diagonalization is done within the program by FORTRAN
functions, encapsulated into a C++ member function of an appropriate class (for example
matrix class).
The central class in HARLEM is the MolSet class (see Figure ?? which is an abstraction
of a set of molecules. There can be several instances of MolSet classes in a running
instantance of the program (HarlemApp class) (for example a molecular system under study
and its fragments or several unrelated molecular systems). Molecular sets (instances of
MolSet class) can be dynamically created and destroyed (for example when one performs
a subdivision of the molecular system under study into molecular fragments). MolSet
class maintains as a list of pointers to fragments, which are also the instances of MolSet
class.
Molecular set (MolSet class) contain several molecules (vector of pointers to instances
of HaMolecule class). Molecules can be independently created (for example being loaded
from the molecular coordinate file) or manipulated (for example rotated and translated in
3D space).
HaMolecule class contain lists of atoms(Instances of HaAtom class), covalent bonds(Instances
of HaBond class) and hydrogen bonds( instances of HaHBond class). In the HaMolecule
class atoms are organized into residues (instances of HaResidue class corresponding for
example to aminoacids in proteins), which in turn belong to chains (instances of HaChain
class). HaMolecule class also maintains the list of chemical groups corresponging to chem-
ical functional groups such as peptide groups of methyl groups (instances of HaChemGroup
class), which provides an alternative classification of atoms in the molecule.
HARLEM makes extensive use of STL (Standard Template Library of C++). Thus col-
lections of the items usually represented by some STL data structure, like list, vector, set,
map or others. Cycles through the collections of the objects such as a set of atoms in the
molecule is done through the standard incrementation of the iterators of the corresponding
STL containers or with the help of simple utility member functions such as GetFirstAtom()
and GetNextAtom() which manipulate these iterators.
MolSet class has member functions which read and write atom coordinates, list of bonds
and other molecular info in the internal format of HARLEM (denoted by the extention
117
"*.hlm") or in some standard molecular description formats such as PDB (Protein Data
Bank) format.
MolSet provides a basic description of the molecular system under study and is assumed
to be "conservative" so changes in the interface of the class as provided by access member
functions will change only infrequently, so other program modules which use the MolSet
will not need to be updated accordingly too often.
Useful simulation work in HARLEM is done in computational modules which are realized
as classes inherited from HaCompMod base class. MolSet class maintains a list of
pointers of instances of HaCompMod class associated with a given molecular set. There
are currently 5 types of the computational modules in HARLEM: 

1. HaQCMod - class
which setup quantum chemical model of the system such as basis set, electronic wave
function type etc.
2. ETCouplMod - class which responsible for donor/acceptor electronic coupling calculation.
3. HaGaussMod - class which manages interaction with GAUSSIAN quantum chemical pack-
age. 
4. HaDaltonMod - class to setup quantum chemical calculations through the DALTON quantum
chemical program.
5. ElectrostMod - class to the numerical solution of 3D Poisson-Boltzmann.
6. HaInterMolMod - class to compute intermolecular interaction energy and perform dockin simulations
                   including computations of bimolecular ET rates
7. HaMolMechMod  - class to perform molecular mechanics calculations ( uses AMBER as an external program)
8. HaScatterMod  - class to simulate electronic scattering experiments 
9. StmMod        - class to simulate STM (Scanning Tunneling Microscopy) experiments
10. NuclAcidMod   - Nucleic Acid modeling module 

The computational modules are created and desroyed dynamically whenever they needed.
Every Computational Module recieves a pointer to the parent instance of MolSet class.
Thus molecular geometry and structural informations became immidiately available for a
created computational module. Computational modules can easily create other modules
and access information in them through the member function GetCompModule() of the
MolSet class. For example when a user initially load a geometry of an electron transfer
system no computational modules exist. When the user start to define donor and acceptor
and perform PATHWAYS calculations, an instance of ETCouplMod class (donor/acceptor cou-
pling module) is created. When the user start to setup Devide-and-Conquer calculation the
instance of HaQCMod class with associated data structures is created. Such a modular
organization permits an easy extension of the program functionality probably by differ-
ent developers working simulatneously but independently of each other. Thus HARLEM
provides a convinient environment for a creation of molecular simulation software.

*/
