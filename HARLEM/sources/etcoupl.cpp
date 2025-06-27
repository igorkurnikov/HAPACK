/*!  \file etcoupl.cpp
 
  Classes to calculate Electron Transfer coupling in HARLEM.

  \author  Igor Kurnikov , University of Pittsburgh 
  \date 1998-2002

*/

#define ETCOUPL_CPP

#include <mpi.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "g94_globals.h"
#include "haio.h"
#include "halinalg.h"
#include "etcoupl.h"
#include "hamatdb.h"
#include "halocorb.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "harlemapp.h"
#include "hamolview.h"
#include "hamultipole.h"
#include "haintengine.h"
#include "math_num.h"
#include "haqchem.h"
#include "moleditor.h"
#include "hasurface.h"
#include "tokens.h" 
#include "tinyxml.h"


/*! \page page_et_calc Electron Transfer Rates calculations

  Subsections: \ref pathways_calc and \ref qchem_da_calc
 
The main module (C++ class) which perform calculations of donor/acceptor coupling
is ETCouplMod class. As all computational modules it is derived from HaCompMod and
HaTextCmdTarget classes.
There are three dialogs associated with ETCouplMod:
ChooseRedox - dialog to interactively choose donor and acceptor groups,
PathwaysDlg - dialog to perform PATHWAYS calculations for donor/acceptor coupling,
ETEffHamDlg - dialog class to manage calculations of donor/acceptor electronic coupling
with divide-and-conquer method.
ETCouplMod class performs PATHWAYS calculations of three types:
1.BestPath calculations finding the best coupling between donor and acceptor atoms using
Diijkstra algorithm
2. Coupling Map calculations when PATHWAYS couplings from the donor to all other
atoms of the system are calculated.
3. Identification of atoms which belong to the PATHWAYS with coupling values within
certain threshold of the best path. For this purpose coupling from the donor and acceptor
to all other atoms in the system are calculated. An atoms is selected if the product of its
couplings to the donor and acceptor exceed a threshold value.
Divide-and-Conquer calculation of the donor/acceptor coupling are performed in ETCouplMod
class and ETHeffHF class, the letter is an abstraction for effective hamiltonian in Hartree-
Fock approximation. An instance of ETHeffHF class is contained in the ETCouplMod class. To
calculate effective hamiltonian on active orbitals first Green function matrix on the active
orbitals is calculated given the molecular orbitlas and molecular energies of the system.
This matrix is then inverted to get an effective hamiltonian matrix. The effective hamil-
tonian matrix is logically disected into interaction submatricies between active orbitals of
different atomic groups.
Submatricies of the effective hamiltonian are stored in a database file using functions of
Berkeley Database library encapsulated in the class HaMatDB. When storing a group-
group submatrix of the effective hamiltonian, the program check if the submatrix between
group with the same IDs already exist in the database (from the calulations on another
overlapping fragment of the system). If such a matrix exists, the products of protection
124
factors of the groups for the submatrix stored in the database and for the submatrix about
to be saved are calculated. The new submatrix is stored if the product of the protection
factors of its atomic groups are larger then that of the stored submatrix. The less the
protection factor of an atomic group the more we expect the electronic structure of the
group is perturbed in the fragment compare to the intact structure of the system.
The effective hamiltonian of the system accumulated in the database is restored and
donor/acceptor coupling is calculated. Two methods of the donor/acceptor electronic
coupling calculations are currently supported. In the first method the energy splitting of
eigen states of the effective hamiltonian is minimized by the application of the electrical
field between donor and acceptor. This acomplished by
1. The matricies of components of electrical dipole operator with an origin in the midpoint
between donor and acceptor are calculated using class HaOperR.
2. Dipole operator matricies are premultipled by cosines of the direction from the donor
to acceptor and the value of the strength of the electrical field.
3. The energy splitting between donor and acceptor localized eigenstates is minimized
using conjugate gradient method, changing the strength of the applied electrical field and
diagonalizing the matrix of the effective hamiltonian. The eignestates localized in the donor
and acceptor can be chosen interactively using 3D contour representation of isodensity
levels of the eigenstates.
The second method of donor/acceptor coupling calculations, calculate Green-Function ma-
trix elements between donor and acceptor localized states and described in the chapter ??.
The user interactively chooses donor/acceptor localized eigenstates of the eective hamil-
tonian. The program truncates these states to be localized only on the donor or acceptor
atoms to obtain localized donor and acceptor states. Green-function matrix elements be-
tween these states are calculated and the obtained 2x2 matrix is inverted to give eective
donor/acceptor interaction. The second method is more than an order of magnitude faster
than the method based on energy splitting calculations. In the nearest future we will
modify the program to derive approximate localized donor and acceptor states from the
calculations on fragments, and will calculate donor/acceptor Green Function elements with
iterative methods instead of using the inversion of the whole eective hamiltonian matrix,
which should the calculations dramatically faster for very large ET system.
Outer sphere reorganization energies in HARLEM are calculated with the help of 
ElectrostMod class which manages electrostatic tinteractions solving 3D
Poisson-Boltzmann equation. These types of calculations are accessed through the dialog
class ElectrostDlg. The program setup two Electrostatic calculations. One with dielectric con-
stants of the media set to its satic value (. = 80 for water and . = 4 for the protein. In both
calculations +1e charge uniformly distributed over atoms of the donor and (-1e) charge
distributed over atoms of the acceptor. This is very crude representation of the changes
of the charge distribution if the system upon the electron transfer. However this approxi-
mation doesn't typically cause as great error in the calculated reorganization energy value
while make calculations very easy for the user of the program. The user need just choose
interactively donor and acceptor and run reorganization energy calculations with a press
of a button. The calculations take just seconds for default 65x65x65 grid.   

\subsection pathways_calc PATHWAYS calculations

Example of PATHWAYS calculation
in Gray's Horse Cyt-C HIS-33-Ru-bpy complex:

1. Start HARLEM and load PDB file of the complex: 
 
       Easiest way  probably will be to drop the icon of PDB file 
       (in this case cytc_33_1.pdb) on harlem.exe icon.
	Another method will be to use File->Open Menu in HARLEM
        (one should choose PDB file option in Open Dialog).

2. Specify donor and acceptor in Edit Redox Centers Dialog:
	
	a)Open dialog by ET->'Edit Redox Centers'
        b)choose DONOR or ACCEPTOR  toggle buttons.
	c)Fill list of atoms of donor and acceptor by clicking on the 
          desired atom and pressing INCLUDE button. Selected atoms
	  are seen as small spheres.
          One also could choose the whole residue by choosing selection 
          level RESIDUE toggle button and clicking on the desired residue 
          (one can also just type in the prompt). 
	  In the example case of cyt one may wish to choose donor 
          as HEM105 and acceptor RBP200(Ru-bpy). 
	  first one can choose whole HEM105 as a donor
          ( choose selection level residue,
          click on some HEM atom and press INCLUDE),
          then when you press ACCEPTOR toggle button and back DONOR toggle 
          button list of donor atoms (all atoms belonging to HEM group will
          appear). One may delete from the donor atoms  of 
          propionates and methyl groups by selecting them in the list 
          and pressing DELETE button. 

	  All atoms of the donor(acceptor) assumed coupled with the value 1.0
          (short-circuited).

3. Run Pathways:

Choose ET->"Run Pathways" which will open PATHWAYS Dialog.
  
  One can do three types of calculations currently:
  
3 a) Best Path calculations( should take about 30s on Pentium-133 IBM-PC)

   Choose BEST_PATH toggle button and press RUN 
   Description of the BEST PATH between specified Donor and 
   acceptor will appear in HARLEM CONSOLE Window. 
   Best Path member atoms will be denoted by spheres on the molecular
   image.  
   
3 b) Coupling Map: Coupling of all atoms to the donor
   
    Choose Coupling Map option and press RUN button.
    Coupling map will be printed to HARLEM Console Window.
    and Protein will be coloured by coupling values (red- max, blue- min)
    
    If one save the molecule in pdb format
    ( currently one need to open HARLEM command line window and type 
      command at HARLEM> prompt like:
      save pdb cyt33_2.pdb
      
      the last name of the file in arbitrary of course)
   the ln of coupling values will be in the temperature factors of the 
   atoms in the pdb file.

3 c) 3-rd type of PATHWAYS calculations is the selection of all atoms which
     belong to the paths which have coupling vlaued within certain threshold 
    value from the BEST PATH say within a facot of 0.5 or 0.1. 
    These calculations give you an estimate which residues are important 
    for mediating ET coupling between Donor and Acceptor. 

    One runs these calculations by selecting threshold value in the prompt
    and press "SELECT COUPLED" button. (In the case of cyt this will be
    atoms of HIS33, Arg38, Leu23)

    Identifiers of atoms on the important paths will be printed to 
    HARLEM Console Window.  
    Atoms on the important paths and 
    donor and acceptor will be plotted as tubes. One can get idea what atoms
    are important by clicking on selected atoms - their identifiers 
    will appear in HARLEM CONSOLE window.


\subsection qchem_da_calc Quantum Chemical Donor-Acceptor coupling calculations

la-la

*/

ETEdge::ETEdge()
{
	SetDefaultParam();
}

ETEdge::ETEdge(int new_inode1, int new_inode2)
{
	SetDefaultParam();
	inode1=new_inode1;
	inode2=new_inode2;
}

void ETEdge::SetDefaultParam()
{
	inode1=-1;
	inode2=-1;
	coupling=0.0;
//	dist=0.0;
}

bool ETEdge::operator==(const ETEdge& rhs) const
{
	if(inode1 != rhs.inode1) return false;
	if(inode2 != rhs.inode2) return false;
	if(coupling != rhs.coupling) return false;
//	if(dist != rhs.dist) return false;
	return true;
}

bool ETEdge::operator<(const ETEdge& rhs) const
{
	if(inode1 != rhs.inode1) return (inode1 < rhs.inode1);
	if(inode2 != rhs.inode2) return (inode2 < rhs.inode2);
	if(coupling != rhs.coupling) return (coupling < rhs.coupling);
//	if(dist != rhs.dist) return (dist < rhs.dist);
	return false;
}

ETPath::ETPath()
{
	clear();
}

bool ETPath::clear()
{
	if(!trace.empty()) trace.clear();
	coupling=0.0;
	return true;
}

PathStep::PathStep()
{
	coupling = 0.0;
	destination = -1;
	source = -1;
}

PathStep::PathStep(double new_coupling, int new_destination, int new_source)
{
	coupling = new_coupling;
	destination= new_destination;
	source = new_source;

}

PathStep::~PathStep()
{

}
	
bool PathStep::operator < (const PathStep & rhs) const
{
	return( coupling < rhs.coupling);
}

ETCouplMod::ETCouplMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_ET_COUPL,new_phost_mset)
{
	if(new_phost_mset == NULL) 
	{
		std::cerr << " Error in ETCouplMod::ETCouplMod() "<< std::endl;
		std::cerr << " The pointer to the new host molecular set is invalid " << std::endl;
	}
	ptr_qc_mod = NULL;
	ptr_qc_mod = new_phost_mset->GetQCMod(true);
	Clear();
	ham_trunc_type = HAM_S_DIP_TR;
	CopyEigVecsFromMO();
}


ETCouplMod::~ETCouplMod()
{
	Clear();	
}

void ETCouplMod::SetDebugLevel(int new_debug_level)
{
	debug_level = new_debug_level;
}

int ETCouplMod::CopyEigVecsFromMO()
{
	if(ptr_qc_mod == NULL) return FALSE;
	if(ptr_qc_mod->ActBas != &ptr_qc_mod->AtBasis) return FALSE;

	int nb = ptr_qc_mod->ActBas->GetNBfunc();
	int nmo = ptr_qc_mod->MO_coef.num_cols();

	if( ptr_qc_mod->MOene.size() != nmo)
	{
		PrintLog(" Error in ETCouplMod::CopyEigVecsFromMO() \n");
		PrintLog(" Number of MO vectors %d doesn't equal to the number of MO eneries %d \n",
		           nmo, ptr_qc_mod->MOene.size());
		return FALSE;
	}

	enel = ptr_qc_mod->MOene;
	eigv.coef = ptr_qc_mod->MO_coef;
	eigv.bas = ptr_qc_mod->ActBas;

	return TRUE;
}

int ETCouplMod::OnDelAtoms(AtomContainer& del_atoms)
{
	Clear();
	pathways_graph_init_flag = false;
	return TRUE;
}

bool ETCouplMod::Clear()
{
	PrintLog(" ETCouplMod::Clear() \n ");
	ClearPathwaysGraph();

	pathways_calc_type=BEST_PATH;
	nb_dist_limit = 5.0;
	DA_field=0.0;

    pw_nb_decay = 1.7;
	pw_nb_min_dist = 1.4;
	pw_hb_decay = 1.7;
	pw_hb_min_dist = 2.8;
	pw_ln_cov_decay = log(0.6);
	pw_nb_decay_intermol = 0.0;

	m_hbond_paths_flag = true;
	
	use_pert_mat = TRUE;

	best_path_coupl = 0.0;
	log_calc_result = 1;
	tun_ene = -0.2;
	db_file_name="ET_HEFF_HF.hdb";
	

	donor_orbs.Clear();
	acc_orbs.Clear();

	heff_mat.clear();
	ssl.clear();  

	ieig_don.newsize(0);
	ieig_acc.newsize(0);
	eigv.Clear();
	enel.clear();
	set_dab_huck_inter  = False;

    heff_pert_mat.clear(); 
    use_pert_mat = FALSE;          

	don_acc_gf.clear(); 

	el_field_min.clear();  
	da_coupl_val.clear();  

	extern_field.clear();
	DA_dipole.clear(); 

	return true;
}

void ETCouplMod::ClearPathwaysGraph()
{
	if(!nodes.empty()) 
	{
		int nn = nodes.size();
		int i;
		for(i = 0; i < nn; i++)
		{
			if( edges[i] != NULL) 
			{
				delete ((ETEdge*)edges[i]);			
			}
		}
		nodes.clear();
	}
	if(!coupl_map.empty()) coupl_map.clear();
	if(!best_path.empty()) best_path.clear();
	pathways_graph_init_flag = false;
}

bool ETCouplMod::InitiatePathwaysGraph()
{
	if(!phost_mset)
	{
		std::cerr << " ETCouplMod::InitiatePathwaysGraph(): " << std::endl;
		std::cerr << " Host Molecule Set is not defined " << std::endl;
		return false;
	}
	MolEditor* p_mol_editor = phost_mset->GetMolEditor(true);
	p_mol_editor->CalcHBonds(phost_mset);

	ClearPathwaysGraph();

	nodes.reserve(phost_mset->GetNAtoms());

	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if( !(aptr->IsHydrogen()) ) nodes.push_back(aptr);
	}
	
    BoxPartition nodes_part;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	nodes.GetMinMaxCrd( xmin, ymin, zmin, xmax, ymax, zmax);
	
	xmin -= 0.5;
	ymin -= 0.5;
	zmin -= 0.5;
	xmax += 0.5;
	ymax += 0.5;
	zmax += 0.5;

	nodes_part.SetBoundaries(xmin, ymin, zmin, xmax, ymax, zmax);
	nodes_part.DistributePointsToCells(nodes);

	nodes_part.SetRegionRad(nb_dist_limit);

	int nn = nodes.size();
	edges.resize(nn);

	PtrIntMap at_idx;
	int i;
	for( i=0; i < nn; i++)
	{
		void* ptr = (void*) nodes[i];
		at_idx[ptr] = i;
	}
	
	AtomGroup local_atoms;
	
	for(i=0; i < nn; i++)
	{
		HaAtom* aptr = nodes[i];
		std::list<ETEdge>* p_edge_list = new std::list<ETEdge>;
		edges[i] = p_edge_list;
				
		nodes_part.GetNeighbors(*aptr, local_atoms );
		
		int j,nloc;

		nloc = local_atoms.size();
		for( j = 0; j < nloc; j++)
		{
			if( local_atoms[j] == aptr ) continue;

			void* ptr = (void*) local_atoms[j];
			int idx_2 = at_idx[ptr];
			p_edge_list->push_back( ETEdge( i, idx_2));
			ETEdge* pcur_edge;
			pcur_edge=&(p_edge_list->back());
			pcur_edge->coupling = CalcAtomContactCoupling(aptr,local_atoms[j]);
		}
	}

	pathways_graph_init_flag = true;
	return true;
}

static const double LN_06= log(0.6);
static const double LN_036= log(0.36);

bool ETCouplMod::path_coupl_calc()
{
	std::priority_queue< PathStep, std::deque<PathStep>> que;

	char buf[256];
	HaAtom* aptr;
	int i_cur_node;

	int i;

	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if(pathways_calc_type == BEST_PATH && !pacc)
	{
		PrintLog(" Error in ETCouplMod::path_coupl_calc() \n"); 
		PrintLog("Acceptor is not set up \n"); 
		return false;
	}

	HaColor path_col(255,255,0);  // yellow

	best_path_coupl = 0.0;
	if(pathways_calc_type == BEST_PATH) best_path.clear();

	if(!pathways_graph_init_flag) InitiatePathwaysGraph();

	for(i=0; i < nodes.size(); i++)
	{
		if(pdon->HasAtom(nodes[i]) )
		{
			que.push(PathStep(0.0, i, -1));	
		}
	}
	
	if(que.empty())
	{ 
		PrintLog(" ETCouplMod::PathCouplCalc() \n"); 
		PrintLog(" No Donor Nodes are found \n");
		return false;
	}

	int num_nodes=nodes.size();

	if(!coupl_map.empty()) coupl_map.clear();
	if(coupl_map.capacity() != num_nodes) coupl_map.reserve(num_nodes);
	
	for(i=0; i < num_nodes; i++)
		coupl_map.push_back(PathStep(-40,-1,-1));

// Dijkstra's algorithm

	int acc_found=-1;
	while(!que.empty())
	{
		// pull nearest city from the queue
		PathStep cur_step=que.top();
	    int cur_node= cur_step.destination;
		que.pop();

		if(coupl_map[cur_node].destination == -1)
		{
			// add to the largest coupling map

			coupl_map[cur_node]= cur_step;
			HaAtom* pnode= nodes[cur_node];

			if(pathways_calc_type == BEST_PATH)
			{
				if( pacc->HasAtom( pnode ) )
				{
					acc_found=cur_node;
					break;
				}
			}
            
			// put coupling to neighbors to the queue

			std::list<ETEdge>::iterator itr;
			std::list<ETEdge>* p_edges_list = (std::list<ETEdge>*) edges[cur_node];

			for(itr= p_edges_list->begin(); itr != p_edges_list->end(); itr++)
			{
				double edge_coupl = (*itr).coupling;
				double DestCoupling= cur_step.coupling + edge_coupl;
				que.push(PathStep(DestCoupling,(*itr).inode2, cur_node) );
			}
		}
	}
	if(pathways_calc_type == BEST_PATH)
	{
		if(acc_found == -1)
		{
			std::cerr << " Error in ETCouplMod::path_coupl_calc() " << std::endl;
			std::cerr << " No path have been found to acceptor " << std::endl;
			return False;
		}
		best_path_coupl = exp(coupl_map[acc_found].coupling);

		if( log_calc_result > 0 )
		{
			PrintLog("\n\n");
			sprintf(buf, " best PATHWAY coupling value between donor and acceptor is %12.6e ",
			             exp(coupl_map[acc_found].coupling) );

			PrintMessage(buf);
			PrintLog("\n\n");

			PrintLog(" structure of the best path \n");
			PrintLog("ACCEPTOR \n");
		}
		int i_cur_node=acc_found;
		for(;;)
		{

		   if(i_cur_node == -1)
		   {
				if( log_calc_result > 0 )
				{
					PrintLog("DONOR \n");;
				}
				return true;
		   }
		   int isrc= coupl_map[i_cur_node].source;
		   best_path.push_back(coupl_map[i_cur_node]);
		   aptr=nodes[i_cur_node];
		   if( log_calc_result > 0)
		   {
				aptr->FillRef(buf);
				PrintLog("%s  %12.6e  ",buf, exp(coupl_map[i_cur_node].coupling)); 
		   }

		   if(isrc != -1)
		   {
				HaAtom* aptr2=nodes[isrc];
				double step_val= coupl_map[i_cur_node].coupling - coupl_map[isrc].coupling;
				if( log_calc_result > 0)
				{
					aptr2->FillRef(buf);
					PrintLog("%s     %9.3e ",buf, exp(step_val));
				}

				if(aptr->IsBonded(*aptr2))
				{
					if( log_calc_result > 0)
						PrintLog(" COVALENT BOND ");
				}
				else if(m_hbond_paths_flag && phost_mset->AreHBonded(aptr,aptr2))
				{
					if( log_calc_result > 0)
						PrintLog(" HYDROGEN BOND ");
				}
				else
				{
					if( log_calc_result > 0)
						PrintLog(" NONBONDED JUMP ");
				}
		   }
		   if( log_calc_result > 0)
				PrintLog("\n");
		   i_cur_node=isrc;
		}
	}

	if( log_calc_result > 0)
		PrintLog(" Atom Coupling Map \n");

	for(i=0; i < num_nodes; i++)
	{
		i_cur_node= coupl_map[i].destination;
		if(i_cur_node == -1) // Fix Nodes which never were accessed
		{
			i_cur_node=i;
			coupl_map[i].destination=i;
		}
		aptr = nodes[i_cur_node];
		aptr->tempf = coupl_map[i].coupling;
		aptr->FillRef(buf);
	}
	HaBond* bptr;
	// assign coupling to Hydrogen atoms equal to PATHWAYS 
	// coupling to heteroatoms they bonded to
	BondIteratorMolSet bitr(phost_mset);
	for(bptr= bitr.GetFirstBond(); bptr; bptr= bitr.GetNextBond())
	{
		if(bptr->srcatom->GetElemNo() == 1) bptr->srcatom->tempf = bptr->dstatom->tempf;
		if(bptr->dstatom->GetElemNo() == 1) bptr->dstatom->tempf = bptr->srcatom->tempf;
	}

	return true;
	
}

bool ETCouplMod::select_important(double thresh)
{
	if(thresh > 0.9) thresh= 0.9;
	if(thresh <= 1.0e-20) thresh= 1.0e-20;
	thresh=log(thresh);

	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");
	
	if(!pdon || !pacc || pdon->size() == 0 || pacc->size() == 0)
	{
		ErrorInMod(" ETCouplMod::select_important() ",
		           " Donor or acceptor are not set ");
		return false;
	}
	pathways_calc_type=COUPL_MAP;
	path_coupl_calc();

	std::vector<PathStep> coupl_map_save;
	coupl_map_save.swap(coupl_map);

	pdon->swap(*pacc);
	
	pathways_calc_type=BEST_PATH;
	path_coupl_calc();
    PathStep* acc_pstep= &best_path.front();
	double bp_coupl=acc_pstep->coupling;
	double coupl_max=bp_coupl+thresh;
	
	PrintLog("best path coupling is %12.6e \n", bp_coupl);
	PrintLog("Threshold coupling is %12.6e \n",coupl_max);
	PrintLog("List of atoms involved in paths with coupling > threshold: \n");

	pathways_calc_type=COUPL_MAP;
	path_coupl_calc();

	phost_mset->UnSelectAtomsAll();
	int n=nodes.size();
	int i;
	HaAtom* aptr;
	for(i=0; i < n; i++)
	{
		if((coupl_map[i].coupling+coupl_map_save[i].coupling) > coupl_max)
		{
			aptr=nodes[i];
			aptr->Select();
			if(pdon->HasAtom(aptr)) continue; // do not print donor and
			if(pacc->HasAtom(aptr)) continue; // acceptor atoms
		}
	}

	pdon->swap(*pacc);
	coupl_map.clear();

	AtomIteratorAtomGroup aitr_d(pdon);
	for(aptr = aitr_d.GetFirstAtom(); aptr; aptr = aitr_d.GetNextAtom() )
	{
		aptr->Select();
	}

	AtomIteratorAtomGroup aitr_a(pacc);
	for(aptr = aitr_a.GetFirstAtom(); aptr; aptr = aitr_a.GetNextAtom())
	{
		aptr->Select();
	}

	HaMolView* pView= phost_mset->GetActiveMolView();

	// Represent Selected atoms as sticks:
	if(pView) 
	{
		pView->DisableSpacefill();
		if( phost_mset->GetNAtoms() < 256 )
		{  
			pView->EnableWireframe(CylinderFlag,0.16);
		} 
		else 
		{
			pView->EnableWireframe(CylinderFlag,0.32);
		}
		pView->SetRibbonStatus(False,0,0.0);
		pView->DisableBackbone();
		phost_mset->RefreshAllViews();
	}
	
	return true;	
}


double ETCouplMod::CalcAtomContactCoupling(HaAtom* aptr1, HaAtom* aptr2)
{
//!  Set up coupling between nodes in PATHWAYS model 

	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	// short circuit donor and acceptor
	if( (pdon->HasAtom(aptr1) &&  pdon->HasAtom(aptr2)) )
	{
		return 0.0;
	}

	if( pacc && (pacc->HasAtom(aptr1) &&  pacc->HasAtom(aptr2)) )
	{
		return 0.0;
	}

	if(aptr1->IsBonded(*aptr2)) 
	{// Covalent bond coupling
		return pw_ln_cov_decay;
	}

	double ss = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U); 

	double coupling;

	if(m_hbond_paths_flag && phost_mset->AreHBonded(aptr1,aptr2)) // H-Bond coupling
	{
		coupling =  (ss < pw_hb_min_dist) ? 2.0*pw_ln_cov_decay : 
	                    2.0*pw_ln_cov_decay - pw_hb_decay*(ss - pw_hb_min_dist);
	}
	else
	{
	    coupling =  (ss < pw_nb_min_dist) ? pw_ln_cov_decay :                   // non-bonded contact 
	                    pw_ln_cov_decay - pw_nb_decay*(ss - pw_nb_min_dist);
		
		if( pw_nb_decay_intermol > 0.0001)
		{
            if(aptr1->GetHostMol() != aptr2->GetHostMol())
			{
	        	coupling= (ss < pw_nb_min_dist) ? pw_ln_cov_decay : 
	                            pw_ln_cov_decay - pw_nb_decay_intermol*(ss - pw_nb_min_dist);

			}
		}
	}
	return coupling;
}

int ETCouplMod::ColorMolSurfETCoupl()
{	
	MolSet* pmset = GetMolSet();
	if( pmset == NULL)
		return FALSE;

	HaDisplayedSurface* sptr = pmset->CalcMolSurface();
	
	if(sptr == NULL)
		return false;

	int ic;
	int ncol = 5; // the number of colors in the scale

	int nmid = ncol/2 + 1;
	int nspan2 = ncol/2;

	HaVec_int col_num(ncol);

	int r,g,b;

	int col;

	for(ic = 1; ic <= ncol; ic++) //  prepare color shades interpolated between red -> white -> blue
	{
		if(ic <= nmid )
		{
			r = 255; 
			g =  (int) ( 255.01* ( (double)(ic - 1.0)/ (double)nspan2) );
			b =  (int) ( 255.01* ( (double)(ic - 1.0)/ (double)nspan2) );
		}
		else
		{
			r =  (int) ( 255.01* ( (double)(ncol - ic )/ (double)nspan2) ); 
			g =  (int) ( 255.01* ( (double)(ncol - ic )/ (double)nspan2) );
			b = 255;
		}
		col_num(ic) = HaColor::RegisterColor(r,g,b);
	}

	int numverts = sptr->GetNumVerts();
	if(sptr->colors.size() < numverts)
	{
			ErrorInMod("ElectrostMod::ColorMolSurfElPot()","Colors array size is smaller that number of verticies");
			sptr->ColourUniform(200,200,200);
			return FALSE;
	}

	
	double ecoupl_high_val = -4.0;
	double ecoupl_low_val  = -16.0;

	double fspan = ecoupl_high_val - ecoupl_low_val;
	
	if( fspan < 1E-10) 
	{
			ErrorInMod("ETCouplMod::ColorMolSurfETCoupl()",
				       "The difference between high and low boundary potential values is too small");
			return FALSE;
	}

	int i;
	for(i = 1; i <= numverts; i++ )
	{
		double x = sptr->verts(1,i);
		double y = sptr->verts(2,i);
		double z = sptr->verts(3,i);

		AtomIteratorMolSet aitr(pmset);
		HaAtom* aptr;
		
		double dist_min2 = 1000.0;

		double fval = 0.0;

		for( aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			double dist2 = (x -  aptr->GetX()) * (x -  aptr->GetX());
			dist2 += (y -  aptr->GetY()) * (y -  aptr->GetY());
			dist2 += (z -  aptr->GetZ()) * (z -  aptr->GetZ());

			if( dist2 < dist_min2) 
			{
				dist_min2 = dist2;
				fval = aptr->tempf;
			}
		}

		if( fval > ecoupl_high_val) fval = ecoupl_high_val;
		if( fval < ecoupl_low_val ) fval = ecoupl_low_val;

		double didx_col = 1.0 + ((double) (ncol-1) ) * ( (ecoupl_high_val - fval)/fspan);

		int idx_col = (int) didx_col;
		if( (didx_col - (double)idx_col ) > 0.5 )
			idx_col++;

		if( idx_col < 1) idx_col = 1;
		if( idx_col > ncol ) idx_col = ncol;

		col = col_num(idx_col);
		sptr->colors(i) = col;
	}

	return true;
}

bool ETCouplMod::DuttonModelCalc()
{
	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if( !pdon || !pacc || pdon->empty() || pacc->empty() )
	{
		ErrorInMod("ETCouplMod::DuttonModelCalc()", "Donor or acceptor is not set");
		return false;
	}
	
	double x_mid = 0.0; // point in the middle between donor an acceptor
	double y_mid = 0.0;
	double z_mid = 0.0;

	HaAtom* aptr; 
	
	int nn= 0;
	AtomIteratorAtomGroup aitr_don(pdon);
	for(aptr = aitr_don.GetFirstAtom(); aptr; aptr= aitr_don.GetNextAtom())
	{
		x_mid+= aptr->GetX();
		y_mid+= aptr->GetY();
		z_mid+= aptr->GetZ();
		nn++;
	}

	AtomIteratorAtomGroup aitr_acc(pacc);
	for(aptr = aitr_acc.GetFirstAtom(); aptr; aptr= aitr_acc.GetNextAtom())
	{
		x_mid+= aptr->GetX();
		y_mid+= aptr->GetY();
		z_mid+= aptr->GetZ();
		nn++;
	}

	x_mid /= nn;
	y_mid /= nn;
	z_mid /= nn;

	double max_rad2 = 0.0; // square of the maximum distance from the middle point 
	                       // to atoms of the donor and acceptor

	double d2; 
	for(aptr = aitr_don.GetFirstAtom(); aptr; aptr= aitr_don.GetNextAtom())
	{
		d2 = 0.0;
		d2 += (aptr->GetX() - x_mid)*(aptr->GetX() - x_mid);
		d2 += (aptr->GetY() - y_mid)*(aptr->GetY() - y_mid);
		d2 += (aptr->GetZ() - z_mid)*(aptr->GetZ() - z_mid);
		
		if( d2 > max_rad2) max_rad2 = d2;
	}
	
	for(aptr = aitr_acc.GetFirstAtom(); aptr; aptr= aitr_acc.GetNextAtom())
	{
		d2 = 0.0;
		d2 += (aptr->GetX() - x_mid)*(aptr->GetX() - x_mid);
		d2 += (aptr->GetY() - y_mid)*(aptr->GetY() - y_mid);
		d2 += (aptr->GetZ() - z_mid)*(aptr->GetZ() - z_mid);
		
		if( d2 > max_rad2) max_rad2 = d2;
	}

// increase radius for choosing active atoms by 2.0 Ang
	max_rad2= sqrt(max_rad2);
	max_rad2 = (max_rad2 + 4.0)*(max_rad2 + 4.0);

	AtomGroup active_atoms; // Non-hydrogen atoms within a sphere with a center 
	                      // in the middle point between the donor and acceptor 
	                      // and enclosing donor and acceptor atoms

	AtomIteratorMolSet aitr_mset(pmset);
	for(aptr= aitr_mset.GetFirstAtom(); aptr; aptr= aitr_mset.GetNextAtom())
	{
		if(aptr->IsHydrogen())
			continue;
		d2 = 0.0;
		d2+= ( aptr->GetX() - x_mid)* ( aptr->GetX() - x_mid);
		d2+= ( aptr->GetY() - y_mid)* ( aptr->GetY() - y_mid);
		d2+= ( aptr->GetZ() - z_mid)* ( aptr->GetZ() - z_mid);

		if( d2 < max_rad2 )
			active_atoms.InsertAtom(aptr);
	}

// Probing straight lines between atoms of the donor and acceptor

	double step = 0.2;

	int npt_tot = 0;
	int npt_in_atoms = 0;

	HaAtom* aptr_d;
	HaAtom* aptr_a;

	for( aptr_d = aitr_don.GetFirstAtom(); aptr_d; aptr_d = aitr_don.GetNextAtom())
	{
		for( aptr_a = aitr_acc.GetFirstAtom(); aptr_a; aptr_a = aitr_acc.GetNextAtom())
		{
			double dist = Vec3D::CalcDistance( aptr_d, aptr_a, BOHR_U);
			int np = (int)(dist/step);

			if( np == 0 )
			{
				PrintLog(" Error in ETCouplMod::DuttonModelCalc() \n%s\n",
					     " Distance between donor and acceptor atoms is zero ");
				continue;
			}

			double stepx = (aptr_a->GetX() - aptr_d->GetX())/np;
			double stepy = (aptr_a->GetY() - aptr_d->GetY())/np;
			double stepz = (aptr_a->GetZ() - aptr_d->GetZ())/np;
			
			double don_rad, acc_rad;

			don_rad = HaAtom::ElemDuttonRadius( aptr_d->GetElemNo());
			acc_rad = HaAtom::ElemDuttonRadius( aptr_a->GetElemNo());

			int ibeg = (int)(don_rad/step);
			int iend = np - (int)( acc_rad/step );

			for( int i = ibeg; i < iend ; i++) 
			{
				double px = aptr_d->GetX() + stepx* i;
				double py = aptr_d->GetY() + stepy* i;
				double pz = aptr_d->GetZ() + stepz* i;
				
				bool inside_atom = false; 
				AtomIteratorAtomGroup aitr(&active_atoms);
				for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
				{
					d2 = 0.0;
					d2+= ( aptr->GetX() - px )* ( aptr->GetX() - px );
					d2+= ( aptr->GetY() - py )* ( aptr->GetY() - py );
					d2+= ( aptr->GetZ() - pz )* ( aptr->GetZ() - pz );
					
					double at_rad = HaAtom::ElemVDWRadius( aptr->GetElemNo(), true);
					if( d2 < at_rad*at_rad )
					{
						inside_atom = true;
						break;
					}
				}
				if( inside_atom )
					npt_in_atoms++;

				npt_tot++;	
			}
			
		}
	}
//	double rho; // density ratio 
	
	if( npt_tot == 0) 
	{
		PrintLog("ETCouplMod::DuttonModelCalc(): No space between donor and acceptor \n"); 
		rho = 1.0;
	}
	else
	{
		rho = ((double) npt_in_atoms) /( (double) npt_tot);
	}
	
	beta = 0.9* rho + (1.0 - rho)* 2.8;
	double edge_dist = calc_edge_dist();
	dim_less_coupling = exp(-0.5*beta*edge_dist);
	
	max_rate = 10E13 * exp(-beta*edge_dist);
	
	PrintLog(" Density of the atoms between donor and acceptor: rho = %12.6f \n", rho);
	PrintLog(" Average decay exponent beta = %12.6f \n", beta);
	PrintLog(" Dimensionless coupling H_DA = %12.6e \n", dim_less_coupling);
	PrintLog(" Predicted maximat ET rate k_max = %12.6e \n",max_rate); 

	return true;
}


bool ETCouplMod::CalcGFDonAccOrb()
{
	int nb=ptr_qc_mod->GetNBfunc();
	int nab = ptr_qc_mod->ActBas->GetNBfunc();

	HaMat_double gfull;
	HaMat_double gactive;

    ptr_qc_mod->CalcEPfromMO(gfull,tun_ene);

//	HaMat_double& ss = ptr_qc_mod->GetOvlpMat();

	int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();

	if(ndon == 0 || nacc == 0)
	{
		return false;
	}

	if(ptr_qc_mod->ActBas == &(ptr_qc_mod->AtBasis))
	{
		gactive = gfull;
	}
	else
	{
		LinCombOrb3D::Eval1eOp((LinCombOrb3D*)donor_orbs.bas,(LinCombOrb3D*)acc_orbs.bas,gfull,gactive);	
	}

	LinCombOrb3D don_acc_orbs;

	don_acc_orbs.AddOrbs(donor_orbs);
	don_acc_orbs.AddOrbs(acc_orbs);

	LinCombOrb3D::Eval1eOp(&don_acc_orbs,&don_acc_orbs,gactive,don_acc_gf);

	PrintLog(" Calculated GF matrix \n"); 

	int i,j;
	PrintLog(" labels of Target orbital "); 
	for( i = 0; i < ndon; i++)
	{
		PrintLog("   %s   ", (donor_orbs.GetLabel(i)).c_str() );
	}
	PrintLog("\n");
	
	for( j=0; j < nacc ; j++)
	{
//		const HaAtom* aptr2= acc_orbs[j]->GetFstAtHost();
		PrintLog(" %s  ",(acc_orbs.GetLabel(j)).c_str());
		for(i= 0; i < ndon; i++)
		{
//			const HaAtom* aptr1= donor_orbs[i]->GetFstAtHost();
//			double dist=HaAtom::CalcDistance(aptr1,aptr2,ANGSTROM_U);
//			PrintLog("  %16.9f  %16.9f ", dist, fmat(i+1,j+1));
			PrintLog("  %16.9f ", don_acc_gf.GetVal_idx0(i,ndon+j));
		}
		PrintLog("\n");
	}

	return true;
}


bool ETCouplMod::SetDAdipoleMat()
{
	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if(!pdon || pdon->empty() )
	{
		PrintLog(" Error in ETCouplMod::SetDAdipoleMat() \n"); 
		PrintLog(" donor is empty \n"); 
		return false;
	}
	if(!pacc || pacc->empty() )
	{
		PrintLog(" Error in ETCouplMod::SetDAdipoleMat() \n"); 
		PrintLog(" acceptor is empty \n");
		return false;		
	}
	
	HaVec_double vd(3,0.0);
	HaVec_double va(3,0.0);
	HaVec_double vda(3,0.0);
	std::vector<HaAtom*>::iterator iatr;
	int i,na=0;

	for(iatr=pacc->begin(); iatr != pacc->end(); iatr++)
	{
		na++;
		va(1)+= (*iatr)->GetX_Bohr();
		va(2)+= (*iatr)->GetY_Bohr();
		va(3)+= (*iatr)->GetZ_Bohr();
	}

	for(i=1; i <=3; i++)
		va(i)=va(i)/na;

	na=0;
	for(iatr= pdon->begin(); iatr != pdon->end(); iatr++)
	{
		na++;
		vd(1)+= (*iatr)->GetX_Bohr();
		vd(2)+= (*iatr)->GetY_Bohr();
		vd(3)+= (*iatr)->GetZ_Bohr();
	}


	for(i=1; i <=3; i++)
		vd(i)=vd(i)/na;

	HaVec_double cntr(3);
	for(i=1; i<=3; i++)
		cntr(i)= (va(i) + vd(i))/2.0; // Point in the center of the space
	                                    // between donor and acceptor	
	double rr=0.0;
	
	for(i=1; i <= 3; i++)
	{
		vda(i)=va(i)-vd(i);
		rr+= vda(i)*vda(i);
	}

	PrintLog(" Donor center %8.3f %8.3f %8.3f ", vd[0]*BOHR_TO_ANG, vd[1]*BOHR_TO_ANG, vd[2]*BOHR_TO_ANG );
	PrintLog(" Acceptor center %8.3f %8.3f %8.3f ", va[0]*BOHR_TO_ANG, va[1]*BOHR_TO_ANG, va[2]*BOHR_TO_ANG );
	PrintLog(" Molecule center %8.3f %8.3f %8.3f ", cntr[0]*BOHR_TO_ANG, cntr[1]*BOHR_TO_ANG, cntr[2]*BOHR_TO_ANG );

	if( rr <= 1.0)
	{
		std::cerr << " Error in  ETCouplMod::SetDAdipoleMat() " << std::endl;   
		std::cerr << " donor and acceptor coincide " << std::endl;
		std::cerr << " Donor-acceptor field dipole matrix is not set " << std::endl;
		return false;
	}

	rr=sqrt(rr);
	for(i=1; i<=3; i++)
		vda(i)=vda(i)/rr;

	int nab= ptr_qc_mod->ActBas->GetNBfunc();

	if(DA_dipole.num_rows() != nab || DA_dipole.num_cols() != nab) DA_dipole.newsize(nab,nab);

	DA_dipole=0.0;
    HaMat_double tmp_mat;
    HaMat_double ss_mat;
    HaMat_double scr;    
        
	HaOperR rop;

	if(ssl.num_rows() != nab)
	{
		HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas, ptr_qc_mod->ActBas, ss_mat);
	}
	else
	{
		ss_mat = ssl;
	}

    HaMat_doubleArr rmats;
	rop.FillMat(ptr_qc_mod->ActBas, rmats);
	
	//for( i = 0; i < nab; i++ )
	//{
	//	Vec3D* phost = ptr_qc_mod->ActBas->GetHostPt(i);
	//	HaAtom* pat_host = (HaAtom*) phost; 
	//	std::string lbl = ptr_qc_mod->ActBas->GetLabel(i);
	//	std::vector<double> rdiag(3);

	//	rdiag[0] = rmats[0].r0(i,i);
	//	rdiag[1] = rmats[1].r0(i,i);
	//	rdiag[2] = rmats[2].r0(i,i);

	//	PrintLog(" %d  at= %s   at_crd= %8.3f %8.3f %8.3f   r_diag= %8.3f %8.3f %8.3f ",
	//		       i, pat_host->GetRef().c_str(), pat_host->GetX_Ang(), pat_host->GetY_Ang(), pat_host->GetZ_Ang(),
	//			   rdiag[0]*BOHR_TO_ANG, rdiag[1]*BOHR_TO_ANG,rdiag[2]*BOHR_TO_ANG);
	//}
	
			
	for(i= 0; i < 3; i++)
	{
		mat_scale(tmp_mat,ss_mat,cntr[i]);
		scr = rmats[i];
		mat_diff(scr,scr,tmp_mat);
		mat_scale(scr, scr,vda[i]);
		mat_add(DA_dipole, DA_dipole, scr);
	}

	return true;
	
}

bool ETCouplMod::SetDAfield(double field)
{
	if(ptr_qc_mod == NULL) return false;
	int nab = ptr_qc_mod->GetNActiveOrb();
	bool result;
	if( DA_dipole.num_rows() == 0 && DA_dipole.num_cols() == 0 && (AbsFun(field) < 1e-12))
	{
		DA_field=field;
		return true;
	}
	if(nab != DA_dipole.num_rows() || nab != DA_dipole.num_cols()) 
	{
		result=this->SetDAdipoleMat();
		if(!result) return false;
	}
	DA_field=field;

	mat_scale(extern_field, DA_dipole,DA_field); 

	return true;

}

class fun_DA_split: public UnivariateFunctor
// Axxiliary function to realize one-parameter function to 
// calculate splitting between donor/acceptor orbitals 
{
public:
	fun_DA_split(ETCouplMod& new_et_mod): et_mod(new_et_mod)
	{
		idon_min = 0;
		iacc_min = 0;
		counter = 0;
		with_sign = 0;
	}
    
	~fun_DA_split() {}
	virtual double operator() (const double x) 
	{  
		bool result;
		result= et_mod.SetDAfield(x);
		if(!result)
		{
			PrintLog(" Error in fun_DA_split() \n "); 
			PrintLog(" Error to apply DA electric field \n");
			return 0.0;
		}
		result=et_mod.DiagHeff();
        HaVec_int idx_da = et_mod.FindDonAccEigVecs(idon_min,iacc_min);

		if(idx_da[0] == -1)
		{
			PrintLog(" Error in fun_DA_split() \n");
			PrintLog(" Error to find donor/acceptor orbitals \n "); 
			return 0.0;
		}
		counter++;
        double split = et_mod.enel[idx_da[0]]- et_mod.enel[idx_da[1]]; 
		if( idx_da[3] > 0 ) split = -split;
		PrintLog(" min step= %4d  field= %20.12f ene_split= %12.6e \n",
                           counter,x,split); 
		if(!with_sign) split = fabs(split);
		return split;
	}

	int idon_min;
	int iacc_min;
	int with_sign;
protected:
	int counter;
	ETCouplMod& et_mod;

};

HaVec_int ETCouplMod::FindDonAccEigVecs(int idx_don, int idx_acc)
{
	HaVec_int idx_max(4);
	idx_max = -1;

	int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();
	
    if( idx_don < 0 || ndon <= idx_don)
	{
		PrintLog(" Error in: ETCouplMod::FindDonAccEigVecs() \n ");
		PrintLog(" Invalid donor orbital index \n");
		return idx_max;
	}

    if( idx_acc < 0 || nacc <= idx_acc)
	{
		PrintLog(" Error in: ETCouplMod::FindDonAccEigVecs() \n ");
		PrintLog(" Invalid acceptor orbital index \n");
		return idx_max;
	}

	HaMat_double ovlp_mat_don;
	HaMat_double ovlp_mat_acc;

	HaBasisSet::CalcOvlpMat(&donor_orbs, &eigv, ovlp_mat_don);
	HaBasisSet::CalcOvlpMat(&acc_orbs, &eigv, ovlp_mat_acc);

//	LinCombOrb3D::Eval1eOp(&donor_orbs, &eigv, ssl, ovlp_mat_don);
//	LinCombOrb3D::Eval1eOp(&acc_orbs, &eigv, ssl, ovlp_mat_acc);

	int j;
	int neig = eigv.GetNOrbs();

	double ov_max_1 = 0.0;
	double ov_max_2 = 0.0;
	
	for(j = 0; j < neig; j++)
	{
		double ov_don = ovlp_mat_don.r0(idx_don,j);
		double ov_acc = ovlp_mat_acc.r0(idx_acc,j);

		int isign = 1;
		if( ov_don*ov_acc < 0) isign = -1;
		
		double ov = sqrt(ov_don*ov_don + ov_acc*ov_acc);
		
		if(ov > ov_max_1) 
		{
			ov_max_2 = ov_max_1; 
			idx_max[1] = idx_max[0];
			idx_max[3] = idx_max[2]; // relative signs of donor and acceptor overlaps 
			idx_max[0] = j;
            idx_max[2] = isign;
			ov_max_1 = ov;
		}
		if( ov > ov_max_2 && idx_max[0] != j)
		{
            ov_max_2 = ov;
			idx_max[1] = j;
			idx_max[3] = isign;
		}
	}
	PrintLog(" Max Ovlp of donor orb %d and acceptor %d \n is with eigen vectors %d and %d = %12.6f and %12.6f\n",
		idx_don+1, idx_acc, idx_max[0]+1,idx_max[1]+1, ov_max_1,ov_max_2 );
	PrintLog(" relative signs %d %d \n", idx_max[2],idx_max[3]);

	return idx_max;
}

bool ETCouplMod::RotateRedoxOrb( const HaMat_double& rot_mat )
{
	return false;
}

bool ETCouplMod::FindRedoxOrbsOvlpEigVecs( StrIntMap& lbl_idx_map, StrDoubleMap& ovlp_val_map, REDOX_ORB_TYPE redox_orb_type  )
{
	HaMat_double ovlp_mat;

	LinCombOrb3D* p_rdx_orb = &donor_orbs;
	if( redox_orb_type == REDOX_ORB_ACCEPTOR ) p_rdx_orb = &acc_orbs;

	HaVec_int* p_ieig_rdx = &ieig_don;
	if( redox_orb_type == REDOX_ORB_ACCEPTOR ) p_ieig_rdx = &ieig_acc;

	int nrdx_orb = p_rdx_orb->GetNOrbs();

	lbl_idx_map.clear();
	ovlp_val_map.clear();

	PrintLog(" ETCouplMod::FindRedoxOrbsOvlpEigVecs() pt 1 \n");

	try
	{
		if( ptr_qc_mod->ActBas == NULL ) throw std::runtime_error(" Active Basis is not set ");
		
		int nb = ptr_qc_mod->ActBas->GetNBfunc();
		HaMat_double ssl_loc = GetActBasOvlpMat();
		
		LinCombOrb3D::Eval1eOp(p_rdx_orb, &eigv, ssl_loc, ovlp_mat);

		int i,j;
		int neig = eigv.GetNOrbs();

		p_ieig_rdx->newsize(nrdx_orb);

		for( i = 0; i < nrdx_orb; i++)
		{
			double ov_max = 0.0;
			std::string lbl = p_rdx_orb->GetLabel(i);
			int imax = -1;
			for(j = 0; j < neig; j++)
			{
				double ov_val = fabs(ovlp_mat.GetVal_idx0(i,j));
				if(ov_val > ov_max) 
				{
					ov_max = ov_val;
					imax = j;
				}
			}
			lbl_idx_map[lbl] = imax+1;
			ovlp_val_map[lbl] = ov_max;
//			PrintLog(" Max Ovlp of donor orb %d lbl = %s is with eigen vector %d  = %12.6f \n", i+1,lbl.c_str(), imax+1, ov_max);
			(*p_ieig_rdx)[i] = imax+1; 
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ETCouplMod::FindRedoxOrbsOvlpEigVecs() \n");
		PrintLog("%s\n",ex.what());
		lbl_idx_map.clear();
		ovlp_val_map.clear();
		return false;
	}
	return true;
}


bool ETCouplMod::FindRedoxOrbSpaceOvlpEigVecs( StrIntMap& lbl_idx_map, StrDoubleMap& ovlp_val_map, HaVec_double& eigv_space_max_ovlp_val, REDOX_ORB_TYPE redox_orb_type ) 
{
	HaMat_double ovlp_mat;
	HaMat_double ovlp_rdx_orb;
	HaMat_double ovlp_msqrt_rdx_orb; // S^(-1/2) for donor orb

	LinCombOrb3D* p_rdx_orb = &donor_orbs;
	if( redox_orb_type == REDOX_ORB_ACCEPTOR ) p_rdx_orb = &acc_orbs;

	HaVec_int* p_ieig_rdx = &ieig_don;
	if( redox_orb_type == REDOX_ORB_ACCEPTOR ) p_ieig_rdx = &ieig_acc;

	int nrdx_orb = p_rdx_orb->GetNOrbs();

	PrintLog(" ETCouplMod::FindRedoxOrbSpaceOvlpEigVecs() pt 1 \n");

	try
	{
		if( ptr_qc_mod->ActBas == NULL ) throw std::runtime_error(" Active Basis is not set ");
		if( nrdx_orb == 0 ) throw std::runtime_error(" Redox orbitals are not set ");
		 
		int nb = ptr_qc_mod->ActBas->GetNBfunc();
		HaMat_double ssl_loc = GetActBasOvlpMat();
		
		LinCombOrb3D::Eval1eOp(p_rdx_orb, p_rdx_orb, ssl_loc, ovlp_rdx_orb);
		ovlp_msqrt_rdx_orb = ovlp_rdx_orb;
		int ires = ovlp_msqrt_rdx_orb.SqRoot(-1);
		if( ires != TRUE ) throw std::runtime_error(" Failed to find Square Root of redox orbitals overlap matrix ");

		LinCombOrb3D orthog_redox_orbs; //!< Orthogonalized redox orbitals
		orthog_redox_orbs.bas  = p_rdx_orb->bas;
		matmult(orthog_redox_orbs.coef, p_rdx_orb->coef, ovlp_msqrt_rdx_orb); 

		LinCombOrb3D::Eval1eOp(&orthog_redox_orbs, &eigv, ssl_loc, ovlp_mat);

		int i,j;
		int neig = eigv.GetNOrbs();

		if( neig == 0 ) throw std::runtime_error(" Eigen Vectors are not set ");
		if( neig <  nrdx_orb ) throw std::runtime_error(" The number of eigenvectors is smaller than the number of redox orvitals  ");

		p_ieig_rdx->clear();  // indexes of eigen vectors with maximal overlap with redox orbitals space   
		eigv_space_max_ovlp_val.clear();  // values of maximal overlap of eigen vectors with redox orbitals space 

		HaVec_double redox_space_eigv_ovlp_all;
		redox_space_eigv_ovlp_all.resize(neig);

		std::vector<int> max_ovlp_eigv_idx_1;

		for(j = 0; j < neig; j++)
		{
			double ov_rdx_space = 0.0;
			for( i = 0; i < nrdx_orb; i++ )
			{
				ov_rdx_space += ovlp_mat.r0(i,j)*ovlp_mat.r0(i,j);
			}
			ov_rdx_space = sqrt( ov_rdx_space );
			redox_space_eigv_ovlp_all[j] = ov_rdx_space;

			if( max_ovlp_eigv_idx_1.size() < nrdx_orb )
			{
				max_ovlp_eigv_idx_1.push_back( j );
				eigv_space_max_ovlp_val.push_back( ov_rdx_space );
				continue;
			}
			
			double ovlp_min = 1000.0;  
			int idx_min = -1;
			for( i = 0; i < nrdx_orb; i++ ) 
			{
				if( eigv_space_max_ovlp_val[i] < ovlp_min )
				{
					ovlp_min = eigv_space_max_ovlp_val[i];
					idx_min = i;
					continue;
				}
			}
			if( ovlp_min < ov_rdx_space )
			{
				eigv_space_max_ovlp_val[idx_min] = ov_rdx_space;
				max_ovlp_eigv_idx_1[idx_min] = j;
//				(*p_ieig_rdx)[idx_min] = j + 1;
			}
		}

		if( max_ovlp_eigv_idx_1.size() != nrdx_orb ) throw std::runtime_error(" wrong number of overlap values ");
		
		std::sort( max_ovlp_eigv_idx_1.begin(), max_ovlp_eigv_idx_1.end() );

		p_ieig_rdx->resize( nrdx_orb );
		eigv_space_max_ovlp_val.resize( nrdx_orb );

		for(i = 0; i < nrdx_orb; i++ )
		{
			int idx = max_ovlp_eigv_idx_1[i];
			(*p_ieig_rdx)[i] = idx + 1;
			eigv_space_max_ovlp_val[i] = redox_space_eigv_ovlp_all[idx];
		}

		lbl_idx_map.clear();
		ovlp_val_map.clear();

		for(i = 0; i < nrdx_orb; i++)
		{
			std::string lbl = p_rdx_orb->GetLabel(i);
			lbl_idx_map[lbl] = 0;
			ovlp_val_map[lbl] = 0.0;
		}

		for(j = 0; j < nrdx_orb; j++)
		{
			int idx_eigv = (*p_ieig_rdx)[j];
			for( i = 0; i < nrdx_orb; i++ )
			{
				std::string lbl = p_rdx_orb->GetLabel(i);
				double ovlp_old = ovlp_val_map[lbl];
				double ovlp = ovlp_mat.r0(i,idx_eigv-1); 
				if( ovlp > ovlp_old )
				{
					lbl_idx_map[lbl]  = idx_eigv;
					ovlp_val_map[lbl] = ovlp;
				}
			}
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ETCouplMod::FindRedoxOrbSpaceOvlpEigVecs() \n");
		PrintLog("%s\n",ex.what());
		p_ieig_rdx->clear();
		eigv_space_max_ovlp_val.clear();
		p_ieig_rdx->clear();
		return false;
	}
	return true;
}


HaMat_double& ETCouplMod::GetActBasOvlpMat()
{
    if(ptr_qc_mod->ActBas == NULL) 
	{
		ssl.newsize(0,0);
		return ssl;
	}
	int nab = ptr_qc_mod->ActBas->GetNBfunc();
	if( ssl.num_cols() == nab && ssl.num_rows() == nab)
	{
		return ssl;
	}
	ssl.newsize(nab,nab);
	
	if( ptr_qc_mod->ActBas == &(ptr_qc_mod->AtBasis))
	{
		if( ptr_qc_mod->ovlp_mat.num_cols() == nab) 
		{
			ssl = ptr_qc_mod->ovlp_mat;
			return ssl;
		}
	}
	HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,ssl);
	return ssl;
}

bool ETCouplMod::CalcHDAEneSplit()
{
	if(ptr_qc_mod->ActBas == NULL)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAEneSplit() \n");
		PrintLog(" Active Basis set is not set \n");
		return false;
	}
	int nab = ptr_qc_mod->ActBas->GetNBfunc();
	if( nab == 0)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAEneSplit() \n");
		PrintLog(" Active Basis set is empty \n");
		return false;		
	}

	if( heff_mat.num_rows() != nab || heff_mat.num_cols() != nab)
	{
		RecalcHeff();
	}

	if( heff_mat.num_rows() != nab || heff_mat.num_cols() != nab)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAEneSplit() \n");
		PrintLog(" Unable to recompute one electron hamiltonian \n");
		return false;		
	}

	int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();

	int id;
	int ia;
	
	da_coupl_val.newsize(ndon,nacc);
	el_field_min.newsize(ndon,nacc);

	for(id = 0; id < ndon; id++)
	{
		for( ia = 0; ia < nacc; ia++)
		{
			double delt=0.00000001;
			double min_field=0.0;
			double min_split=0.0;

			fun_DA_split fsplit(*this);
			fsplit.idon_min = id;
			fsplit.iacc_min = ia;
	
			double x0,x1,x2,f1,f2;

			x0= 0.0;

			f1=fsplit(x0);
			f2=fsplit(x0+delt);

			if(fabs(f1-f2) < DBL_EPSILON)
			{
				PrintLog(" Error in ETCouplMod::CalcHDAEneSplit() \n"); 
				PrintLog(" f1 - f2 == 0 \n"); 
				return false;
			}

			if(f2 > f1)
			{
				x1= x0 - 2*delt*f1/(f2-f1);
				x2= x0;
			}
			else
			{
				x1=x0;
				x2=x0 + 2*delt*f1/(f1-f2);
			}
	
			min_field = fminbr(x1,x2,fsplit);
            fsplit.with_sign = TRUE;
			min_split = fsplit(min_field);

			el_field_min.r0(id,ia) = min_field;
			da_coupl_val.r0(id,ia) = min_split/2.0;
		}
	}
		
	PrintLog(" Couplings between donor and acceptor orbitals: \n");
	for(id=0; id < ndon; id++)
	{
		for(ia=0; ia < nacc; ia++)
		{
			std::string lbld = donor_orbs.GetLabel(id);
			std::string lbla = acc_orbs.GetLabel(ia);
			PrintLog(" %s %s  %12.6e  at field %12.9f \n", 
				lbld.c_str(), lbla.c_str(), da_coupl_val.GetVal_idx0(id,ia), el_field_min.GetVal_idx0(id,ia));
		}
	}

	return true;
}

int ETCouplMod::AddRedoxOrbFromEigVec(const HaVec_int& mo_idx)
{
	MolSet* pmset = GetMolSet();

	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
    AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	bool has_donor = false;
	bool has_acc = false;

	if( pdon && !pdon->empty() ) has_donor = true;
	if( pacc && !pacc->empty() ) has_acc = true;

	try
	{  
		if( !has_donor && !has_acc ) throw std::runtime_error(" donor and acceptor groups are not set ");
		int i,j;
		int idx;
		int nmo_idx = mo_idx.size();

		PrintLog(" eigv.GetNOrbs() = %d \n", eigv.GetNOrbs() );

		for( i = 0; i < mo_idx.size(); i++)
		{
			idx = mo_idx[i];
			if( idx < 1 || idx > eigv.GetNOrbs()  ) throw std::runtime_error("Invalid index of the donor or acceptor orbital " + harlem::ToString(idx));
		}

		int nab = ptr_qc_mod->ActBas->GetNBfunc();

		if( nab != eigv.bas->GetNBfunc())
		{
			throw std::runtime_error(" The number of rows in Eigen Vector matrix " + harlem::ToString(eigv.bas->GetNBfunc()) + 
			       "\n Doesn't correspond to the number of active orbitals " + harlem::ToString(nab) ); 
		}

		LinCombOrb3D tmp_rdx_orb;
		tmp_rdx_orb.bas = ptr_qc_mod->ActBas;
		tmp_rdx_orb.coef.newsize(nab,nmo_idx*2);
		tmp_rdx_orb.coef = 0.0;

		HaVec_int don_orb_idx(nab);
		HaVec_int acc_orb_idx(nab);

		int nb_don = 0;
		int nb_acc = 0;

		HaMat_double eigv_f(nab,nmo_idx); // expansion coefficients for trial mos

		for(i=0; i < nab; i++) 
		{
			for(j = 0 ; j < nmo_idx; j++)
			{
				idx = mo_idx[j];
				eigv_f.r0(i,j) = eigv.coef.r0(i, idx-1);
			}
			Vec3D* hpt = ptr_qc_mod->ActBas->GetHostPt(i);
			const HaAtom* pat_host= (const HaAtom*) hpt;
			if( has_donor && pdon->HasAtom(pat_host))
			{		
				don_orb_idx[nb_don] = i;
				nb_don++;
			}

			if( has_acc && pacc->HasAtom(pat_host))
			{		
				acc_orb_idx[nb_acc] = i;
				nb_acc++;
			}
		}

		HaMat_double don_ss_mat(nb_don,nb_don); // overlap matrix of atomic orbitals residing on donor 
		HaMat_double acc_ss_mat(nb_acc,nb_acc); // overlap matrix of atomic orbitals residing on acceptor

		HaMat_double& ss = GetActBasOvlpMat();

	//	PrintLog(" ETCouplMod::AddRedoxOrbFromEigVec() ss.r0(0,100) = %12.6f \n", ss.r0(0,100) );

		for(i=0; i < nb_don; i++)
		{
			for(j = 0; j < nb_don; j++)
			{
				don_ss_mat.r0(i,j) = ss.r0( don_orb_idx[i],don_orb_idx[j]);
			}
		}

		for(i=0; i < nb_acc; i++)
		{
			for(j = 0; j < nb_acc; j++)
			{
				acc_ss_mat.r0(i,j) = ss.r0( acc_orb_idx[i],acc_orb_idx[j]);
			}
		}

		HaMat_double ss_don_tr(nb_don,nmo_idx);
		HaMat_double ss_acc_tr(nb_acc,nmo_idx);
		HaMat_double tmp;

		matmult(tmp,ss,eigv_f);  // S C

		for(j = 0; j < nmo_idx; j++)
		{
			for(i = 0; i < nb_don; i++)
			{
				ss_don_tr.r0(i,j) = tmp.r0(don_orb_idx[i],j); // S C_tr
			}
		}

		for( j = 0; j < nmo_idx; j++)
		{
			for(i = 0; i < nb_acc; i++)
			{
				ss_acc_tr.r0(i,j) = tmp.r0(acc_orb_idx[i],j);
			}
		}

		if( nb_don > 0 ) HaMat_double::solv_lin_syst_1(don_ss_mat,ss_don_tr); //  C' = S_DD(-1) S C_tr
		if( nb_acc > 0 ) HaMat_double::solv_lin_syst_1(acc_ss_mat,ss_acc_tr); //  C' = S_AA(-1) S C_tr

		for( j = 0; j < nmo_idx; j++)
		{
			for(i = 0; i < nb_don; i++)
			{
				tmp_rdx_orb.coef.r0(don_orb_idx[i],j) =  ss_don_tr.r0(i,j);
			}
		}

		for( j = 0; j < nmo_idx; j++)
		{
			for(i = 0; i < nb_acc; i++)
			{
				tmp_rdx_orb.coef.r0(acc_orb_idx[i],j+nmo_idx) = ss_acc_tr.r0(i,j);
			}
		}

		HaMat_double self_ovlp;

		HaBasisSet::CalcOvlpMat(&tmp_rdx_orb,&tmp_rdx_orb,self_ovlp);

		PrintLog(" Norms of truncated Eigen Functions residing on the Redox Centers: \n");  
		for( idx = 0; idx < nmo_idx; idx++)
		{
			double snorm_d = sqrt(self_ovlp.GetVal_idx0(idx,idx));
			double snorm_a = sqrt(self_ovlp.GetVal_idx0(idx+nmo_idx,idx+nmo_idx));

			int use_d = FALSE;
			int use_a = FALSE;

			if( snorm_d > 0.6) use_d = TRUE;
			if( snorm_a > 0.6) use_a = TRUE;

			PrintLog("  %12.6f %12.6f \n", 
				snorm_d, snorm_a);

			if(snorm_d < 0.6 && snorm_a < 0.6) 
			{
				PrintLog(" Possible Error setting Donor or Acceptor  orbital: \n");
				PrintLog(" Truncated one-el eigen function %d has small norms on donor %12.6f and acceptor %12.6f\n",
					mo_idx[idx], snorm_d, snorm_a);
				continue;
			}

			for(i=0; i< nab; i++)
			{
				if(use_d) tmp_rdx_orb.coef.SetVal_idx0(i,idx, tmp_rdx_orb.coef.GetVal_idx0(i,idx)/snorm_d);
				if(use_a) tmp_rdx_orb.coef.SetVal_idx0(i,idx+nmo_idx, tmp_rdx_orb.coef.GetVal_idx0(i,nmo_idx+idx)/snorm_a);
			}

			if(use_d || use_a)
			{
				int k;
				for( k = 0; k < 2; k++)
				{
					LinCombOrb3D* prdx_orb = NULL;
					int ishift = 0;

					if(k==0) 
					{
						if(!use_d) continue;
						prdx_orb = &donor_orbs;
					}
					if(k==1)
					{
						if(!use_a) continue;
						prdx_orb = &acc_orbs;
						ishift = nmo_idx;
					}
					if(prdx_orb->bas != ptr_qc_mod->ActBas)
					{
						prdx_orb->Clear();
						prdx_orb->bas = ptr_qc_mod->ActBas;
					}

					int n_old = prdx_orb->GetNOrbs();
					HaMat_double tmp2(nab, n_old+1);

					for(i = 0; i < n_old; i++)
					{
						for(j = 0; j < nab; j++)
						{
							tmp2.SetVal_idx0(j,i,prdx_orb->coef.GetVal_idx0(j,i));
						}
					}
					for(j = 0; j < nab; j++)
					{
						tmp2.SetVal_idx0(j,n_old,tmp_rdx_orb.coef.GetVal_idx0(j,ishift + idx));
					}
					prdx_orb->coef = tmp2;
					char buf[128];
					sprintf(buf,"%d",mo_idx[idx]);
					std::string s_n = buf;
					boost::trim(s_n);
					std::string lbl = "TRUNC_EIGV_" + s_n;
					prdx_orb->ids.push_back(lbl);
				}
			}
		}
	}
	catch(const std::exception& ex )
	{
		PrintLog("Error in ETCouplMod::AddRedoxOrbFromEigVec() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;	
}

int ETCouplMod::GetRedoxOrbsFromFrag( MolSet* pfrag )
{
	MolSet* pmset = this->GetMolSet();
	
	try
	{
		if( pfrag == NULL ) throw std::runtime_error(" no fragment specified ");
		ETCouplMod* p_etmod_frag = pfrag->GetETCouplMod( false );
		if( p_etmod_frag == NULL ) throw std::runtime_error(" fragment ET module is not specified ");

		AtomGroup* donor = pmset->GetAtomGroupByID("DONOR");
		AtomGroup* acceptor = pmset->GetAtomGroupByID("ACCEPTOR");

		AtomGroup* donor_frag = pfrag->GetAtomGroupByID("DONOR");
		AtomGroup* acceptor_frag = pfrag->GetAtomGroupByID("ACCEPTOR");

		bool set_donor_orb = false;
		bool set_acceptor_orb = false;

		if( donor == NULL && acceptor == NULL ) throw std::runtime_error(" main molecule does not have donor and acceptor groups ");
		if( donor_frag == NULL && acceptor_frag == NULL ) throw std::runtime_error(" fragment does not have donor and acceptor groups ");
		if( donor_frag != NULL ) 
		{
			if( donor == NULL ) throw std::runtime_error(" fragment has a donor group but main molecule doesn't ");
			set_donor_orb = true;
		}
		if( acceptor_frag != NULL ) 
		{
			if( acceptor == NULL ) throw std::runtime_error(" fragment has an acceptor group but main molecule doesn't ");
			set_acceptor_orb = true;
		}

		int irdx;
		for( irdx = 0; irdx < 2; irdx++ )
		{
			if( irdx == 0 && !set_donor_orb ) continue;
			if( irdx == 1 && !set_acceptor_orb ) continue;

			AtomGroup* rdx_grp = donor;
			AtomGroup* rdx_grp_frag = donor_frag;
			if( irdx == 1 )
			{
				rdx_grp = acceptor;
				rdx_grp_frag = acceptor_frag;
			}

			AtomAtomMap frag_atom_map; 
		
			int i;
			int na = rdx_grp->GetNAtoms();
			if( na == 0 )  throw std::runtime_error(" redox group of the main molecule is empty ");
			if( rdx_grp_frag->GetNAtoms() != na )  throw std::runtime_error(" fragment and main molecule redox groups have different numbers of atoms ");

			//for( i = 0; i < na; i++ )
			//{
			//	frag_atom_map[ rdx_grp->at(i) ]      = rdx_grp_frag->at(i); 
			//	frag_atom_map[ rdx_grp_frag->at(i) ] = rdx_grp->at(i); 
			//}

			pmset->BuildFragmentAtomMap(pfrag, frag_atom_map);

			HaMat_double rot_mat;
			HaVec_double trans_vec;
//			HaVec_double crd_frag = rdx_grp_frag->GetCrdArray();
			
			HaVec_double crd_frag( 3*na, 0.0);

			for( i= 0; i < na; i++)
			{
				HaAtom* aptr_f = rdx_grp_frag->at(i);
				if ( frag_atom_map.count( aptr_f ) == 0 )  throw std::runtime_error( " fragment atom " + aptr_f->GetRef() + " is not mapped ");  
			}

			for( i= 0; i < na; i++)
			{
				HaAtom* aptr = rdx_grp->at(i);
				if ( frag_atom_map.count( aptr ) == 0 )  throw std::runtime_error( " redox group atom  " + aptr->GetRef() + " is not mapped ");  
				HaAtom* aptr_f = frag_atom_map[aptr];
				crd_frag.r0(3*i  ) = aptr_f->GetX_Ang();
				crd_frag.r0(3*i+1) = aptr_f->GetY_Ang();
				crd_frag.r0(3*i+2) = aptr_f->GetZ_Ang();
			}

			double eps = harlem::geom::GetSuperimposeMat( crd_frag, *rdx_grp, rot_mat, trans_vec );
			PrintLog(" ETCouplMod::GetRedoxOrbsFromFrag() RMSD for fit of a fragment redox group  eps= %12.6f \n ", eps );
			if( eps < 0.0 ) throw std::runtime_error(" Failed to fit atom coordinates of redox group of the fragment " );


			if( ptr_qc_mod->ActBas == NULL ) throw std::runtime_error( "  Active Basis is not set for main molecule " );
			if( p_etmod_frag->ptr_qc_mod->ActBas == NULL ) throw std::runtime_error( "  Active Basis is not set for the fragment " );

			int nb = ptr_qc_mod->ActBas->GetNBfunc();
			int nb_f = p_etmod_frag->ptr_qc_mod->ActBas->GetNBfunc();

			AtomIntMap low_idx;
			AtomIntMap high_idx;
			AtomIntMap low_idx_f;
			AtomIntMap high_idx_f;

			for( i = 0; i < nb; i++ )
			{
				Vec3D* pt = ptr_qc_mod->ActBas->GetHostPt(i);
				HaAtom* aptr = dynamic_cast<HaAtom*>( pt );
				if( !rdx_grp->HasAtom( aptr ))  continue;
				if( low_idx.count( aptr) == 0 )
				{
					low_idx[ aptr ] = i;
					high_idx[ aptr ] = i;
					continue;
				}
				high_idx[ aptr ] = i;
			}

			for( i = 0; i < nb_f; i++ )
			{
				Vec3D* pt_f =  p_etmod_frag->ptr_qc_mod->ActBas->GetHostPt(i);
				HaAtom* aptr_f = dynamic_cast<HaAtom*>( pt_f );
				if( !rdx_grp_frag->HasAtom( aptr_f ))  continue;
				if( low_idx_f.count( aptr_f ) == 0 )
				{
					low_idx_f[ aptr_f ] = i;
					high_idx_f[ aptr_f ] = i;
					continue;
				}
				high_idx_f[ aptr_f ] = i;
			}

			LinCombOrb3D* redox_orbs      = &this->donor_orbs;
			LinCombOrb3D* redox_orbs_frag = &p_etmod_frag->donor_orbs;
			if( irdx == 1 )
			{
				redox_orbs      = &this->acc_orbs;
				redox_orbs_frag = &p_etmod_frag->acc_orbs;
			}

			int norb = redox_orbs_frag->GetNOrbs();
			if( norb == 0 ) throw std::runtime_error(" No redox orbitals are set on the fragment ");

			HaMat_double cf_mat;
			cf_mat.newsize( nb, norb );
			cf_mat = 0.0;

			PrintLog(" ETCouplMod::GetRedoxOrbsFromFrag() pt 1  nb = %d nb_f = %d \n", nb, nb_f); 

			int j;
			int i_f;
			HaAtom* aptr_fst = NULL;

			for( i_f = 0; i_f < nb_f; i_f++ )
			{
				HaAtom* aptr_f = dynamic_cast<HaAtom*>( p_etmod_frag->ptr_qc_mod->ActBas->GetHostPt(i_f) );
				if( !rdx_grp_frag->HasAtom( aptr_f ))  continue;
				HaAtom* aptr  = frag_atom_map[aptr_f];

				if( aptr_fst == NULL ) aptr_fst = aptr;

				int idx_l = low_idx [ aptr ];
				int idx_h = high_idx[ aptr ];
				int idx_l_f = low_idx_f [ aptr_f ];
				int idx_h_f = high_idx_f[ aptr_f ];

				if( (idx_h - idx_l) != (idx_h_f - idx_l_f) )
				{
					throw std::runtime_error("Number of basis functions is different for mapped atoms " + aptr_f->GetRef() + " and " + aptr->GetRef() );
				}

				i = i_f - idx_l_f + idx_l;
				for( j = 0; j < norb; j++ )
				{
					cf_mat.r0( i, j ) =  redox_orbs_frag->coef.r0( i_f, j );  
				}
			}	

			PrintLog("Set Redox Orbitals from fragment \n");
			redox_orbs->Clear();	
			redox_orbs->bas = ptr_qc_mod->ActBas;
			redox_orbs->coef = cf_mat;
			redox_orbs->ids.resize(norb);
			redox_orbs->at_ptr.resize(norb);
			for( j = 0; j < norb; j++ )
			{
				//				std::string lbl = "DON_ORB_" + harlem::ToString(j);
				//				this->donor_orbs.ids[j] = lbl;
				std::string lbl = redox_orbs_frag->GetLabel(j);
				redox_orbs->ids[j] = lbl;
				redox_orbs->at_ptr[j] = aptr_fst;
			}	

			redox_orbs->TrCoefRot( rot_mat );
		}
		PrintLog(" ETCouplMod::GetRedoxOrbsFromFrag() Number of donor orbitals: %d ", this->donor_orbs.GetNOrbs() );
		PrintLog(" ETCouplMod::GetRedoxOrbsFromFrag() Number of acceptor orbitals: %d ", this->acc_orbs.GetNOrbs() );
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in ETCouplMod::GetRedoxOrbsFromFrag() \n%s\n", ex.what() );
		return TRUE;
	}
	return FALSE;
}

bool ETCouplMod::PrintEigVecCoef(int idx)
{
	HaQCMod* ptr_qc_mod = GetQCMod();
	
	if(ptr_qc_mod == NULL)
	{
		std::cerr << " Error in ETCouplMod::PrintEigVecCoef() " << std::endl;
		std::cerr << " Qchem module is not set for ET Coupling module " << std::endl;
		return false;
	}
	
	if(idx > eigv.GetNOrbs() && idx <= 0)
	{
		std::cerr << " Error in ETCouplMod::PrintEigVecCoef() " << std::endl;
		std::cerr << " Index of the effective orbital " << idx << " Is out of range " << std::endl;
		return false;
	}
	
	int nao= ptr_qc_mod->GetNActiveOrb();

	if(eigv.bas->GetNBfunc() != nao)
	{
		std::cerr << " Error in ETCouplMod::PrintEigVecCoef() " << std::endl;
		std::cerr << " The Number of rows in eigenvector matrix " << eigv.bas->GetNBfunc() << 
			    " is not equal to the number of Active orbitals " << nao << std::endl;
		return false;
	}

	PrintLog(" EigenVector Expansion Coef: \n");
	for(int i=0; i < nao; i++)
	{
		std::string lbl = ptr_qc_mod->ActBas->GetLabel(i);
		PrintLog(" %d %s %12.6f \n", i, lbl.c_str(), eigv.coef.GetVal_idx0(i,idx-1));
	}
	PrintLog("*************************** \n");

	return true;
}

bool ETCouplMod::ScanEigEneField(int first_eig_val, int last_eig_val, double ifield_val, double ffield_val, double step_val)
{
	int i;

	PrintLog(" Scan eigen energies of the effective hamiltonian \n");
	PrintLog(" Applying electrical field in the direction from the Donor to Acceptor \n");
	PrintLog(" Range of eigenvalues from %6d to %6d \n",first_eig_val, last_eig_val); 
	PrintLog(" Init Field = %16.9f  Finale field = %16.9f  step= %16.9f \n", 
		       ifield_val, ffield_val,step_val);
 
	if( fabs(step_val) < DBL_EPSILON ) 
	{
		ErrorInMod("ETCouplMod::ScanEigEneField()",
			       " step value for electric field is zero ");
		return false;
	}

	if( (ffield_val - ifield_val)/step_val  < DBL_EPSILON )
	{
		ErrorInMod("ETCouplMod::ScanEigEneField()",
			       " difference between final and init field is zero or different sign from step value ");
		return false;
	}

	FILE* fp_ene = fopen("scan_ene_field.dat","w");
	if( fp_ene == NULL)
	{
		ErrorInMod("ETCouplMod::ScanEigEneField()",
			       " Can't open create file scan_ene_field.dat \n will write energies only to log ");
		return false;
	}

	PrintLog("                ");
	if(fp_ene) { fprintf(fp_ene,"                "); }

	for(i= first_eig_val; i <=last_eig_val ; i++)
	{
		PrintLog("      %5d      ",i);
		if(fp_ene) { fprintf(fp_ene,"      %5d      ",i); }
	}
	PrintLog("\n");
    if(fp_ene) { fprintf(fp_ene,"\n"); }
	 
	double efield;
	for( efield = ifield_val; efield <= ffield_val; efield+= step_val)
	{
		SetDAfield(efield);
		DiagHeff();
		
		PrintLog(" %14.9f ", efield);
		if(fp_ene) { fprintf(fp_ene," %14.9f ", efield); }
             
		for(i= first_eig_val; i <=last_eig_val ; i++)
		{
			 PrintLog(" %14.9f ", enel(i));
             if(fp_ene) { fprintf(fp_ene," %14.9f ", enel(i)); }
		}
		PrintLog("\n");
		if(fp_ene) { fprintf(fp_ene,"\n"); }
	}

    if(fp_ene) { fclose(fp_ene); }
    return true;
}

bool ETCouplMod::CreateEigVecContour(int idx, double flvl, int grid_size)
{
	HaQCMod* ptr_qc_mod = GetQCMod();
	char buf[128];
	if(ptr_qc_mod == NULL)
	{
		std::cerr << " Error in ETCouplMod::CreateEigVecContour() " << std::endl;
		std::cerr << " Qchem module is not set for ET Coupling module " << std::endl;
		return false;
	}

	if( idx > enel.size() )
	{
		std::cerr << " Error in ETCouplMod::CreateEigVecContour() " << std::endl;
		std::cerr << " Selected Eigen Vec index is larger than the number of Eigen Vec " << std::endl;
		return false;
	}
	
	if(idx > eigv.GetNOrbs() && idx <= 0)
	{
		PrintLog(" Error in ETCouplMod::CreateEigVecContour() \n");
		PrintLog(" Index of the effective orbital %d Is out of range \n", idx );
		return false;
	}
	
	int nao= ptr_qc_mod->GetNActiveOrb();

	if(eigv.bas->GetNBfunc() != nao)
	{
		PrintLog(" Error in ETCouplMod::CreateEigVecContour() \n");
		PrintLog(" The Number of basis functions for the Heff eigenvector %d",
			       eigv.bas->GetNBfunc());
	    PrintLog(" is not equal to the number of Active orbitals %d \n", nao);
		return false;
	}

	HaVec_double cf(nao);
	int i;
	for(i = 0; i < nao; i++)
	{
		cf[i] = eigv.coef.GetVal_idx0(i,idx-1);
	}

	VecPtr surfs_ptr = ptr_qc_mod->CreateOrbContour(cf, *(ptr_qc_mod->ActBas),flvl, grid_size);
	int nsurf = surfs_ptr.size();
	for(i = 0; i < nsurf; i++)
	{
		sprintf(buf,"%d",idx);
		std::string str_idx = buf;
		boost::trim(str_idx);
		
		std::string cname = "HAM_EIGV_" + str_idx;
		if( i == 0) cname += "_POS";
		if( i == 1) cname += "_NEG";
		HaDisplayedSurface* pcnt = (HaDisplayedSurface*)surfs_ptr[i];
		pcnt->SetObjName(cname.c_str());
	}
	return true;
}

double ETCouplMod::calc_edge_dist()
{
	double dist_min=1000.0;
	double dist;

	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if(!pdon || !pacc || pdon->empty() || pacc->empty())
	{
		ErrorInMod(" ETCouplMod::calc_edge_dist() ",
		          " Donor or acceptor are not set ");
		return 0.0;
	}

	HaAtom* aptr1;
	HaAtom* aptr2;
	AtomIteratorAtomGroup aitr_don(pdon);
	AtomIteratorAtomGroup aitr_acc(pacc);
	for(aptr1= aitr_don.GetFirstAtom(); aptr1; aptr1=aitr_don.GetNextAtom())
	{
		for(aptr2= aitr_acc.GetFirstAtom(); aptr2; aptr2= aitr_acc.GetNextAtom() )
		{
			dist=HaAtom::CalcDistance(aptr1,aptr2,ANGSTROM_U);
			if(dist < dist_min) dist_min = dist;
		}
	}

	PrintLog(" Donor - Acceptor edge to edge distance is %12.3f Ang \n", dist_min); 
	return dist_min;
}

bool ETCouplMod::RecalcHeff()
{
	if(ptr_qc_mod == NULL) 
	{
		ErrorInMod(" ETCouplMod::RecalcHeff() ",
		" Host molecule is not set ");
		return false;
	}
	int nb= ptr_qc_mod->GetNBfunc();
	
	if( nb <= 0 ) 
	{
		ErrorInMod(" Error in ETCouplMod::RecalcHeff()",
		   " Basis set is not set for the host molecule ");
		return false;
	}
	
	int nab= ptr_qc_mod->GetNActiveOrb();
	
	if( nab <= 0 ) 
	{
		ErrorInMod("ETCouplMod::RecalcHeff()",
		  " Active Local Orbitals are not set for the host molecule ");
		return false;
	}

	if( ptr_qc_mod->MO_coef.num_cols() == 0 || ptr_qc_mod->MOene.size() == 0)
	{
		ErrorInMod("ETHeffHF::Recalc()",
        " MOs or/and their energies are not set for the host QChem Module ");
		return false;
	}
	
	int iresult;

	heff_mat.newsize(nab,nab);

	bool bres;

	if( ptr_qc_mod->IsLocOrbFullBasis())
	{
		bres = ptr_qc_mod->BuildFockMatFromMOs(heff_mat);
		ssl = ptr_qc_mod->GetOvlpMat();
		return bres;
	}

	HaMat_double gfull(nb,nb);
	HaMat_double gl(nab,nab); // Green Function matrix on local active Orbitals

	HaMat_double& ss = ptr_qc_mod->GetOvlpMat();
	GetActBasOvlpMat();
//	ptr_qc_mod->ProjMatToActBas(ss,ssl);
	
//    ptr_qc_mod->CalcEPfromMO(gfull,tun_ene);
//	ptr_qc_mod->ProjMatToActBas(gfull,gl);

	HaMat_double ss_12, ss_c, ss_c_em;

	HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas, &(ptr_qc_mod->AtBasis),ss_12);
	matmult(ss_c,ss_12,ptr_qc_mod->MO_coef);
	
	ss_c_em = ss_c;
	
	int nmo = ptr_qc_mod->MO_coef.num_cols();
	int imo,j;

//	HaVec_double bnorm(nab,0.0);

	for(imo = 0; imo < nmo; imo++)
	{
		double scale = 1.0/(tun_ene - ptr_qc_mod->MOene[imo]);
		for(j=0; j < nab; j++)
		{
			double ssv = ss_c_em.GetVal_idx0(j,imo);
//			bnorm[j] += ssv*ssv;
			ss_c_em.SetVal_idx0(j,imo, ssv*scale);
		}
	}
	matmult_T2(gl,ss_c,ss_c_em);

//	PrintLog(" Projected Norms of Active orbital basis functions: \n");
//	for(j = 0; j < nab; j++)
//	{
//        bnorm[j] = 1/sqrt(bnorm[j]);
//		PrintLog(" %5d %10.6f \n", j+1,bnorm[j]);
//	}

//	for(i=0; i < nab; i++)
//	{
//		for(j = 0; j < nab; j++)
//		{
//			gl.SetVal_idx0(i,j, gl.GetVal_idx0(i,j)*bnorm[j]*bnorm[j]);
//		}
//	}


//		PrintLog(" Electron Propagator on local orbitals \n");
//		gl.Print_format(cout, " %10.5f ");
	
	iresult = HaMat_double::mat_inverse(gl);
	if(!iresult) 
		return false;


	HaMat_double scr;
	matmult(scr,ssl,gl);
	matmult(gl,scr,ssl);
	
	mat_scale(scr,ssl, tun_ene);
	
	mat_diff(heff_mat,scr,gl); // H_eff = E*S_PP - S_PP G_PP^-1 S_PP  
	

	return true;
}

bool ETCouplMod::PutSubMatToDB()
{
	if(!ptr_qc_mod)
	{
		PrintLog(" Error in ETHeffHF::PutSubMatToDB(): \n"); 
		PrintLog(" Host QCMod is not set \n");
		return false;
	}
	
	MolSet* phmol_set= ptr_qc_mod->GetMolSet();

	int nab=ptr_qc_mod->GetNActiveOrb();
	if(nab <= 0)
	{
		PrintLog(" Error in ETHeffHF::PutSubMatToDB(): \n"); 
		PrintLog(" Active orbitals are not set \n"); 
		return false;
	}
	
	int ngrp=phmol_set->GetNChemGroups();
	if( ngrp <= 0)
	{
		PrintLog(" Error in ETHeffHF::PutSubMatToDB(): \n"); 
		PrintLog(" Structural groups are not set for the host molecule \n"); 
		return false;
	}

	HaMatDB db_file(db_file_name.c_str(),"a");

	HaMat_double sub_mat; 

	for( int ig1=1; ig1 <= ngrp; ig1++)
	{
		for( int ig2=1; ig2 <= ngrp; ig2++) 
		{
			ChemGroup& g1= phmol_set->GetChemGroupByIdx(ig1);
			ChemGroup& g2= phmol_set->GetChemGroupByIdx(ig2); 
			HaGrpOperID gmat_id("HEFF",g1,g2);
			ptr_qc_mod->ExtractLocOrbSubMat(g1.GetID(),g2.GetID(),heff_mat,sub_mat);		
			db_file.put(gmat_id,sub_mat);
		}
	}
	db_file.close();
	return true;

}

bool ETCouplMod::GetSubMatFromDB()
{
	HaMatDB db_file(db_file_name.c_str(),"a");
	if(ptr_qc_mod == NULL)return false;

	MolSet* phmol_set= ptr_qc_mod->GetMolSet();

	int ngrp= phmol_set->GetNChemGroups(); 
	
	int nab= ptr_qc_mod->GetNActiveOrb();
	heff_mat.newsize(nab,nab);
	heff_mat=0.0;

	HaMat_double sub_mat; 

	protect_mat.newsize(ngrp,ngrp);
	protect_mat=0.0;
        HaMat_double tmp_protect;
		std::string key_prot;

	for( int ig1=1; ig1 <= ngrp; ig1++)
	{
		for( int ig2=1; ig2 <= ngrp; ig2++) 
		{
			ChemGroup& g1= phmol_set->GetChemGroupByIdx(ig1);
			ChemGroup& g2= phmol_set->GetChemGroupByIdx(ig2); 
			HaGrpOperID gmat_id("HEFF",g1,g2);
			
			if(!db_file.get(gmat_id,sub_mat))
			{
				if(debug_level > 10)
				{
					std::cout << "not found: " << gmat_id.GetStr() << std::endl;
				}
				continue;
			}
			if(debug_level > 10)
			{
				std::cout << "found: " << gmat_id.GetStr() << std::endl;
			}
			ptr_qc_mod->InsertLocOrbSubMat(g1.GetID(),g2.GetID(),heff_mat,sub_mat);	
                        key_prot = gmat_id.GetStr();
                        key_prot += "_prot";
                        db_file.GetMat(key_prot.c_str(),tmp_protect);
                        if(tmp_protect.num_rows() == 1 && tmp_protect.num_cols() == 1)
                        {
			   protect_mat(ig1,ig2)= tmp_protect.GetVal_idx0(0,0);
                        }
                        else
                        {
                           protect_mat(ig1,ig2)=0.0;
                        }
		}
	}
	db_file.close();
	return true;
}

bool ETCouplMod::DiagHeff()
{
	if( ptr_qc_mod == NULL)
	{
		ErrorInMod("ETCouplMod::DiagHeff()",
		   " QChem module is not set for ETCouplMod module ");
		return false;
	}

	int nab= ptr_qc_mod->GetNActiveOrb();
	if( nab <= 0 )
	{
		ErrorInMod("ETCouplMod::DiagHeff()",
		 " Active Local Orbitals are not set");
		return false;		
	}
	if( nab != heff_mat.num_rows() ||  nab != heff_mat.num_cols())
	{
		PrintLog(" Error in ETCouplMod::DiagHeff() \n");
		PrintLog(" The size of local orbitals basis: %d \n", nab);
		PrintLog(" doesn't correspond to the size of Heff %d X %d \n ",
			     heff_mat.num_rows(), heff_mat.num_cols() );
		return false;		
	}

	HaMat_double  hm;

	if(extern_field.num_rows()  == heff_mat.num_rows() && 
		extern_field.num_cols() == heff_mat.num_cols())
	{
		mat_add(hm, heff_mat, extern_field);
	}
	else
	{
		PrintLog(" ETHeffHF::Diagonalize(): \n"); 
		PrintLog(" Dimension of External field matrix doesn't match that of heff_mat \n");
		PrintLog(" External field in not added \n");
		hm= heff_mat;
	}

	HaMat_double ss_empty;
	
	if( ptr_qc_mod->wave_fun_type != harlem::qc::NDO)
	{
		HaMat_double& ssl_1 = GetActBasOvlpMat();
		HaMat_double::DiagMat(hm,ssl_1,eigv.coef,enel);
	}
	else
	{
		HaMat_double::DiagMat(hm,ss_empty,eigv.coef,enel);
	}

	eigv.bas = ptr_qc_mod->ActBas;
			
	return true;
	
}

int ETCouplMod::ZeroLongInter(double cutoff)
{
	if( ptr_qc_mod == NULL) return FALSE;
	int nab= ptr_qc_mod->GetNActiveOrb();
	if( nab == 0 || nab != heff_mat.num_rows() ||  nab != heff_mat.num_cols())
	{
		return FALSE;
	}

	Vec3DValArray orb_at(nab);

	int i,j;
	for( i = 0 ; i < nab; i++)
	{
		Vec3D* pt = ptr_qc_mod->ActBas->GetHostPt(i);
		if(pt != NULL)
			orb_at[i] = *pt;
	}

	if(nab != DA_dipole.num_rows() || nab != DA_dipole.num_cols()) 
	{
		bool result=this->SetDAdipoleMat();
		if(!result) return false;
	}

	for( i = 0; i < nab; i++)
	{
		for( j = 0; j < nab; j++)
		{
			double dist = Vec3D::CalcDistance(&orb_at[i-1],&orb_at[j-1]);
			if( dist > cutoff) 
			{
				heff_mat.SetVal_idx0(i,j,0.0);
				if(ham_trunc_type == 0)
				{
					ssl.SetVal_idx0(i,j,0.0);
					DA_dipole.SetVal_idx0(i,j,0.0);
				}
			}
		}		
	}
	
	return TRUE;
}

int ETCouplMod::SaveHeffXml(FILE* file_out)
{
	if( file_out == NULL ) 
	{
		PrintLog("Error in ETCouplMod::SaveHeffXml(FILE* file_out) \n");
		PrintLog("file_out == NULL \n");
		return FALSE;
	}

	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );

	TiXmlElement* root_element = new TiXmlElement("HARLEM_DATA");
		
	doc.LinkEndChild(root_element);

	heff_mat.AddXml( root_element, "heff_mat" );

	if(ptr_qc_mod)
	{
		ptr_qc_mod->ActBas->AddXml(root_element, "basis" );
	}
	
	bool bres = doc.SaveFile( file_out );
	
	if( !bres ) return FALSE;

	return TRUE;
}

int ETCouplMod::LoadFragmHeffXml(FILE* file_inp)
{
	TiXmlDocument doc;
	bool bres = doc.LoadFile(file_inp);

	const TiXmlElement* root_element;

	if(bres)
	{
		root_element = doc.FirstChildElement();
		if( root_element == NULL )
		{
			PrintLog("Error in ArrayOrb3D::LoadXml() \n");
			PrintLog("No ROOT Element \n");
			return FALSE;
		}
	}	
	const TiXmlElement* heff_element  = root_element->FirstChildElement("heff_mat");
	const TiXmlElement* basis_element = root_element->FirstChildElement("basis");

	if( heff_element == NULL )
	{
		PrintLog("Error in ETCouplMod::LoadFragmHeffXml() \n");
		PrintLog("No heff_mat element \n");
		return FALSE;
	}

	if( basis_element == NULL )
	{
		PrintLog("Error in ETCouplMod::LoadFragmHeffXml() \n");
		PrintLog("No heff_mat element \n");
		return FALSE;
	}
	
	HaMat_double heff_frag;
	heff_frag.LoadXml(heff_element);

	GauBasisSet basis_frag;
	basis_frag.LoadXml(basis_element);

	int nb_frag = basis_frag.GetNBfunc();

	if( heff_frag.num_cols() != nb_frag || heff_frag.num_rows() != nb_frag )
	{
		PrintLog("\nError in ETCouplMod::LoadFragmHeffXml() \n");
		PrintLog("Dimentions of the the fragment hamiltonian matrix doesn't correspond to the basis set \n");
		return FALSE;
	}
	
	if( this->ptr_qc_mod->ActBas == NULL)
	{
		PrintLog("Error in ETCouplMod::LoadFragmHeffXml() \n");
		PrintLog("No Active Basis Set \n");
		return FALSE;
	}

	std::string bas_type = this->ptr_qc_mod->ActBas->GetClassName();
	
	if( bas_type != "GauBasisSet")
	{
		PrintLog("Error in ETCouplMod::LoadFragmHeffXml() \n");
		PrintLog("Active Basis is not GauBasisSet - not supported yet \n");
		return FALSE;
	}

	GauBasisSet* act_bas = (GauBasisSet*) this->ptr_qc_mod->ActBas;
    int nb = act_bas->GetNBfunc();

	if(heff_mat.num_rows() != nb || heff_mat.num_cols() != nb ) 
	{
		heff_mat.resize(nb,nb);
		heff_mat.set(0.0);
	}

	basis_frag.RecompFstBasVec();
	act_bas->RecompFstBasVec();

	IntIntMap frag_bas_fun_map;

	HaVec_double frag_bas_pert_vec;

	HaVec_double*  ptr_bas_pert_vec = NULL;

	if( use_pert_mat)
	{
		ptr_bas_pert_vec = &frag_bas_pert_vec;
		if( heff_pert_mat.num_rows() != nb  || heff_pert_mat.num_cols() != nb )
		{
			heff_pert_mat.resize(nb,nb);
			heff_pert_mat = 1.0;
		}
	}

	act_bas->MatchBasisSet(&basis_frag,frag_bas_fun_map,ptr_bas_pert_vec );
	
	IntIntMap::iterator itr1,itr2;

	for( itr1 = frag_bas_fun_map.begin(); itr1 != frag_bas_fun_map.end(); itr1++)
	{
		int j1 = (*itr1).first;
		int i1 = (*itr1).second;
		
		for( itr2 = frag_bas_fun_map.begin(); itr2 != frag_bas_fun_map.end(); itr2++)
		{
			int j2 = (*itr2).first;
			int i2 = (*itr2).second;
			
			int assign_elem = 1;

			if( use_pert_mat )
			{
				double pert = sqrt( frag_bas_pert_vec[j1]* frag_bas_pert_vec[j2] );
				if( pert < heff_pert_mat.r0(i1,i2)*1.001 )
				{
					heff_pert_mat.r0(i1,i2) = pert;
					heff_pert_mat.r0(i2,i1) = pert;
				}
				else
				{
					assign_elem = 0;
				}
			}

			if( assign_elem) heff_mat.r0( i1, i2 ) = heff_frag.r0(j1,j2);
		}
	}

	return TRUE;
}

bool ETCouplMod::CalcGFDonAccOrbHeff()
{
	int i,j;
	int nab = ptr_qc_mod->GetNActiveOrb();
    int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();
	int nv = ndon + nacc;

	if( nab != heff_mat.num_rows() || nab != heff_mat.num_cols() )
	{
		PrintLog(" Error in ETHeffHF::CalcGFDonAccOrbHeff() \n"); 
		PrintLog(" The dimensions of Heff matrix %4d X %4d \n",
			       heff_mat.num_rows(),heff_mat.num_cols());
		PrintLog(" Doesn't correspond to number of active orbitals %5d \n",nab); 
		return false;
	}

	HaMat_double es(nab,nab,0.0);
	HaMat_double& ss = GetActBasOvlpMat();

	if(ptr_qc_mod->wave_fun_type != harlem::qc::NDO)
	{
		mat_scale(es,ss,tun_ene);	
	}
	else
	{
		mat_add_unit(es,es,tun_ene);
	}

	HaMat_double es_h;
	mat_diff(es_h,es,heff_mat);  // gf = (ES-H)

// Compute Green Function solving System of linear equations:
	
  
	LinCombOrb3D don_acc_orbs;

	don_acc_orbs.CreateEmptyOrbs(nv,ptr_qc_mod->ActBas);

	HaMat_double tmp;

	donor_orbs.ProjectToBasis(tmp,ptr_qc_mod->ActBas);

	for(j=0; j < ndon; j++)
	{
		for(i=0; i < nab; i++)
		{
			don_acc_orbs.coef.r0(i,j) = tmp.r0(i,j); 
		}
	}
    
	acc_orbs.ProjectToBasis(tmp,ptr_qc_mod->ActBas);

	for(j=0; j < nacc; j++)
	{
		int jacc = j + ndon;
		for(i=0; i < nab; i++)
		{
			don_acc_orbs.coef.r0(i,jacc) = tmp.r0(i,j); 
		}
	}

//   don_acc_orbs.AddOrbs(donor_orbs);
//	don_acc_orbs.AddOrbs(acc_orbs);

	if(ptr_qc_mod->wave_fun_type != harlem::qc::NDO)
	{
		matmult(tmp,ss,don_acc_orbs.coef);  // S*(D,A)
        don_acc_orbs.coef= tmp;
	}

	HaMat_double vec_da = don_acc_orbs.coef; // vec_da = S*(D,A) 

	HaMat_double::solv_lin_syst_1(es_h,vec_da); // vec_da = (ES-H)^(-1)*S*(D,A)

	matmult_T1(don_acc_gf,don_acc_orbs.coef,vec_da);    // don_acc_gf = (D,A)^T*S*(ES-H)^(-1)*S*(D,A)

	PrintLog(" Green Function elements between donor/acceptor orbitals: \n");
	PrintLog(" labels of Target orbital ");
	for( i = 0; i < ndon; i++)
	{
		PrintLog("   %s    ", (donor_orbs.GetLabel(i)).c_str());
	}
	PrintLog("\n");
	
	for( j=0; j < nacc ; j++)
	{
		PrintLog(" %s  ", (acc_orbs.GetLabel(j)).c_str() );
		for(i= 0; i < ndon; i++)
		{
//			double dist=HaAtom::CalcDistance(aptr1,aptr2,ANGSTROM_U);
//			PrintLog("  %16.9f  %16.9f \n", dist, gf_mat(idx_src[i],idx_tgt[j]) );
            PrintLog("  %16.9f ", don_acc_gf.GetVal_idx0(i,j+ndon) );
//            PrintLog("  %16.9f ", don_acc_gf.GetVal_idx0(i,j) );
		}
		PrintLog("\n");
	}

	return true;
}

int ETCouplMod::PrintOvlpElem()
{
	int nab = ptr_qc_mod->GetNActiveOrb();
	int n1 = donor_orbs.GetNOrbs();
	int n2 = acc_orbs.GetNOrbs();

	int i,j;
	HaMat_double ovlp_mat;
	HaBasisSet::CalcOvlpMat(&donor_orbs,&acc_orbs,ovlp_mat);

	PrintLog(" Active Orbitals Overlap Matrix elements: \n\n"); 
	PrintLog("            ");
	for(i = 0; i < n1; i++)
	{
		PrintLog( "  %s  ", (donor_orbs.GetLabel(i)).c_str() );
	}
	PrintLog("\n");

	for(j = 0; j < n2; j++)
	{
		for(i = 0; i < n1; i++)
		{
			if( i == 0)
			{
				PrintLog( "%s ",(donor_orbs.GetLabel(j)).c_str() );
			}
			PrintLog("  %12.6e ", ovlp_mat.GetVal_idx0(i,j) );			
		}
		PrintLog("\n");
	}

	return TRUE;
}

int ETCouplMod::PrintHeffElem()
{
	int nab = ptr_qc_mod->GetNActiveOrb();
	if( heff_mat.num_rows() != nab || heff_mat.num_cols() != nab )
	{
		PrintLog(" Error in: ETCouplMod::PrintHeffElem() \n");
		PrintLog(" H_eff matrix on active orbitals is not set \n");
		return FALSE;
	}

	return TRUE;
}

bool ETCouplMod::CalcHDAfromGF()
{
	if(ptr_qc_mod->ActBas == NULL) return false;

	int nab = ptr_qc_mod->ActBas->GetNBfunc();

	int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();
	int nv = ndon + nacc;
	
	if(ndon == 0 || nacc == 0)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAfromGF() \n");
		PrintLog(" No donor or acceptor orbitals set \n");
		return false;
	}

	if(nab == heff_mat.num_rows() && nab == heff_mat.num_cols())
	{
		CalcGFDonAccOrbHeff();
	}
	else
	{
		CalcGFDonAccOrb();
	}

	if(don_acc_gf.num_rows() != nv || don_acc_gf.num_rows() != nv)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAfromGF() \n");
		PrintLog(" Error computing Green function between donor and acceptor orbitals \n");
		return false;
	}

	HaMat_double ss_da_1;

	HaBasisSet::CalcOvlpMat(&donor_orbs,&acc_orbs,ss_da_1);

	HaMat_double ss_da(nv,nv);
	ss_da = 0.0;
	int i,j;
	for( i = 0; i < ndon; i++)
	{
		ss_da.SetVal_idx0(i,i,1.0);
	}
	for( i = ndon; i < nv; i++)
	{
		ss_da.SetVal_idx0(i,i,1.0);
	}
	for(i = 0; i < ndon; i++)
	{
		for(j = 0; j < nacc; j++)
		{
			double sij = ss_da_1.GetVal_idx0(i,j);
			ss_da.SetVal_idx0(i,j+ndon,sij);
			ss_da.SetVal_idx0(j+ndon,i,sij);
		}
	}

	HaMat_double gf = don_acc_gf;

	int iresult=HaMat_double::mat_inverse(gf);
	if(!iresult)
	{
		ErrorInMod("ETCouplMod::CalcHDAfromGF()",
			       " Failed to invert GF matrix on donor and acceptor vectors");
		return false;
	}

	HaMat_double hda,scr;

	matmult(scr,ss_da,gf);
	matmult(gf,scr,ss_da);
	
	HaMat_double es;
	mat_scale(es,ss_da,tun_ene);
	
	mat_diff(hda,es,gf); // H_DA = E*S_DA - S_DA G_DA^-1 S_DA
	
	PrintLog(" Diagonal Energies of the donor orbitals: \n"); 
	for(i = 0; i < ndon; i++)
	{
		PrintLog(" %d %12.6f \n",i+1,hda.GetVal_idx0(i,i));
	}

	PrintLog(" Diagonal Energies of the acceptor orbitals: \n"); 
	for(i = 0; i < nacc; i++)
	{
		PrintLog(" %d %12.6f \n",i+1,hda.GetVal_idx0(i+ndon,i+ndon));
	}
	
	da_coupl_val.newsize(ndon,nacc);

	PrintLog(" Effective coupling between donor and acceptor orbitals (a.u.) is \n");
	for(j = 0; j < nacc; j++)
	{
		for(i = 0; i < ndon; i++)
		{	
			double tda = hda.GetVal_idx0(i,j+ndon);
			da_coupl_val.SetVal_idx0(i,j,tda);
			PrintLog(" %9.3e ",tda);

		}
		PrintLog("\n");
	}
	
	return true;
	
}

bool ETCouplMod::CalcHDAPert(double &hda_coupl)
{
	if(ptr_qc_mod == NULL)
		return false;

	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if(!pdon || !pacc || pdon->empty() || pacc->empty() )
	{
		ErrorInMod("ETCouplMod::CalcHDAPert()", 
		           " Donor or acceptor atom groups are not defined ");
		return false;
	}

	int ndon = donor_orbs.GetNOrbs();
	int nacc = acc_orbs.GetNOrbs();

	if( ndon == 0 || nacc == 0)
	{
		PrintLog(" Error in ETCouplMod::CalcHDAPert() \n");
		PrintLog(" Donor or acceptor orbitals are not set \n");
		return false;
	}

	int nab = ptr_qc_mod->GetNActiveOrb();

	if( nab != heff_mat.num_rows() || nab != heff_mat.num_cols() )
	{
		PrintLog(" Error in ETHeffHF::CalcHDAPert() ");
		PrintLog(" The dimensions of Heff matrix %d X %d \n",
			      heff_mat.num_rows(), heff_mat.num_cols() );
		PrintLog( " Doesn't correspond to the number of active orbitals %d ", nab);		return false;
	}

	std::vector<int> bridge_orb_idx;
	std::vector<int> donor_orb_idx;
	std::vector<int> acceptor_orb_idx;

	int i,j;
	for(i=1; i <= nab; i++) 
	{
		Vec3D* hpt = ptr_qc_mod->ActBas->GetHostPt(i-1);
		const HaAtom* pat_host= (const HaAtom*) hpt;
		if(pdon->HasAtom(pat_host))
		{		
			donor_orb_idx.push_back(i);
		}
		else if(pacc->HasAtom(pat_host))
		{		
			acceptor_orb_idx.push_back(i);
		}
		else
		{
			bridge_orb_idx.push_back(i);
		}
	}

	int nbr_orb = bridge_orb_idx.size();
	int ndon_orb = donor_orb_idx.size();
	int nacc_orb = acceptor_orb_idx.size(); 

	HaMat_double don_vec_short(ndon_orb,ndon);
	HaMat_double acc_vec_short(nacc_orb,nacc);

	HaMat_double es_h_br(nbr_orb,nbr_orb);
	HaMat_double es_h_dbr(ndon_orb,nbr_orb);
	HaMat_double es_h_abr(nacc_orb,nbr_orb);

	int id,ia;

	for(id = 0; id < ndon; id++)
	{
		for(i = 0; i < ndon_orb; i++)
		{
			don_vec_short.SetVal_idx0(i,id, donor_orbs.coef.GetVal_idx0(donor_orb_idx[i]-1,id));
		}
	}

	for(ia = 0; ia < nacc; ia++)
	{
		for(i = 0; i < nacc_orb; i++)
		{
			acc_vec_short.SetVal_idx0(i,ia, acc_orbs.coef.GetVal_idx0(acceptor_orb_idx[i]-1,ia));
		}
	}

    HaMat_double ss_et(nab,nab,0.0);
	if(ptr_qc_mod->wave_fun_type != harlem::qc::NDO)
	{
		mat_scale(ss_et,ssl,tun_ene);
	}
	else
	{
		mat_add_unit(ss_et,ss_et,tun_ene);
	}
		
	for( i = 0; i < nbr_orb; i++ )
	{
		for( j = 0; j < nbr_orb; j++ )
		{
			es_h_br(i+1,j+1) = ss_et(bridge_orb_idx[i],bridge_orb_idx[j])-
			                   heff_mat(bridge_orb_idx[i],bridge_orb_idx[j]);
		}
		for( j = 0; j < ndon_orb; j++)
		{
			if(set_dab_huck_inter)
			{
				es_h_dbr(j+1,i+1) = -0.2*ss_et(donor_orb_idx[j],bridge_orb_idx[i]);
			}
			else
			{
				es_h_dbr(j+1,i+1) = ss_et(donor_orb_idx[j],bridge_orb_idx[i])-
			                        heff_mat(donor_orb_idx[j],bridge_orb_idx[i]);
			}
		}
		for( j = 0; j < nacc_orb; j++)
		{
			if(set_dab_huck_inter)
			{
			   es_h_abr(j+1,i+1) = -0.2*ss_et(acceptor_orb_idx[j],bridge_orb_idx[i]);
			}
			else
			{
			   es_h_abr(j+1,i+1) = ss_et(acceptor_orb_idx[j],bridge_orb_idx[i])-
			                       heff_mat(acceptor_orb_idx[j],bridge_orb_idx[i]);
			}
		}
	}

	int iresult = HaMat_double::mat_inverse(es_h_br); // (ES_b-H_b)^(-1) 
	if(!iresult)
	{
		ErrorInMod("ETCouplMod::CalcHDAPert()",
			       " Failed to invert (ES-H)_br matrix ");
		return false;
	}

    HaMat_double scr;
    HaMat_double tmp_mat;

    matmult_T1(scr,don_vec_short,es_h_dbr); 
	matmult(tmp_mat,scr,es_h_br);  // D^T X (ES-H)_db X (ES_b-H_b)^(-1)
	matmult_T2(scr,tmp_mat,es_h_abr); 
	matmult(tmp_mat,scr,acc_vec_short); // D^T X (ES-H)_db X (ES_b-H_b)^(-1) X
	                                    //       (ES-H)_ab X A 

	da_coupl_val.newsize(ndon,nacc);

	PrintLog(" D/A coupling by the perturbation formula are:\n");
	for(j = 0; j < nacc; j++)
	{
		for(i = 0; i < ndon; i++)
		{	
			double tda = tmp_mat.GetVal_idx0(i,j);
			da_coupl_val.SetVal_idx0(i,j,tda);
			PrintLog(" %9.3e ",tda);

		}
		PrintLog("\n");
	}

	return true;
}

bool ETCouplMod::PrintProtectMat() const
{
	MolSet* phmol_set= ptr_qc_mod->GetMolSet();
	if(phmol_set == NULL)
	{
		ErrorInMod(" ETCouplMod::PrintProtectMat()",
		   " Host Molecule is not set ");
		return false;
	}
	int ngrp= phmol_set->GetNChemGroups(); 
	if( ngrp == 0)
	{
		ErrorInMod(" ETCouplMod::PrintProtectMat()",
		   "Number of groups in the host Molecule is equal zero ");
		return false;
	}

	if( protect_mat.num_rows() != ngrp || protect_mat.num_cols() != ngrp)
	{
		ErrorInMod(" ETCouplMod::PrintProtectMat()", 
		" Size of the protection matrix doesn't correspond to the number of groups ");
		return false;
	}

	int ig1;

	for( ig1=1; ig1 <= ngrp; ig1++)
	{
		ChemGroup& g1= phmol_set->GetChemGroupByIdx(ig1);
		PrintLog(" %s ",g1.GetID());
		PrintLog(" %6.4f \n", protect_mat(ig1,ig1));
	}
	PrintLog("\n");

// Print Diagonal:

//	for( ig1=1; ig1 <= ngrp; ig1++)
//	{
//		sprintf(buf," %6.4f ", protect_mat(ig1,ig1));
//		sout << buf;
//	}
//	sout << std::endl;


//	for( ig1=1; ig1 <= ngrp; ig1++)
//	{
//		ChemGroup& g1= phmol->GetChemGroupByIdx(ig1);
//		sprintf(buf," %5d ",g1.GetID());
//		sout << buf;
		
//		for( int ig2=1; ig2 <= ngrp; ig2++)
//		{
//			sprintf(buf," %6.4f ",protect_mat(ig1,ig2));
//			sout << buf;
//		}
//		sout << std::endl;
//	}
	return true;
}

bool ETCouplMod::calc_intermol_path_coupl()
{
	MolSet* pmset = GetMolSet();
	AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

	if(!pdon || !pacc || pdon->empty() || pacc->empty() )
	{
		ErrorInMod("ETCouplMod::calc_intermol_path_coupl()", 
		           " Donor or acceptor are not defined ");
		return false;
	}

	int nmol = pmset->GetNMol();

	if(nmol < 2)
	{
		ErrorInMod("ETCouplMod::calc_intermol_path_coupl()",
			       "Number of molecules in the set is less than 2");
		return false;
	}

	HaAtom* aptr1 = (*pdon)[0];
	HaAtom* aptr2 = (*pacc)[0];

	HaMolecule* pmol1 = aptr1->GetHostMol();
	HaMolecule* pmol2 = aptr2->GetHostMol();
	
	if(pmol1 == pmol2)
	{
		ErrorInMod("ETCouplMod::calc_intermol_path_coupl()",
			       "Both donor and acceptor belong to the same molecule");
		return false;
	}

	PtrDoubleMap::iterator mitr1; 
	PtrDoubleMap::iterator mitr2; 
	int ires;

    if(rebuild_mol_coupl_map || mol1_coupl_map.empty() || mol2_coupl_map.empty())
	{
		if( pw_nb_decay_intermol < 0.0001)
		{
			pw_nb_decay_intermol = pw_nb_decay;
		}	
		mol1_coupl_map.clear();
		mol2_coupl_map.clear();

		pathways_calc_type=COUPL_MAP;
		path_coupl_calc();
		
		int n,i;
		
		n =nodes.size();	
		HaAtom* aptr;
		for(i=0; i < n; i++)
		{
			aptr=nodes[i];
			if(aptr->GetHostMol() == pmol1)
			{
				aptr->tempf = coupl_map[i].coupling;
				mol1_coupl_map[aptr] = coupl_map[i].coupling;
			}
		}
		
		pdon->swap(*pacc);
		
		pathways_calc_type = COUPL_MAP;
		path_coupl_calc();
		for(i=0; i < n; i++)
		{
			aptr=nodes[i];
			if(aptr->GetHostMol() == pmol2)
			{
				aptr->tempf = coupl_map[i].coupling;
				mol2_coupl_map[aptr] = coupl_map[i].coupling;
			}
		}

		for( mitr1 = mol1_coupl_map.begin(); mitr1 != mol1_coupl_map.end() ; mitr1++ )
		{
			aptr1 = (HaAtom*) mitr1.GetKey();
			aptr1->tempf =    mitr1.GetVal();
		}
		
		pdon->swap(*pacc);
		rebuild_mol_coupl_map = FALSE;

		// assign coupling to Hydrogen atoms equal to PATHWAYS 
		// coupling to heteroatoms they bonded to
		HaBond* bptr;
		BondIteratorMolSet bitr(phost_mset);
		for(bptr= bitr.GetFirstBond(); bptr; bptr= bitr.GetNextBond())
		{
			if(bptr->srcatom->GetElemNo() == 1) bptr->srcatom->tempf = bptr->dstatom->tempf;
			if(bptr->dstatom->GetElemNo() == 1) bptr->dstatom->tempf = bptr->srcatom->tempf;
		}
	}

	double coupl1,coupl2;
	
	double avg_crd_1[3]={0.0,0.0,0.0};
	double avg_crd_2[3]={0.0,0.0,0.0};

	AtomIteratorAtomGroup aitr1(pdon);
	AtomIteratorAtomGroup aitr2(pacc);

	for(aptr1 = aitr1.GetFirstAtom(); aptr1; aptr1 = aitr1.GetNextAtom())
	{
		avg_crd_1[0] += aptr1->GetX();
		avg_crd_1[1] += aptr1->GetY();
		avg_crd_1[2] += aptr1->GetZ();
	}
	avg_crd_1[0] /= pdon->size();
	avg_crd_1[1] /= pdon->size();
	avg_crd_1[2] /= pdon->size();

	for(aptr2 = aitr2.GetFirstAtom(); aptr2; aptr2 = aitr2.GetNextAtom())
	{
		avg_crd_2[0] += aptr2->GetX();
		avg_crd_2[1] += aptr2->GetY();
		avg_crd_2[2] += aptr2->GetZ();
	}
	avg_crd_2[0] /= pacc->size();
	avg_crd_2[1] /= pacc->size();
	avg_crd_2[2] /= pacc->size();

	int i;
	double dist2_da = 0.0;

	for(i = 0; i < 3; i++)
	{
		dist2_da += (avg_crd_2[i] - avg_crd_1[i])*(avg_crd_2[i] - avg_crd_1[i]);
	}

//        PrintLog("dist2_da = %12.6f \n",dist2_da);

	double coupl_max = pw_ln_cov_decay - pw_nb_decay_intermol* ( sqrt(dist2_da) - pw_nb_min_dist);

	HaAtom* at_cont_1 = NULL;
	HaAtom* at_cont_2 = NULL;
	double coupl_jump = 0.0;

	for( mitr1 = mol1_coupl_map.begin(); mitr1 != mol1_coupl_map.end() ; mitr1++ )
	{
		aptr1 = (HaAtom*) mitr1.GetKey();
		coupl1 = mitr1.GetVal();
//		coupl1 = aptr1->tempf;

		double dist2_to_acc = 0.0;

		dist2_to_acc += (aptr1->GetX() - avg_crd_2[0])*(aptr1->GetX() - avg_crd_2[0]);
		dist2_to_acc += (aptr1->GetY() - avg_crd_2[1])*(aptr1->GetY() - avg_crd_2[1]);
		dist2_to_acc += (aptr1->GetZ() - avg_crd_2[2])*(aptr1->GetZ() - avg_crd_2[2]);

		if( dist2_to_acc > dist2_da) continue;

		for( mitr2 = mol2_coupl_map.begin(); mitr2 != mol2_coupl_map.end() ; mitr2++ )
		{
			aptr2 = (HaAtom*) mitr2.GetKey();
			coupl2 = mitr2.GetVal();
//			coupl2 = aptr2->tempf;

			double cur_coupl = coupl1 + coupl2;

			if( cur_coupl < coupl_max) continue;

			double c_jump = pw_ln_cov_decay - pw_nb_decay_intermol* ( Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U) - pw_nb_min_dist); ;
			cur_coupl += c_jump;
			if( cur_coupl > coupl_max) 
			{
				at_cont_1 = aptr1;
				at_cont_2 = aptr2;
				coupl_jump = c_jump;
				coupl_max = cur_coupl;
			}
		}
	}

	best_path_coupl = exp(coupl_max);

	char buf[256];

	if(debug_level >=5 )
	{
		PrintLog(" Best PATHWAYS coupling between two molecules %12.6e \n", best_path_coupl);
		if(aptr1 != NULL && aptr2 != NULL)
		{
			PrintLog("Best Contact is between atoms: \n ");
			at_cont_1->FillRef(buf); 
			PrintLog(" %s coupling1= %12.6e  \n ",
			                              buf,exp(mol1_coupl_map[at_cont_1]));
			at_cont_2->FillRef(buf); 
  			PrintLog(" %s coupling2= %12.6e  \n ",
			                              buf,exp( mol2_coupl_map[at_cont_2]));
			PrintLog("Coupling through the jump %12.6e ",exp(coupl_jump));
		}
	}
	
	return true;
}
