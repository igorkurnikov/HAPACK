/*! \file moleditor.cpp
 
    Classes to Edit Molecular Structures 
   
    \author Igor Kurnikov  
    \date 2008- 
  
*/

#include "harlemapp.h"
#include "habond.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "haresdb.h"
#include "hamolview.h"

#include <math.h>
#include <float.h>

#include <stdexcept>

const int NoLadder   = 0x00;
const int ParaLadder = 0x01;
const int AntiLadder = 0x02;

const double Cos70Deg = 0.34202014332567;

static int MemSize;

MolEditor::MolEditor()
{
	Init();
}

void MolEditor::Init()
{
	max_bond_length = 2.4;
	max_hbond = 3.8;
	max_da_dist_no_acc = 3.4; // max D/A dist for H -Bond including N,O acceptor atoms
	max_da_dist_s_acc =  3.7; // max D/A dist for H -bond including S atoms
	max_ha_dist_no_acc = 2.5; // max H/A dist for H -Bond including N,O acceptor atoms
	max_ha_dist_s_acc  = 2.8; // max H/A dist for H -Bond including S acceptor atoms
	max_hda_angle = 40.0 * DEG_TO_RAD;    // max H - donor - Acceptor angle in radians (40 Degrees)
	
	m_calc_s_hbonds_flag = false;

	solv_name = "water";
    solv_buffer_dist = 6.0;
}

MolEditor::~MolEditor()
{

}




int MolEditor::FindHBondsAtomCollection(AtomContainer* p_at_coll, vector<HaHBond>& hbonds )
{	
	hbonds.clear();

	AtomGroup don_atoms;
	vector< AtomGroup > h_don_atoms;

	int na_tot = p_at_coll->GetNAtoms();

	don_atoms.reserve(na_tot);
	h_don_atoms.reserve(na_tot);

	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	

	BoxPartition acc_part; // Acceptor atoms distributed on a grid
	double xmin, xmax, ymin, ymax, zmin, zmax;
	p_at_coll->GetMinMaxCrd( xmin, ymin, zmin, xmax, ymax, zmax);
	
	xmin -= 0.5;
	ymin -= 0.5;
	zmin -= 0.5;
	xmax += 0.5;
	ymax += 0.5;
	zmax += 0.5;

	acc_part.SetBoundaries(xmin, ymin, zmin, xmax, ymax, zmax);

	AtomIteratorGen aitr(p_at_coll);
	HaAtom* aptr;
	MolEditor mol_editor;

	vector<HaAtom> add_hydrogens; // Additional hydrogen attached to donor atoms to use in Hydrogen bond calculations
	AtomGroup empty_at_array;

	int i;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{	
		if( aptr->GetElemNo() == 1 || aptr->GetElemNo() == 6 ) continue;
		bool is_donor = false;
		bool is_acc   = false;
		HaAtom* aptr_templ = p_res_db->GetTemplateForAtom(aptr);
		if( aptr->IsHBDonor() )  is_donor = true;
		if( aptr_templ != NULL && aptr_templ->IsHBDonor() ) is_acc = true;
		if( aptr->IsHBAcceptor() )  is_acc = true;
		if( aptr_templ != NULL && aptr_templ->IsHBAcceptor() ) is_acc = true;

		if( is_acc ) acc_part.AddPoint( aptr );
		
		if( is_donor) 
		{	
			don_atoms.push_back( aptr );
			h_don_atoms.push_back(empty_at_array);
			AtomGroup& h_at = h_don_atoms.back(); //!< hydrogen atoms that belong to the current donor atom
			
			AtomGroup bonded_atoms;
			aptr->GetBondedAtoms( bonded_atoms );
			vector<HaAtom> extra_bonded_h = mol_editor.FindMissingHydrogens(aptr);
			
			for(i = 0; i < bonded_atoms.size(); i++)
			{
				if( bonded_atoms[i]->GetElemNo() == 1)
				{
					h_at.push_back(bonded_atoms[i]); 
				}
			}

			for(i = 0; i < extra_bonded_h.size(); i++)
			{
				add_hydrogens.push_back( extra_bonded_h[i] );
				h_at.push_back( &(add_hydrogens.back()) );
			}
		}
	}
	
	acc_part.SetRegionRad(this->max_hbond);

	for(i = 0; i < don_atoms.size(); i++ )
	{
		HaAtom* aptr_d = don_atoms[i];
		AtomGroup acc_atoms;
		AtomIteratorAtomGroup aitr_acc(&acc_atoms);
		acc_part.GetNeighbors(*aptr_d, acc_atoms );
		acc_atoms.DeleteAtom(aptr_d);
		if(acc_atoms.empty())
			continue;
		
		HaAtom* acc_ptr;
		for(acc_ptr=aitr_acc.GetFirstAtom(); acc_ptr; acc_ptr = aitr_acc.GetNextAtom())
		{	
			double da_dist_criteria = this->max_da_dist_no_acc;
			if(acc_ptr->GetElemNo() == 16)
				da_dist_criteria = this->max_da_dist_s_acc;

			double da_dist = Vec3D::CalcDistance(acc_ptr,aptr_d);

			if( da_dist > da_dist_criteria) 
				continue;

            HaAtom* haptr;
			int jh;
			int nh = h_don_atoms[i].size();

			for( jh = 0; jh < nh; jh++)
			{
				haptr = h_don_atoms[i][jh];
//				double ha_dist = Vec3D::CalcDistance(acc_ptr,haptr);
//				
//				double ha_dist_criteria = this->max_ha_dist_no_acc;
//				if( acc_ptr->GetElemNo() == 16 )
//					ha_dist_criteria = this->max_ha_dist_s_acc;
//				
//				if( ha_dist >  ha_dist_criteria)
//					continue;
				
				double hda_angle = Vec3D::CalcAngle(haptr,aptr_d,acc_ptr);

				if( hda_angle < this->max_hda_angle)
				{
					hbonds.push_back( HaHBond(aptr_d,acc_ptr) );
					break;
				}
			}
		}
	}
	return TRUE;
}

bool MolEditor::CalcHBonds(MolSet* pmset, bool recalc)
{
	bool do_calc = false;

	if( recalc || !pmset->HBonds_found )
	{
		do_calc = true;
	}
	
	if(!do_calc)
		return false;
	
	pmset->HBonds.clear();
	pmset->HBonds_found = false;
	
	SetStdAtomicParams(pmset,ATOM_HBOND_DA_STATUS);

	ResidueIteratorMolSet ritr(pmset);
	HaResidue* rptr;
	for(rptr = ritr.GetFirstRes(); rptr ; rptr = ritr.GetNextRes())
	{
		rptr->AddMissingAtoms(ADD_POLAR_HYDROGENS);
	}

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->AnnounceGeomChange();

	BoxPartition acc_part;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	pmset->GetMinMaxCrd( xmin, ymin, zmin, xmax, ymax, zmax);
	
	xmin -= 0.5;
	ymin -= 0.5;
	zmin -= 0.5;
	xmax += 0.5;
	ymax += 0.5;
	zmax += 0.5;

	acc_part.SetBoundaries(xmin, ymin, zmin, xmax, ymax, zmax);
	
	HaAtom* aptr;

	AtomIteratorMolSet aitr(pmset);

	AtomGroup acc_atoms;
	for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
	{
		if(aptr->IsHBAcceptor())
		{
			acc_atoms.InsertAtom(aptr);
		}
	}

	acc_part.DistributePointsToCells(acc_atoms);
	acc_part.SetRegionRad(max_hbond);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
	{
		if( !aptr->IsHBDonor() )
			continue;

		AtomGroup bonded_atoms;
		acc_part.GetNeighbors(*aptr, acc_atoms );
		acc_atoms.DeleteAtom(aptr);
		if(acc_atoms.empty())
			continue;

		aptr->GetBondedAtoms(bonded_atoms);
		
		HaAtom* acc_ptr;
		int nacc = acc_atoms.size();
		int iacc;
		for(iacc = 0; iacc < nacc; iacc++)
		{
			acc_ptr = acc_atoms[iacc];
			double da_dist_criteria = max_da_dist_no_acc;
			if(acc_ptr->GetElemNo() == 16)
				da_dist_criteria = max_da_dist_s_acc;

			double da_dist = Vec3D::CalcDistance(acc_ptr,aptr);
			if( da_dist > da_dist_criteria) 
				continue;
//			if( acc_ptr->GetHostMol() != aptr->GetHostMol())
//				continue;

            HaAtom* haptr;
			int nb = bonded_atoms.size();
			int ib;
			for( ib = 0; ib < nb; ib++ )
			{
				haptr = bonded_atoms[ib];
				if( !haptr->IsHydrogen())
					continue;

				double ha_dist = Vec3D::CalcDistance(acc_ptr,haptr);
				
				double ha_dist_criteria = max_ha_dist_no_acc;
				if( acc_ptr->GetElemNo() == 16 )
					ha_dist_criteria = max_ha_dist_s_acc;
				
				if( ha_dist >  ha_dist_criteria)
					continue;

				double hda_angle = Vec3D::CalcAngle(haptr,aptr,acc_ptr);

				if( hda_angle < max_hda_angle)
				{
					pmset->CreateHydrogenBond(aptr,acc_ptr,0,0);
					break;
				}
			}
		}
	}

	pmset->HBonds_found=true;

	return true;
}

int MolEditor::FindBackbone( MolSet* pmset )
{
	pmset->ClearBackbone();
	ChainIteratorMolSet chitr(pmset);
	HaChain* pch;
	for( pch = chitr.GetFirstChain(); pch ; pch = chitr.GetNextChain() )
	{
		if( pch->GetNRes() < 2 ) continue;
		ResidueIteratorChain ritr(pch);
		HaResidue* pres1 = ritr.GetFirstRes();
		HaResidue* pres2 = ritr.GetNextRes();
		for( ; pres2; pres1 = pres2, pres2 = ritr.GetNextRes() )
		{
			if( pres1->IsProtein() && pres2->IsProtein() )
			{
				HaAtom* aptr1 = pres1->GetAtomByName("CA");
				HaAtom* aptr2 = pres2->GetAtomByName("CA");
				if( aptr1 == NULL || aptr2 == NULL) continue;
				if( Vec3D::CalcDistanceSq(aptr1,aptr2) > 7.0*7.0 ) 
				{
					aptr1->flag |= BreakFlag;
					continue;
				}
				HaBond* p_bond = new HaBond( aptr1,aptr2 );
				pmset->BackboneBonds.push_back(p_bond);
			}		
			if( pres1->IsNucleo() && pres2->IsNucleo() )
			{
				HaAtom* aptr1 = pres1->GetAtomByName("P");
				HaAtom* aptr2 = pres2->GetAtomByName("P");
				if( aptr1 == NULL || aptr2 == NULL) continue;
				HaBond* p_bond = new HaBond(aptr1,aptr2);
				pmset->BackboneBonds.push_back(p_bond);
			}
		}
	}
	pmset->to_find_backb = false;
	return TRUE;
}

int MolEditor::UpdateBackBone(MolSet* pmset)
{
	if( !pmset->to_find_backb ) return TRUE;
	return FindBackbone(pmset);
}

int MolEditor::IsValidHBond(HaHBond* p_hb)
{
	AtomGroup bonded_atoms;
	HaAtom* p_don = p_hb->GetDonorAtom(); 
	p_don->GetBondedAtoms(bonded_atoms);
	HaAtom* p_acc = p_hb->GetAcceptorAtom();
	
	double da_dist_criteria = max_da_dist_no_acc;
	if(p_acc->GetElemNo() == 16) da_dist_criteria = max_da_dist_s_acc;

	double da_dist = Vec3D::CalcDistance(p_acc,p_don);
	if( da_dist > da_dist_criteria) return FALSE; 

    HaAtom* haptr;
	int nb = bonded_atoms.size();
	int ib;
	for( ib = 0; ib < nb; ib++ )
	{
		haptr = bonded_atoms[ib];
		if( !haptr->IsHydrogen())
			continue;

		double ha_dist = Vec3D::CalcDistance(p_acc,haptr);
				
		double ha_dist_criteria = max_ha_dist_no_acc;
		if( p_acc->GetElemNo() == 16 ) ha_dist_criteria = max_ha_dist_s_acc;
				
		if( ha_dist >  ha_dist_criteria) continue;

		double hda_angle = Vec3D::CalcAngle(haptr,p_don,p_acc);

		if( hda_angle < max_hda_angle) return TRUE;
	}
	return FALSE;
}

void MolEditor::CalcHydrogenBonds(MolSet* pmset)
{
    HaChain  *chn1;

    pmset->HBonds.clear();
	pmset->HBonds_found=true;

	ChainIteratorMolSet chitr(pmset);

	for(chn1= chitr.GetFirstChain(); chn1; chn1= chitr.GetNextChain())
	{
		if( !chn1->res_map.empty() )
		{   
			if( (*chn1->res_map.begin()).second->IsProtein() )
			{   
				CalcProteinHBonds(chn1);
			} 
			else if( (*chn1->res_map.begin()).second->IsDNA() )
			{
				CalcNucleicHBonds(chn1);
			}
		}
	}
	
	const int IRES_FSF = 1;
	const int IRES_CYS = 2;

	ResidueIteratorMolSet ritr1(pmset);
	ResidueIteratorMolSet ritr2(pmset);

	HaResidue* pres1;
	HaResidue* pres2;

	if(m_calc_s_hbonds_flag)
	{
		for( pres1 = ritr1.GetFirstRes(); pres1; pres1 = ritr1.GetNextRes())
		{
			int itres =0;
			if( stricmp_loc(pres1->GetName(),"FSF") == 0 )
				itres = IRES_FSF;
			else if( stricmp_loc(pres2->GetName(),"CYS") == 0 )
				itres = IRES_CYS;
			else
				continue;

			ritr2 = ritr1;
			pres2 = ritr2.GetNextRes();
			for( ; pres2; pres2 = ritr2.GetNextRes() )
			{
				if( !pres2->IsProtein() )
					continue;

				HaAtom* ptr_n= pres2->GetAtomByName("N");

				if(ptr_n == NULL)
					continue;

				HaAtom* ptr_h = pres2->GetAtomByName("H");
				if(ptr_h == NULL)
					continue;

				if(!ptr_n->IsBonded(*ptr_h))
					continue;

				HaAtom* ptr_s;
				AtomIteratorAtomGroup aitr_pres1(pres1);
				int ir_1;
				for( ptr_s = aitr_pres1.GetFirstAtom(); ptr_s; ptr_s = aitr_pres1.GetNextAtom())
				{
					if(ptr_s->GetElemNo() != 16) // atom is not Sulfur 
						continue;

					double dist= HaAtom::CalcDistance(ptr_s, ptr_n, ANGSTROM_U);
					double angle= HaAtom::CalcAngle(ptr_s,ptr_h,ptr_n);

					if(dist < 4.0 && angle > 2.0)
					{
						HaHBond* phbond= pmset->AddHBond(ptr_n,ptr_s);

						if(phbond == NULL)
						{
							ErrorInMod("HaMolecule::CalcHydrogenBonds()",
								"Failed to Add a blank HBond");
							return;
						}

						phbond->srcCA = NULL;
						phbond->dstCA = NULL;
						phbond->energy = 0;
						phbond->col = 0;
					}
				}
			}
		}
	}
     
    PrintLog("Number of H-Bonds ... %d\n",pmset->GetNHBonds());
}

static HaAtom* FindCysSulfur(HaResidue* group )
{
    const char *elem;

	AtomIteratorResidue aitr(group);
	HaAtom* aptr;

    for( aptr = aitr.GetFirstAtom();aptr; aptr = aitr.GetNextAtom() )
    {   
		elem = aptr->GetName();
        if( (elem[1]=='S') && (elem[0]==' ')  )
            return( aptr );
    }
    return( NULL );
}

static int TestDisulphideBridge(HaResidue*  group1, HaResidue* group2, HaAtom* cys1 )
{
    HaBond  *ptr;
    HaAtom  *cys2;
    double  dx, dy, dz;
    double  max,dist;

    if( !(cys2=FindCysSulfur(group2)) )
        return FALSE;

    max = 3.0*3.0;

    dx = cys1->GetX() - cys2->GetX();   if( (dist=  dx*dx) > max ) return FALSE;
    dy = cys1->GetY() - cys2->GetY();   if( (dist+= dy*dy) > max ) return FALSE;
    dz = cys1->GetZ() - cys2->GetZ();   if( (dist+= dz*dz) > max ) return FALSE;

	MolSet* pmset = cys1->GetHostMolSet();

    ptr = pmset->AddBond(cys1,cys2);
    if( !ptr ) 
	{
		ErrorInMod("TestDisulphideBridge()",
			       "Error to add SS bond");
		return FALSE;
	}

    ptr->flag = 0;
    ptr->col = 0;

    group1->flag |= CystineFlag;
    group2->flag |= CystineFlag;

	return TRUE;
}




void 
MolEditor::FindDisulphideBridges(MolSet* pmset)
{
    HaResidue  *group1;
    HaResidue  *group2;
    HaAtom  *cys;

    

	ResidueIteratorMolSet ritr1(pmset);
	ResidueIteratorMolSet ritr2(pmset);
	
	int i1=0;
	for( group1 = ritr1.GetFirstRes(); group1; group1= ritr2.GetNextRes())
	{
		if( !group1->IsCysteine() || (cys=FindCysSulfur(group1)) )
		{
			i1++; continue;
		}

		int i2= 0;
		for( group2 = ritr2.GetFirstRes(); group2; group2= ritr2.GetNextRes())
		{
			if(i2 < i1 || !group2->IsCysteine()) 
			{
				i2++; continue;
			}
			i2++;
			TestDisulphideBridge(group1,group2,cys);	
		}
		i1++;
	}

	pmset->SSBonds_found= true;
    PrintLog("Number of SS Bridges ... %d\n\n", pmset->GetNSSBonds());
}

/* Protein Donor HaAtom Coordinates */
static double hxorg,hyorg,hzorg;
static double nxorg,nyorg,nzorg;
static HaAtom  *best1CA;
static HaAtom  *best2CA;
static HaAtom  *best1;
static HaAtom  *best2;
static HaAtom  *optr;
static int res1,res2;
static int off1,off2;

/*=========================================*/
/* Kabsch & Sander Structure Determination */
/*=========================================*/


// Coupling constant for Electrostatic Energy   
// QConst = -332 * 0.42 * 0.2 * 1000.0
#define QConst (-27888.0) // Distance in Ang
#define MaxHDist (9.0*9.0)
#define MinHDist (0.5*0.5)

const double MinBondDist    = 0.16;

// Axxiliary finction to compute the energy of the hydrogen bond - consider editing:
static int CalculateHBondEnergy(HaResidue* group )
{
    double dho,dhc;
    double dnc,dno;
	
    HaAtom  *cptr;
    double dx,dy,dz;
    double dist;
    int result;
	
    if( !(cptr=group->GetAtomByName("C")) )  return(0);
    if( !(optr=group->GetAtomByName("O")) )  return(0);
	
    dx = hxorg - optr->GetX();  
    dy = hyorg - optr->GetY();  
    dz = hzorg - optr->GetZ();
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dho = sqrt(dist);
	
    dx = hxorg - cptr->GetX();  
    dy = hyorg - cptr->GetY();  
    dz = hzorg - cptr->GetZ();
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dhc = sqrt(dist);
	
    dx = nxorg - cptr->GetX();  
    dy = nyorg - cptr->GetY();  
    dz = nzorg - cptr->GetZ();
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dnc = sqrt(dist);
	
    dx = nxorg-optr->GetX();  
    dy = nyorg-optr->GetY();  
    dz = nzorg-optr->GetZ();
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dno = sqrt((double)dist);
	
    result = (int)(QConst/dho - QConst/dhc + QConst/dnc - QConst/dno);
	
    if( result<-9900 ) 
    {   
		return( -9900 );
    } 
	else if( result>-500 ) 
        return( 0 );
	
    return( result );
}

int MolEditor::CalcNucleicHBonds(HaChain* chn1 )
{
    HaResidue  *group1;
    HaResidue  *group2;
    HaResidue  *best;
    HaAtom  *ca1;
    HaAtom  *ca2;
    HaAtom  *n1;
    double max,dist;
    double dx,dy,dz;
    int refno;

	if( chn1 == NULL) return FALSE;
	HaMolecule* pmol = chn1->GetHostMol();
	if( pmol == NULL) return FALSE;
	MolSet* pmset  = pmol->GetHostMolSet();
	if( pmset == NULL) return FALSE;
	
	ResidueIteratorMolSet ritr1(pmset);
	ResidueIteratorMolSet ritr2(pmset);


    for(group1 = ritr1.GetFirstRes(); group1; group1 = ritr1.GetNextRes() )
    {   
		if( !group1->IsPurine() ) continue;
        /* Find N1 of Purine Group */
        if( !(n1=group1->GetAtomByName("N1")) )
            continue;

        /* Maximum N1-N3 distance 5A */
        refno = (group1->refno)^3; // I dont'see it's used
        max = (5.0*5.0);
        best = NULL;

		ritr2 = ritr1;
		group2 = ritr2.GetNextRes();

	    for( ; group2; group2 = ritr2.GetNextRes() ) 
        { 
			if( !group2->IsDNA()) continue;

                if( group2->refno == refno )
                {   /* Find N3 of Pyramidine Group */
                    if( !(ca1=group2->GetAtomByName("N3")) )
                        continue;

                    dx = ca1->GetX() - n1->GetX();
                    if( (dist=dx*dx) >= max ) 
                        continue;

                    dy = ca1->GetY() - n1->GetY();
                    if( (dist+=dy*dy) >= max ) 
                        continue;

                    dz = ca1->GetZ() - n1->GetZ();
                    if( (dist+=dz*dz) >= max )
                        continue;

                    best1 = ca1;
                    best = group2;
                    max = dist;
                }
        }

        if( best )
        {   
            pmset->CreateHydrogenBond( n1, best1, 0, 0 );
            if( group1->IsGuanine() )
            {   // Guanine-Cytosine 
                if( (ca1=group1->GetAtomByName("N2")) &&  /* G.N2 */
                    (ca2=best->GetAtomByName("O2")) )     /* C.O2 */
                    pmset->CreateHydrogenBond( ca1, ca2, 0, 0 );

                if( (ca1=group1->GetAtomByName("O6")) &&  /* G.O6 */
                    (ca2=best->GetAtomByName("N4")) )     /* C.N4 */
                    pmset->CreateHydrogenBond( ca1, ca2, 0, 0 );

            } 
			else // Adenine-Thymine 
                if( (ca1=group1->GetAtomByName("N6")) &&  /* A.N6 */
                    (ca2=best->GetAtomByName("O4")) )     /* T.O4 */
                    pmset->CreateHydrogenBond( ca1, ca2, 0, 0 );
        }
    }
	return TRUE;
}

int MolEditor::CalcProteinHBonds(HaChain* chn1 )
{
    int energy, offset;
    HaChain  *chn2;
    HaResidue  *group1;
    HaResidue  *group2;
	set<HaResidue, less<HaResidue> >::iterator ritr1;
	set<HaResidue, less<HaResidue> >::iterator ritr2;
    HaAtom* ca1;
    HaAtom* ca2;
    HaAtom* pc1;
    HaAtom* po1;
    HaAtom* n1;
	HaAtom* hn1;
    int pos1,pos2;
    double dx,dy,dz;
    double dco;
    double dist;

    pos1 = 0;
    pc1 = po1 = NULL;

	if( chn1 == NULL) return FALSE;
	HaMolecule* pmol = chn1->GetHostMol();
	if( pmol == NULL) return FALSE;
	MolSet* pmset  = pmol->GetHostMolSet();
	if( pmset == NULL) return FALSE;

	ResidueIteratorChain ritr_ch_1(chn1);
    for( group1 = ritr_ch_1.GetFirstRes(); group1 ; group1 = ritr_ch_1.GetNextRes() )
    {   
		pos1++;
		if( pc1 && po1 )
		{   
			dx = pc1->GetX() - po1->GetX();
			dy = pc1->GetY() - po1->GetY();
			dz = pc1->GetZ() - po1->GetZ();
		} 
		else
		{   
			pc1 = group1->GetAtomByName("C");
            po1 = group1->GetAtomByName("O");
            continue;
        }

        pc1 = group1->GetAtomByName("C");
        po1 = group1->GetAtomByName("O");

        if( !group1->IsAmino() || group1->IsProline() )
            continue;

        if( !(ca1=group1->GetAtomByName("CA")) ) continue;
        if( !(n1=group1->GetAtomByName("N")) )  continue;

        dist = dx*dx + dy*dy + dz*dz;
        dco = sqrt( dist );

        nxorg = n1->GetX();   
        nyorg = n1->GetY();   
        nzorg = n1->GetZ();  

		if( (hn1= group1->GetAtomByName("HN")) != NULL )
		{
			hxorg = n1->GetX();   
			hyorg = n1->GetY();   
			hzorg = n1->GetZ();  
		}
		else
		{
			hxorg = nxorg + dx/dco;
			hyorg = nyorg + dy/dco;
			hzorg = nzorg + dz/dco;
		}
		
        res1 = res2 = 0;

        /* Only Hydrogen Bond within a single chain!       */

        chn2 = chn1;
        {   /* Only consider non-empty peptide chains! */

            pos2 = 0;
//			cout << " res # " << group1->serno << endl;
			ResidueIteratorChain ritr_ch_2(chn2);
			for( group2 = ritr_ch_2.GetFirstRes(); group2; group2 = ritr_ch_2.GetNextRes() )
            {   
				pos2++;
                if( (group2==group1) || (group2->GetNextResInChain() == group1) )
                    continue;

                if( !group2->IsAmino() ) 
                    continue;
                if( !(ca2=group2->GetAtomByName("CA")) ) 
                    continue;

                dx = ca1->GetX() - ca2->GetX();
                if( (dist= dx*dx) > MaxHDist )
                    continue;

                dy = ca1->GetY() - ca2->GetY();
                if( (dist+= dy*dy) > MaxHDist )
                    continue;

                dz = ca1->GetZ() - ca2->GetZ();
                if( (dist+= dz*dz) > MaxHDist )
                    continue;

                if( (energy = CalculateHBondEnergy(group2)) )
                {   
					if( chn1 == chn2 )
                    {   
						offset = pos1 - pos2;
                    } 
					else 
						offset = 0;

                    if( energy<res1 )
                    {   
						best2CA = best1CA;  best1CA = ca2;
                        best2 = best1;      best1 = optr;
                        res2 = res1;        res1 = energy;
                        off2 = off1;        off1 = offset;
                    } 
					else if( energy<res2 )
                    {   
						best2CA = ca2;
                        best2 = optr;
                        res2 = energy;
                        off2 = offset;
                    }
                }
            }  /* group2 */
        }      /* chn2 */

        if( res1 ) 
        {   
			if( res2 ) 
				pmset->CreateHydrogenBond(n1,best2,res2,off2);
            pmset->CreateHydrogenBond(n1,best1,res1,off1);
        }
    }
	return TRUE;
}


vector<HaAtom> MolEditor::FindMissingHydrogens( HaAtom* aptr)
{
	vector<HaAtom> extra_h;
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	

	HaAtom* aptr_templ = p_res_db->GetTemplateForAtom(aptr);
	if( aptr_templ == NULL) return extra_h;

	AtomGroup bonded_atoms, bonded_atoms_templ;
	aptr->GetBondedAtoms( bonded_atoms );
	aptr_templ->GetBondedAtoms( bonded_atoms_templ );

	vector<int> from_templ_atmap( bonded_atoms_templ.size(), -1);
	vector<int> to_templ_atmap  ( bonded_atoms.size(), -1);

	int nh = 0;
	int nh_t = 0;

	int i,j;
	for( i = 0; i < bonded_atoms.size(); i++)
	{
		if( bonded_atoms[i]->GetElemNo() == 1) nh++;
	}

	for( j = 0; j < bonded_atoms_templ.size(); j++)
	{
		if( bonded_atoms_templ[j]->GetElemNo() == 1) nh_t++;
	}
	
	if( nh >= nh_t) return extra_h;

// build maps of bonded atoms for an atom and its template

	for( j = 0; j < bonded_atoms_templ.size(); j++)
	{
		if( bonded_atoms_templ[j]->GetHostRes() != aptr_templ->GetHostRes()) 
		{
			for( i = 0; i < bonded_atoms.size(); i++)
			{
				if( bonded_atoms[i]->GetHostRes() != aptr->GetHostRes())
				{
					from_templ_atmap[j] = i;
					to_templ_atmap[i] = j;
				}
			}
		}
		else
		{
			std::string at_name = bonded_atoms_templ[j]->GetName();
			for( i = 0; i < bonded_atoms.size(); i++)
			{
				if(!stricmp_trunc( bonded_atoms[i]->GetName(), at_name.c_str()) )
				{
					from_templ_atmap[j] = i;
					to_templ_atmap[i] = j;
				}
			}
		}
	}

	for( j = 0; j < bonded_atoms_templ.size(); j++)
	{
		if( bonded_atoms_templ[j]->GetElemNo() != 1) continue;
		if( from_templ_atmap[j] >= 0) continue;
	
		HaAtom* hptr_t = bonded_atoms_templ[j];

		HaAtom* aptr1 = NULL;
		HaAtom* aptr2 = NULL;
		HaAtom* aptr1_t = NULL;
		HaAtom* aptr2_t = NULL;

		for(i = 0; i < bonded_atoms.size(); i++)
		{
		   if( aptr2 != NULL) break;
		   if( to_templ_atmap[i] >= 0 ) 
		   {
			   if( aptr1 == NULL )
			   {
					aptr1 = bonded_atoms[i];
					aptr1_t = bonded_atoms_templ[ to_templ_atmap[i] ];
			   }
			   else
			   {
					aptr2 = bonded_atoms[i];
					aptr2_t = bonded_atoms_templ[ to_templ_atmap[i] ];
			   }
		   }
		}
		if( aptr1 == NULL) continue; // Can't place hydrogen atoms
		HaAtom h_atom;
		double dist      = Vec3D::CalcDistance(hptr_t,aptr_templ);
		double val_angle = Vec3D::CalcAngle(hptr_t,aptr_templ,aptr1);
		double torsion_angle;
		if( aptr2 != NULL) 
		{	
			torsion_angle = Vec3D::CalcTorsion(hptr_t,aptr_templ,aptr1,aptr2);
		}
		else
		{
			torsion_angle = 0.0;
		}
		Vec3D::SetAtomPos(&h_atom, aptr, aptr1, aptr2, dist, val_angle, torsion_angle );
		h_atom.SetHostRes( aptr->GetHostRes());
		h_atom.SetParamFrom( *aptr_templ);
		extra_h.push_back(h_atom);
		bonded_atoms.push_back( &(extra_h.back()) );
	}
	return extra_h;
}

int MolEditor::SetBondDist(HaAtom* aptr1, HaAtom* aptr2, double new_dist)
//!  Atoms moved are aptr2 and all atoms reachable from aptr2 
//!  by a graph of covalent bonds
{
	if(aptr1 == NULL || aptr2 == NULL) 
	{
		PrintLog(" Error in MolSet::SetBondDist() \n");
		PrintLog(" atom pointer is zero \n");
		return FALSE;
	}

	double old_dist = Vec3D::CalcDistance(aptr1,aptr2);
	double diff = new_dist - old_dist;

	if(fabs(diff) < DBL_EPSILON) return TRUE;

	Vec3D dn; 
    dn[0] = aptr2->GetX() - aptr1->GetX();
    dn[1] = aptr2->GetY() - aptr1->GetY();
    dn[2] = aptr2->GetZ() - aptr1->GetZ();

	double len = dn.length();
	if(len < DBL_EPSILON) 
	{
		PrintLog(" Error in MolSet::SetBondDist() \n");
		PrintLog(" coordinates of atoms of the bond to change coincide \n");
		return FALSE;
	}

	dn.normalize();

	AtomGroup mov_atoms;
	AtomGroup block_atoms;
	int loop;

	block_atoms.InsertAtom(aptr1);
	int ires = HaAtom::GetReachableAtoms(block_atoms, aptr2, mov_atoms,loop);

	if(loop) 
	{
		PrintLog(" Error in MolSet::SetBondDist() \n");
		PrintLog(" atoms are not moved because of the loop \n");
		return FALSE;
	}

	AtomIteratorAtomGroup aitr(&mov_atoms);

	HaAtom* aptr = aitr.GetFirstAtom();
	for(; aptr; aptr = aitr.GetNextAtom())
	{
		aptr->SetX( aptr->GetX() + diff * dn.GetX());
		aptr->SetY( aptr->GetY() + diff * dn.GetY());
		aptr->SetZ( aptr->GetZ() + diff * dn.GetZ());
	}
	return TRUE;
}

int MolEditor::SetAngle(HaAtom* aptr1, HaAtom* aptr2, HaAtom* aptr3, double ang_new)
//!  Atoms moved are aptr3 and all atoms reachable from aptr3 
//!  by a graph of covalent bonds
//!  Also moved atoms conected to aptr2 by 1/2 of ang_new
{
	if(aptr1 == NULL || aptr2 == NULL || aptr3 == NULL) 
	{
		PrintLog(" Error in MolSet::SetAng() \n");
		PrintLog(" One of atom pointers is NULL \n");
		return FALSE;
	}

	double old_ang = Vec3D::CalcAngle(aptr1,aptr2,aptr3);
	double diff = ang_new - old_ang;

    if(fabs(diff) < DBL_EPSILON) return TRUE;

	int ires;
	double len;
	Vec3D dn1,dn2,dn3;
	
    dn1[0] = aptr1->GetX() - aptr2->GetX();
    dn1[1] = aptr1->GetY() - aptr2->GetY();
    dn1[2] = aptr1->GetZ() - aptr2->GetZ();

    dn2[0] = aptr3->GetX() - aptr2->GetX();
    dn2[1] = aptr3->GetY() - aptr2->GetY();
    dn2[2] = aptr3->GetZ() - aptr2->GetZ();

	if(dn1.length2() < DBL_EPSILON || dn2.length2() < DBL_EPSILON ) 
	{
		PrintLog(" Error in MolSet::SetAng() \n");
		PrintLog(" Coordinates of atoms of the angle to change coincide \n");
		return FALSE;
	}

	Vec3D::VecProduct(dn3,dn1,dn2);

	len = dn3.length();
	if(len < DBL_EPSILON) // bonds are collinear
	{
		if( fabs(dn1[0]) < DBL_EPSILON) 
			dn2[0] = 1.0;
		else
			dn2[0] = 0.0;

		if( fabs(dn1[1]) < DBL_EPSILON) 
			dn2[1] = 1.0;
		else
			dn2[1] = 0.0;

		Vec3D::VecProduct(dn3,dn1,dn2);
	}

	dn3.normalize();

	AtomGroup mov_atoms;
	AtomGroup block_atoms;
	int loop;

	block_atoms.InsertAtom(aptr2);
	ires = HaAtom::GetReachableAtoms(block_atoms, aptr3, mov_atoms,loop);

	if(loop) 
	{
		PrintLog(" Error in MolSet::SetAng() \n");
		PrintLog(" atoms are not moved because of the loop \n");
		return FALSE;
	}

	AtomIteratorAtomGroup aitr(&mov_atoms);
	double cosa = cos(diff);
	double sina = sin(diff);

	HaAtom* aptr = aitr.GetFirstAtom();
	for(; aptr; aptr = aitr.GetNextAtom())
	{
		aptr->RotatePt(*aptr2,dn3,cosa,sina);
	}
	return TRUE;
}

int MolEditor::SetTorsion(HaAtom* aptr1, HaAtom* aptr2, HaAtom* aptr3, HaAtom* aptr4, double tors_new)
//!  Atoms moved are aptr4 and all atoms reachable from aptr4 
//!  by a graph of covalent bonds
{
	if(aptr1 == NULL || aptr2 == NULL || aptr3 == NULL || aptr4 == NULL) 
	{
		PrintLog(" Error in MolSet::SetTorsion() \n");
		PrintLog(" One of atom pointers is NULL \n");
		return FALSE;
	}

	double old_tors = Vec3D::CalcTorsion(aptr1,aptr2,aptr3,aptr4);
	double diff = tors_new - old_tors;

    if(fabs(diff) < DBL_EPSILON) return TRUE;

	int ires;
	Vec3D dn;
	
    dn[0] = aptr3->GetX() - aptr2->GetX();
    dn[1] = aptr3->GetY() - aptr2->GetY();
    dn[2] = aptr3->GetZ() - aptr2->GetZ();


	if(dn.length2() < DBL_EPSILON  ) 
	{
		PrintLog(" Error in MolSet::SetAng() \n");
		PrintLog(" Coordinates of atoms of the torsion to change coincide \n");
		return FALSE;
	}

	dn.normalize();

	AtomGroup mov_atoms;
	AtomGroup block_atoms;
	int loop;

	block_atoms.InsertAtom(aptr1);
	block_atoms.InsertAtom(aptr2);

	ires = HaAtom::GetReachableAtoms(block_atoms, aptr3, mov_atoms,loop);

	if(loop) 
	{
		PrintLog(" Error in MolSet::SetTorsion() \n");
		PrintLog(" atoms are not moved because of the loop \n");
		return FALSE;
	}

	AtomIteratorAtomGroup aitr(&mov_atoms);
	double cosa = cos(diff);
	double sina = sin(diff);

	HaAtom* aptr = aitr.GetFirstAtom();
	for(; aptr; aptr = aitr.GetNextAtom())
	{
		aptr->RotatePt(*aptr2,dn,cosa,sina);
	}
	return TRUE;
}


int MolEditor::DeleteExtraAtoms( MolSet* pmset)
{
	AtomGroup extra_atoms;
	ResidueIteratorMolSet ritr(pmset);
	HaAtom* aptr;
	HaAtom* aptr_templ;
	
	HaResidue* pres; 
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		if(!pres->HasSelectedAtoms()) continue;
		HaResidue* res_templ = pres->GetTemplate();
//		AtomIteratorAtomGroup aitr_t(res_templ);
//		for( aptr_templ = aitr_t.GetFirstAtom(); aptr_templ; aptr_templ = aitr_t.GetNextAtom())
//		{
//			PrintLog("MolEditor::DeleteExtraAtoms() res_templ_name %s =  at_templ_name=%s \n",
//				      res_templ->GetName(), aptr_templ->GetName() );
//		}
		if(res_templ == NULL) continue;
		AtomIteratorAtomGroup aitr(pres);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			std::string atname = aptr->GetName();
			aptr_templ = res_templ->GetAtomByName(atname.c_str());
			if( aptr_templ == NULL)
			{
				extra_atoms.InsertAtom(aptr);
			}
		}
	}
	char buf[256];
	AtomIteratorAtomGroup aitr(&extra_atoms);
	PrintLog("Delete %d Atoms not found in residue templates \n",extra_atoms.size()); 

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->FillRef(buf);
		PrintLog("%s \n",buf);
	}
	pmset->DeleteAtoms(extra_atoms);
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );
	return TRUE;
}

int MolEditor::AddMissingAtoms(MolSet* pmset ) 
{
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		if(pres->size() == 0) continue;
		HaAtom* aptr = (*pres)[0];
		if(aptr->Selected()) pres->AddMissingAtoms(ADD_ALL_ATOMS);
	}
	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );	
	return TRUE;
}

int MolEditor::FixBondsUsingTempl(MolSet* pmset)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		if( pres->size() == 0) continue;
		HaAtom* aptr = (*pres)[0];
		if(aptr == NULL) continue;
		if(!aptr->Selected()) continue;

		HaResidue* pres_next = pres->GetNextResInChain();
		HaResidue* pres_prev = pres->GetPrevResInChain();
		
		HaChain* chain = pres->GetHostChain();

		HaResidue* res_templ = pres->GetTemplate();
		if(res_templ == NULL) continue;
		if(pres->GetNAtoms() < 2) continue;

		if( pres_next != NULL && pres_next->GetHostChain() != chain ) pres_next = NULL;
		if( pres_prev != NULL && pres_prev->GetHostChain() != chain ) pres_prev = NULL;

		PtrPtrMap at_map;
		int na = pres->GetNAtoms();
		int i,j;
		HaAtom* aptr2;
		HaAtom* atempl;
		HaAtom* atempl2;

		for(i = 0; i < na; i++)
		{
			aptr = (*pres)[i];
			atempl = p_res_db->GetTemplateForAtom(aptr);
			if( atempl == NULL) 
			{
				PrintLog("Warning in MolEditor::FixBondsUsingTempl() \n");
				PrintLog("Atom Does not have a match in the template: %s \n",(aptr->GetRef()).c_str());
				PrintLog("Remove Extra Atoms \n");
				at_map[aptr] = NULL;
			}
			else
			{
				at_map[aptr] = atempl;
			}
		}

		for( i = 0; i < na; i++)
		{
			aptr = (*pres)[i];
			atempl = (HaAtom*) at_map[aptr];
			if( atempl == NULL) continue;

			for( j = 0; j < i; j++)
			{
				aptr2 = (*pres)[j];
				atempl2 = (HaAtom*) at_map[aptr2];
				if( atempl2 == NULL) continue;
				
				bool bonded_r = aptr->IsBonded(*aptr2);
				bool bonded_t = atempl->IsBonded(*atempl2);
				
				if( !bonded_r && !bonded_t ) continue;
				if( bonded_r && !bonded_t )
				{
					PrintLog(" Delete Bond between atoms %s %s \n", aptr->GetRef().c_str(), aptr2->GetRef().c_str() );
					HaAtom::DeleteBond(aptr,aptr2);
					continue;
				}
				if( !bonded_r && bonded_t )
				{
					PrintLog(" Create Bond between atoms %s %s \n", aptr->GetRef().c_str(), aptr2->GetRef().c_str() );
					HaAtom::CreateBond(aptr,aptr2);
				}
			}
			std::vector<HaBond*>& bonds_r = aptr->GetBonds();
			std::vector<HaBond*>& bonds_t = atempl->GetBonds();
			int k;
			for( k = 0; k < bonds_r.size(); k++)
			{
				HaBond* pbnd_r = bonds_r[k];
				HaAtom* aptr2 = pbnd_r->srcatom;
				if( aptr2 == aptr ) aptr2 = pbnd_r->dstatom;
				atempl2 = (HaAtom*) at_map[aptr2];
				if( atempl2 == NULL) continue;
				HaBond* pbnd_t = NULL;
				for( j = 0; j < bonds_t.size(); j++)
				{
					if( bonds_t[j]->srcatom == atempl && bonds_t[j]->dstatom == atempl2 )
					{
						pbnd_t = bonds_t[j]; 
						break;
					}
					if( bonds_t[j]->dstatom == atempl && bonds_t[j]->srcatom == atempl2 )
					{
						pbnd_t = bonds_t[j]; 
						break;
					}
				}
				if( pbnd_t ) 
				{
					if( pbnd_t->GetTypeString() != pbnd_r->GetTypeString() )
					{
						PrintLog(" Change Bond type between atoms %s %s  from %s to %s \n", aptr->GetRef().c_str(), aptr2->GetRef().c_str(), pbnd_r->GetTypeString().c_str(), pbnd_t->GetTypeString().c_str() );
						pbnd_r->SetTypeFrom( *pbnd_t);
					}
				}
			}
			AtomGroup bonded_atoms;
			aptr->GetBondedAtoms(bonded_atoms);
			std::string at_name = aptr->GetName();
			AtomIteratorAtomGroup aitr_b(&bonded_atoms);
			HaAtom* aptr_b;
			for( aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom() )
			{
				if( aptr_b->GetHostRes() == pres ) continue;

				std::string at_name_b = aptr_b->GetName();
				
				if( aptr_b->GetHostRes() == pres_next )
				{
					if( at_name == "C" && at_name_b == "N" ) continue;
				}

				if( aptr_b->GetHostRes() == pres_prev )
				{
					if( at_name == "N" && at_name_b == "C" ) continue;
				}
				if( at_name == "SG" && at_name_b == "SG" ) continue;

				PrintLog(" Delete Bond between atoms %s %s \n", aptr->GetRef().c_str(), aptr_b->GetRef().c_str() );
				HaAtom::DeleteBond(aptr,aptr_b);
			}
		}
	}
	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );	
	return TRUE;
}

int MolEditor::OrderAtomsInRes(MolSet* pmset)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		try
		{
			if(pres->size() == 0) continue;
			HaAtom* aptr = (*pres)[0];
			if(!aptr->Selected()) continue;
			HaResidue* res_templ = pres->GetTemplate();
			if(res_templ == NULL) continue;

			int na = res_templ->GetNAtoms();
			int n_nproxy = res_templ->GetNAtomsNonProxy();
			
			if( pres->GetNAtoms() != n_nproxy ) throw std::runtime_error( "Number of atoms in Residue " + pres->GetRef() + " is not equal to that of the template ");
			
			int i,j;
			for(i = 0; i < na; i++)
			{
				HaAtom* aptr_t = (*res_templ)[i];
				if( aptr_t->IsProxy() ) continue;

				HaAtom* aptr_r = (*pres)[i];
				
				std::string atn_t = aptr_t->GetName(); 

				if( stricmp_loc(aptr_r->GetName(), atn_t.c_str() ) == 0) continue;

				aptr_r = pres->GetAtomByName( atn_t.c_str() );
				if( aptr_r == NULL)  throw std::runtime_error( "Can not find atom " + atn_t + " in the residue " + pres->GetRef() + " found in its template ");
				
				int atom_found = FALSE;
				for(j = na-1; j >= i; j--)
				{
					HaAtom* aptr2 = (*pres)[j];
					if( atom_found ) (*pres)[j+1] = aptr2;
					if( aptr2 == aptr_r) atom_found = TRUE;
				}
				(*pres)[i] = aptr_r;
			}
		}
		catch(std::exception& ex)
		{
			PrintLog("Error in MolEditor::OrderAtomsInRes() \n");
			PrintLog("%s\n",ex.what());
		}
	}
	return TRUE;
}

int MolEditor::RenameAtomsToAmber(MolSet* pmset)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		try
		{
			if( pres->size() == 0) continue;
			HaAtom* aptr = (*pres)[0];
			if(!aptr->Selected()) continue;
			if( !pres->IsAmino() ) continue;
			std::string res_name = pres->GetName();
			std::string name_mod = pres->GetNameModifier();
			HaAtom* phb1 = pres->GetAtomByName("HB1");
			HaAtom* phb2 = pres->GetAtomByName("HB2");
			HaAtom* phb3 = pres->GetAtomByName("HB3");
			if( phb1 && phb2 && !phb3 )
			{
				phb2->SetName("HB3");
				phb1->SetName("HB2");
				PrintLog("In Residue %s  Rename HB2->HB3, HB1->HB2 \n",pres->GetRef().c_str());
			}
			if( res_name == "ILE")
			{
				HaAtom* pcd = pres->GetAtomByName("CD");
				HaAtom* phd1 = pres->GetAtomByName("HD1");
				HaAtom* phd2 = pres->GetAtomByName("HD2");
				HaAtom* phd3 = pres->GetAtomByName("HD3");

				HaAtom* phg11 = pres->GetAtomByName("HG11");
				HaAtom* phg12 = pres->GetAtomByName("HG12");
				if( pcd ) 
				{
					pcd->SetName("CD1");
					PrintLog("In Residue %s  Rename CD->CD1 \n",pres->GetRef().c_str());
				}
				if( phd1 ) 
				{
					phd1->SetName("HD11");
					PrintLog("In Residue %s  Rename HD1->HD11 \n",pres->GetRef().c_str());
				}
				if( phd2 ) 
				{
					phd2->SetName("HD12");
					PrintLog("In Residue %s  Rename HD2->HD21 \n",pres->GetRef().c_str());
				}
				if( phd3 ) 
				{
					phd3->SetName("HD13");
					PrintLog("In Residue %s  Rename HD3->HD31 \n",pres->GetRef().c_str());
				}
				if( phg12 ) 
				{
					phg12->SetName("HG13");
					PrintLog("In Residue %s  Rename HG12->HG13 \n",pres->GetRef().c_str());
				}
				if( phg11 ) 
				{
					phg11->SetName("HG12");
					PrintLog("In Residue %s  Rename HG11->HG12 \n",pres->GetRef().c_str());
				}
			}
			
			if( res_name == "LYS" || res_name == "PRO" || res_name == "ARG")
			{
				HaAtom* phd1 = pres->GetAtomByName("HD1");
				HaAtom* phd2 = pres->GetAtomByName("HD2");
				if( phd2 ) 
				{
					phd2->SetName("HD3");
					PrintLog("In Residue %s  Rename HD2->HD3 \n",pres->GetRef().c_str());
				}
				if( phd1 ) 
				{
					phd1->SetName("HD2");
					PrintLog("In Residue %s  Rename HD1->HD2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "GLY" )
			{
				HaAtom* pha1 = pres->GetAtomByName("HA1");
				HaAtom* pha2 = pres->GetAtomByName("HA2");
				
				if( pha2 ) 
				{
					pha2->SetName("HA3");
					PrintLog("In Residue %s  Rename HA2->HA3 \n",pres->GetRef().c_str());
				}
				if( pha1 ) 
				{
					pha1->SetName("HA2");
					PrintLog("In Residue %s  Rename HA1->HA2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "LYS" )
			{
				HaAtom* phe1 = pres->GetAtomByName("HE1");
				HaAtom* phe2 = pres->GetAtomByName("HE2");
				
				if( phe2 ) 
				{
					phe2->SetName("HE3");
					PrintLog("In Residue %s  Rename HE2->HE3 \n",pres->GetRef().c_str());
				}
				if( phe1 ) 
				{
					phe1->SetName("HE2");
					PrintLog("In Residue %s  Rename HE1->HE2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "GLN" || res_name == "GLU" || res_name == "LYS" || res_name == "PRO" || res_name == "MET" || res_name == "ARG")
			{
				HaAtom* phg1 = pres->GetAtomByName("HG1");
				HaAtom* phg2 = pres->GetAtomByName("HG2");
				if( phg2 ) 
				{
					phg2->SetName("HG3");
					PrintLog("In Residue %s  Rename HG2->HG3 \n",pres->GetRef().c_str());
				}
				if( phg1 ) 
				{
					phg1->SetName("HG2");
					PrintLog("In Residue %s  Rename HG1->HG2 \n",pres->GetRef().c_str());
				}
			}

			if( name_mod == "CT")
			{
				HaAtom* poc1 = pres->GetAtomByName("OC1");
				HaAtom* poc2 = pres->GetAtomByName("OC2");
				if( poc1 ) 
				{
					poc1->SetName("O");
					PrintLog("In Residue %s  Rename OC1->OC \n",pres->GetRef().c_str());
				}
				if( poc2 ) 
				{
					poc2->SetName("OXT");
					PrintLog("In Residue %s  Rename OC2->OXT \n",pres->GetRef().c_str());
				}
			}



		}
		catch(std::exception& ex)
		{
			PrintLog("Error in MolEditor::RenameAminoAtomsToAmber() \n");
			PrintLog("%s\n",ex.what());
		}
	}
	return TRUE;
}

int MolEditor::RenameAtomsToGromacs(MolSet* pmset)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		try
		{
			if( pres->size() == 0) continue;
			HaAtom* aptr = (*pres)[0];
			if(!aptr->Selected()) continue;
			if( !pres->IsAmino() ) continue;
			std::string res_name = pres->GetName();
			std::string name_mod = pres->GetNameModifier();
			HaAtom* phb1 = pres->GetAtomByName("HB1");
			HaAtom* phb2 = pres->GetAtomByName("HB2");
			HaAtom* phb3 = pres->GetAtomByName("HB3");
			if( phb2 && phb3 && !phb1 )
			{
				phb2->SetName("HB1");
				phb3->SetName("HB2");
				PrintLog("In Residue %s  Rename HB2->HB1, HB3->HB2 \n",pres->GetRef().c_str());
			}
			if( res_name == "ILE")
			{
				HaAtom* pcd1 = pres->GetAtomByName("CD1");
				HaAtom* phd11 = pres->GetAtomByName("HD11");
				HaAtom* phd12 = pres->GetAtomByName("HD12");
				HaAtom* phd13 = pres->GetAtomByName("HD13");

				HaAtom* phg12 = pres->GetAtomByName("HG12");
				HaAtom* phg13 = pres->GetAtomByName("HG13");
				if( pcd1 ) 
				{
					pcd1->SetName("CD");
					PrintLog("In Residue %s  Rename CD1->CD \n",pres->GetRef().c_str());
				}
				if( phd11 ) 
				{
					phd11->SetName("HD1");
					PrintLog("In Residue %s  Rename HD11->HD1 \n",pres->GetRef().c_str());
				}
				if( phd12 ) 
				{
					phd12->SetName("HD2");
					PrintLog("In Residue %s  Rename HD21->HD2 \n",pres->GetRef().c_str());
				}
				if( phd13 ) 
				{
					phd13->SetName("HD3");
					PrintLog("In Residue %s  Rename HD31->HD3 \n",pres->GetRef().c_str());
				}
				if( phg12 ) 
				{
					phg12->SetName("HG11");
					PrintLog("In Residue %s  Rename HG12->HG11 \n",pres->GetRef().c_str());
				}
				if( phg13 ) 
				{
					phg13->SetName("HG12");
					PrintLog("In Residue %s  Rename HG13->HG12 \n",pres->GetRef().c_str());
				}
			}
			
			if( res_name == "LYS" || res_name == "PRO" || res_name == "ARG")
			{
				HaAtom* phd2 = pres->GetAtomByName("HD2");
				HaAtom* phd3 = pres->GetAtomByName("HD3");
				if( phd2 ) 
				{
					phd2->SetName("HD1");
					PrintLog("In Residue %s  Rename HD2->HD1 \n",pres->GetRef().c_str());
				}
				if( phd3 ) 
				{
					phd3->SetName("HD2");
					PrintLog("In Residue %s  Rename HD3->HD2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "GLY" )
			{
				HaAtom* pha2 = pres->GetAtomByName("HA2");
				HaAtom* pha3 = pres->GetAtomByName("HA3");
				
				if( pha2 ) 
				{
					pha2->SetName("HA1");
					PrintLog("In Residue %s  Rename HA2->HA1 \n",pres->GetRef().c_str());
				}
				if( pha3 ) 
				{
					pha3->SetName("HA2");
					PrintLog("In Residue %s  Rename HA3->HA2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "LYS" )
			{
				HaAtom* phe2 = pres->GetAtomByName("HE2");
				HaAtom* phe3 = pres->GetAtomByName("HE3");
				
				if( phe2 ) 
				{
					phe2->SetName("HE1");
					PrintLog("In Residue %s  Rename HE2->HE1 \n",pres->GetRef().c_str());
				}
				if( phe3 ) 
				{
					phe3->SetName("HE2");
					PrintLog("In Residue %s  Rename HE3->HE2 \n",pres->GetRef().c_str());
				}
			}

			if( res_name == "GLN" || res_name == "GLU" || res_name == "LYS" || res_name == "PRO" || res_name == "MET" || res_name == "ARG")
			{
				HaAtom* phg2 = pres->GetAtomByName("HG2");
				HaAtom* phg3 = pres->GetAtomByName("HG3");
				if( phg2 ) 
				{
					phg2->SetName("HG1");
					PrintLog("In Residue %s  Rename HG2->HG1 \n",pres->GetRef().c_str());
				}
				if( phg3 ) 
				{
					phg3->SetName("HG2");
					PrintLog("In Residue %s  Rename HG3->HG2 \n",pres->GetRef().c_str());
				}
			}

			if( name_mod == "CT")
			{
				HaAtom* po = pres->GetAtomByName("O");
				HaAtom* poxt = pres->GetAtomByName("OXT");
				if( po ) 
				{
					po->SetName("OC1");
					PrintLog("In Residue %s  Rename OC1->O \n",pres->GetRef().c_str());
				}
				if( poxt ) 
				{
					poxt->SetName("OC2");
					PrintLog("In Residue %s  Rename OXT->OC2 \n",pres->GetRef().c_str());
				}
			}



		}
		catch(std::exception& ex)
		{
			PrintLog("Error in MolEditor::RenameAminoAtomsToAmber() \n");
			PrintLog("%s\n",ex.what());
		}
	}
	return TRUE;
}

int MolEditor::ConvertWaterArrowVB(MolSet* pmset)
{
	ResidueIteratorMolSet ritr(pmset);
	
	HaResidue* pres;
	for (pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		bool is_water = false;
		HaMolecule* pMol = pres->GetHostMol();
		std::string mol_name = pMol->GetObjName();
		int nres = pMol->GetNRes();
		if (nres == 1 && mol_name == "HOH" || nres == 1 && mol_name == "WAT")
			is_water = true;
		if (pres->IsWater()) is_water = true;
		if (!is_water) continue;

		if ( pres->GetNAtoms() != 3 )
		{
			PrintLog("Warning in MolEditor::ConvertWaterArrowVB() : ");
			PrintLog("Water Residue %s doesn't have 3 atoms - skip", pres->GetRef().c_str() );
			continue;
		}
		HaAtom* pox = nullptr;
		HaAtom* ph1 = nullptr;
		HaAtom* ph2 = nullptr;

		AtomIteratorResidue aitr(pres);
		HaAtom* aptr;
		for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if (aptr->GetElemNo() == 8) pox = aptr;
			if (aptr->GetElemNo() == 1)
			{
				if (ph1 == nullptr)
					ph1 = aptr;
				else
					ph2 = aptr;
			}
		}	
		if (pox == nullptr || ph1 == nullptr || ph2 == nullptr )
		{
			PrintLog("Warning in MolEditor::ConvertWaterArrowVB() : ");
			PrintLog("Can not assign O,H1,H2 atoms for Water Residue %s - skip", pres->GetRef().c_str());
			continue;
		}
		pox->SetName("O");
		ph1->SetName("H1");
		ph2->SetName("H2");
		
		pox->SetFFSymbol("OW");
		ph1->SetFFSymbol("HW");
		ph2->SetFFSymbol("HW");

		if(!pox->IsBonded(*ph1)) HaAtom::CreateBond(pox, ph1);
		if(!pox->IsBonded(*ph2)) HaAtom::CreateBond(pox, ph2);
		if(!ph1->IsBonded(*ph2)) HaAtom::CreateBond(ph1, ph2);

		HaAtom::BondIterator bitr = ph1->Bonds_begin();
		for (; bitr != ph1->Bonds_end(); bitr++)
		{
			HaBond* bptr = *bitr;
			HaAtom* at2 = nullptr;
			if (bptr->GetFirstAtom() == ph1) 
				at2 = bptr->GetSecondAtom();
			else
				at2 = bptr->GetFirstAtom();
			if (at2 == ph2) bptr->SetVirtual();
		}
	//	pres = wat_res[i];
	//	HaMolecule* old_mol = pres->GetHostMol();
	//	int nr_mol = old_mol->GetNRes();
	//	if (nr_mol == 1) continue;
	//	HaMolecule* pmol = pmset->AddNewMolecule();
	//	pmol->SetObjName("HOH");
	//	HaChain* pchain = pmol->AddChain(' ');
	//	HaResidue* pres2 = pchain->AddResidue(1);
	}


	PrintLog("MolEditor::ConvertWaterArrowVB() \n");
	return TRUE;
}

int MolEditor::AddHydrogens(MolSet* pmset)
{
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		if( pres->size() == 0) continue;
		HaAtom* aptr = (*pres)[0];
		if(aptr->Selected()) pres->AddMissingAtoms(ADD_HYDROGENS);
	}

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );
	return TRUE;
}

int MolEditor::AddPolarHydrogens(MolSet* pmset) 
{
	if( pmset == NULL) return FALSE;
	SetStdAtomicParams(pmset,ATOM_HBOND_DA_STATUS);
	ResidueIteratorMolSet ritr(pmset);
	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres ; pres = ritr.GetNextRes())
	{
		if(pres->size() == 0) continue;
		HaAtom* aptr = (*pres)[0];
		if(aptr->Selected()) pres->AddMissingAtoms(ADD_POLAR_HYDROGENS);
	}

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );	
	return TRUE;
}

int MolEditor::AddHydrogensHybrid(MolSet* pmset)
{
    AtomGroup selected_nonh_atoms;

	double std_h_sp3_dist = 1.08;
	double std_h_sp2_dist = 1.05;
	double std_h_sp_dist  = 1.00;

	double dist, val_angle, dih_angle;

    HaAtom* aptr;
	HaBond* new_bond;

    AtomIteratorMolSet aitr(pmset);

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( !aptr->IsHydrogen() && aptr->Selected())
		{
           selected_nonh_atoms.InsertAtom(aptr);
		}
	}

	AtomGroup bonded_atoms;
	AtomIteratorAtomGroup aitr_sel_nonh_atoms(&selected_nonh_atoms);

	for(aptr = aitr_sel_nonh_atoms.GetFirstAtom(); aptr; aptr = aitr_sel_nonh_atoms.GetNextAtom() )
	{
		PrintLog(" atom = %s %d\n", aptr->GetRef().c_str(), aptr->GetHybrid());
		if(aptr->GetHybrid() == NO_HYBRID )
			continue;
		
		int nb_at = aptr->GetBondedAtoms(bonded_atoms);

		if(aptr->GetHybrid() == SP3_HYBRID)
		{
			dist =  1.08;
			if(nb_at >= 4) continue;
		}
		else if(aptr->GetHybrid() == SP2_HYBRID)
		{
            dist = 1.05;
            if(nb_at >= 3) continue;
		}
		else if(aptr->GetHybrid() == SP_HYBRID)
		{
			dist = 1.0;
            if(nb_at >= 2) continue;
		}
		
		HaResidue* pres = aptr->GetHostRes();
		HaMolecule* pmol = aptr->GetHostMol();

		AtomIteratorAtomGroup aitr_batoms(&bonded_atoms);
        HaAtom* aptr1; HaAtom* aptr2; HaAtom* aptr3; HaAtom* aptr4; 
		
		aptr1 =  aitr_batoms.GetFirstAtom();
		aptr2 = NULL;
		aptr3 = NULL;
		aptr4 = NULL;

		if(aptr1 != NULL)
		{
		   aptr2 =  aitr_batoms.GetNextAtom();
		}
		if(aptr2 != NULL)
		{
		   aptr3 =  aitr_batoms.GetNextAtom();
		}
		
		std::string h_name; 
		
		Vec3D dum1,dum2;
		
		if(aptr1 == NULL)
		{
			aptr1=pres->AddNewAtom();
			aptr1->SetElemNo(1);
    		new_bond = pmset->AddBond(aptr1, aptr );
	    	new_bond->DrawWire();
			
			h_name = pres->GetUniqueAtomName(1);
			
			if(!h_name.empty())
			{
				aptr1->SetName(h_name.c_str());
			}
			else
			{
				aptr1->SetName("H999");
			}
			
			dum1.SetX(aptr->GetX());
			dum1.SetY(aptr->GetY() + 1.0);
			dum1.SetZ(aptr->GetZ());
			
			dum2.SetX(aptr->GetX());
			dum2.SetY(aptr->GetY() + 1.0);
			dum2.SetZ(aptr->GetZ() + 1.0);
			
			Vec3D::SetAtomPos(aptr1, aptr, &dum1, &dum2, dist, PI/2.0, 0.0);
		}
		
		if(aptr2 == NULL)
		{
			aptr2=pres->AddNewAtom();
			aptr2->SetElemNo(1);
     		new_bond = pmset->AddBond(aptr2, aptr );
	    	new_bond->DrawWire();
			
			h_name = pres->GetUniqueAtomName(1);
			
			if(!h_name.empty())
			{
				aptr2->SetName(h_name.c_str());
			}
			else
			{
				aptr2->SetName("H999");
			}
			
			dum1.SetX(aptr->GetX() + 1.0);
			dum1.SetY(aptr->GetY() + 1.0);
			dum1.SetZ(aptr->GetZ() + 1.0);
			
			if( aptr->GetHybrid() == SP3_HYBRID)
			{
				val_angle = DEG_TO_RAD*109.5;
			}
			else if( aptr->GetHybrid() == SP2_HYBRID)
			{
				val_angle = PI/3.0;
			}
			else if( aptr->GetHybrid() == SP_HYBRID)
			{
				val_angle = PI;
			}
			
			Vec3D::SetAtomPos(aptr2, aptr, aptr1, &dum2, dist, val_angle, 0.0);
		}
		
		if( aptr->GetHybrid() == SP_HYBRID)
			continue;
		
		
		if(aptr3 == NULL)
		{
			aptr3=pres->AddNewAtom();
			aptr3->SetElemNo(1);
		    new_bond = pmset->AddBond(aptr3, aptr );
		    new_bond->DrawWire();
			
			h_name = pres->GetUniqueAtomName(1);
			
			if(!h_name.empty())
			{
				aptr3->SetName(h_name.c_str());
			}
			else
			{
				aptr3->SetName("H999");
			}
			
			if( aptr->GetHybrid() == SP3_HYBRID)
			{
				val_angle = DEG_TO_RAD*109.5;
				dih_angle = PI*(2.0/3.0);
			}
			else if( aptr->GetHybrid() == SP2_HYBRID)
			{
				val_angle = PI*(2.0/3.0);
				dih_angle = PI;
			}
			
			Vec3D::SetAtomPos(aptr3, aptr, aptr1, aptr2, dist, val_angle, dih_angle);
		}
		
		if( aptr->GetHybrid() == SP2_HYBRID)
		{
			continue;
		}
		
		aptr4=pres->AddNewAtom();
		aptr4->SetElemNo(1);
		new_bond = pmset->AddBond(aptr4, aptr );
		new_bond->DrawWire();
		
		h_name = pres->GetUniqueAtomName(1);
		
		if(!h_name.empty())
		{
			aptr4->SetName(h_name.c_str());
		}
		else
		{
			aptr4->SetName("H999");
		}
		
		if( aptr->GetHybrid() == SP3_HYBRID)
		{
			val_angle = DEG_TO_RAD*109.5;
			dih_angle = -PI*(2.0/3.0);
		}

		double dih_3= Vec3D::CalcDihedral(aptr3, aptr, aptr1, aptr2);

		if( dih_3 < 0 || dih_3 > PI) dih_angle = -dih_angle;
		
		Vec3D::SetAtomPos(aptr4, aptr, aptr1, aptr2, dist, val_angle, dih_angle); 
	}

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFApply );	
	return TRUE;
}

int MolEditor::SetHBondDonAccStatus(AtomContainer* p_at_coll)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	

	AtomIteratorGen aitr(p_at_coll);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		aptr->SetHBDonor(false);
		aptr->SetHBAcceptor(false);
		HaAtom* atempl = p_res_db->GetTemplateForAtom(aptr);
		if(atempl != NULL)
		{
			if(atempl->IsHBDonor()) aptr->SetHBDonor(true);
			if(atempl->IsHBAcceptor()) aptr->SetHBAcceptor(true);
		}		
	}
	return TRUE;
}

int MolEditor::SetFormalAtChrgFromTempl(MolSet* pmset)
{
	HaResidue* pres;
	HaAtom* aptr;

	ResidueIteratorMolSet ritr(pmset);
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes()) 
	{
		std::string res_name = pres->GetName();
		std::string name_mod = pres->GetNameModifier();
		AtomIteratorAtomGroup aitr(pres);

		if( pres->IsAmino() )
		{
			for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
			{
				aptr->SetCharge(0.0);
			}

			if( name_mod == "NT" || name_mod == "EPSILON_NT" || name_mod == "PROT_NT" )
			{
				HaAtom* aptr_n  = pres->GetAtomByName("N");
				if( aptr_n ) aptr_n->SetCharge(1.0);
			}

			if( name_mod == "CT" || name_mod == "EPSILON_CT" || name_mod == "PROT_CT" )
			{
				HaAtom* aptr_o  = pres->GetAtomByName("O");
				HaAtom* aptr_oxt  = pres->GetAtomByName("OXT");
				if( aptr_o && aptr_oxt ) 
				{
					aptr_o->SetCharge(-0.5);
					aptr_oxt->SetCharge(-0.5);
				}
			}

			if( (res_name == "HIS" && name_mod == "PROT") || res_name == "HIP" )
			{
				HaAtom* aptr_ne2  = pres->GetAtomByName("NE2");
				if( aptr_ne2 ) aptr_ne2->SetCharge(1.0);
			}

			if( res_name == "ASP" )
			{
				HaAtom* aptr_od1  = pres->GetAtomByName("OD1");
				HaAtom* aptr_od2  = pres->GetAtomByName("OD2");
				if( aptr_od1 && aptr_od2 )
				{
					aptr_od1->SetCharge(-0.5);
					aptr_od2->SetCharge(-0.5);
				}
			}

			if( res_name == "GLU" )
			{
				HaAtom* aptr_oe1  = pres->GetAtomByName("OE1");
				HaAtom* aptr_oe2  = pres->GetAtomByName("OE2");
				if( aptr_oe1 && aptr_oe2 )
				{
					aptr_oe1->SetCharge(-0.5);
					aptr_oe2->SetCharge(-0.5);
				}
			}

			if( res_name == "LYS" )
			{
				HaAtom* aptr_nz  = pres->GetAtomByName("NZ");
				if( aptr_nz ) aptr_nz->SetCharge(1.0);
			}

			if( res_name == "ARG" )
			{
				HaAtom* aptr_ne  = pres->GetAtomByName("NE");
				HaAtom* aptr_nh1  = pres->GetAtomByName("NH1");
				HaAtom* aptr_nh2  = pres->GetAtomByName("NH2");
				if( aptr_ne && aptr_nh1 && aptr_nh2 )
				{
					aptr_ne->SetCharge(0.333333);
					aptr_nh1->SetCharge(0.333333);
					aptr_nh2->SetCharge(0.333333);
				}
			}
		}
	}
	return TRUE;
}


int MolEditor::SetStdAtomicParams(MolSet* pmset, int at_params_type)
{
	HaResidue* group;
	HaChain*   chain;
	HaAtom* aptr;

	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	

	if(at_params_type & BACKBONE_CHRG)
	{
		ResidueIteratorMolSet ritr(pmset);
		for(group = ritr.GetFirstRes(); group; group = ritr.GetNextRes()) 
		{
			if(group->IsAmino()) // protein backbone
			{
				HaAtom* aptr_nh = group->GetAtomByName("HN");
				HaAtom* aptr_n  = group->GetAtomByName("N");
				HaAtom* aptr_c  = group->GetAtomByName("C");
				HaAtom* aptr_o  = group->GetAtomByName("O");
				HaAtom* aptr_ca = group->GetAtomByName("CA");
				if(aptr_nh != NULL)
				{
					aptr_nh->SetCharge(0.248);
					if(aptr_n != NULL) aptr_n->SetCharge(-0.520);
					if(aptr_ca != NULL) aptr_ca->SetCharge(0.246);
				}
				else if( !stricmp_loc(group->GetName(),"pro") )
				{
					if(aptr_n != NULL) aptr_n->SetCharge(-0.257);
					if(aptr_ca != NULL) aptr_ca->SetCharge(0.231);
				}
				else
				{
					if(aptr_n != NULL) aptr_n->SetCharge(-0.272);
					if(aptr_ca != NULL) aptr_ca->SetCharge(0.246);
				}
				if(aptr_c != NULL) aptr_c->SetCharge(0.526);
				if(aptr_o != NULL) aptr_o->SetCharge(-0.500);
			}
		}
	}	

	if(at_params_type & PROT_CHARGED_GROUPS_CHRG )
	{
		ResidueIteratorMolSet ritr(pmset);
		for(group = ritr.GetFirstRes(); group; group = ritr.GetNextRes()) 
		{
			if(!stricmp_loc(group->GetName(),"lys") ) // charged residues
			{
				aptr= group->GetAtomByName("NZ");
				if(aptr != NULL) aptr->SetCharge(+1.00);
			}
			else if(!stricmp_loc(group->GetName(),"arg") )
			{
				aptr= group->GetAtomByName("NH1");
				if(aptr != NULL) aptr->SetCharge(+0.5);
				aptr= group->GetAtomByName("NH2");
				if(aptr != NULL) aptr->SetCharge(+0.5);
			}
			else if(!stricmp_loc(group->GetName(),"glu") )
			{
				aptr= group->GetAtomByName("OE1");
				if(aptr != NULL) aptr->SetCharge(-0.5);
				aptr= group->GetAtomByName("OE2");
				if(aptr != NULL) aptr->SetCharge(-0.5);
			}
			else if(!stricmp_loc(group->GetName(),"asp") )
			{
				aptr= group->GetAtomByName("OD1");
				if(aptr != NULL) aptr->SetCharge(-0.5);
				aptr= group->GetAtomByName("OD2");
				if(aptr != NULL) aptr->SetCharge(-0.5);
			}
			else if(!stricmp_loc(group->GetName(),"hem") )
			{
				aptr= group->GetAtomByName("O1A");
				if(aptr != NULL) aptr->SetCharge(-0.5);
				aptr= group->GetAtomByName("O2A");
				if(aptr != NULL) aptr->SetCharge(-0.5);
				aptr= group->GetAtomByName("O1D");
				if(aptr != NULL) aptr->SetCharge(-0.5);
				aptr= group->GetAtomByName("O2D");
				if(aptr != NULL) aptr->SetCharge(-0.5);
			}
		}
	}

	AtomIteratorMolSet aitr(pmset);

	if(at_params_type & ZERO_CHRG )
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
		{
			if(aptr->Selected())
			{
				aptr->charge = 0.0;
			}
		}
	}

	if(at_params_type & ATOM_MASSES_ELEMENT )
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
		{
			if(aptr->Selected())
			{
				int elem = aptr->GetElemNo();
				aptr->SetMass( HaAtom::StdElemMass(elem) );
			}
		}
	}
	
	if(at_params_type & (AMBER_ALL_ATOM_CHRGS | AMBER_ALL_ATOM_FF_SYMBOLS | AMBER_ALL_ATOM_MASSES 
		                 | ATOM_HBOND_DA_STATUS | ATOM_ELEM_FROM_TEMPL ) )
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
		{
			if(aptr->Selected())
			{
				HaAtom* atempl = p_res_db->GetTemplateForAtom(aptr);
				if(atempl != NULL)
				{
					if( at_params_type & AMBER_ALL_ATOM_CHRGS )
						aptr->charge = atempl->GetCharge();
					if(at_params_type &  AMBER_ALL_ATOM_FF_SYMBOLS )
						aptr->SetFFSymbol(atempl->GetFFSymbol());
					if(at_params_type &  AMBER_ALL_ATOM_MASSES )
						aptr->SetMass(atempl->GetMass());
					if(at_params_type &  ATOM_ELEM_FROM_TEMPL )
						aptr->SetElemNo(atempl->GetElemNo());
					if(at_params_type &  ATOM_HBOND_DA_STATUS )
					{
						if(atempl->IsHBDonor())
							aptr->SetHBDonor(true);
						if(atempl->IsHBAcceptor())
							aptr->SetHBAcceptor(true);
					}
				}		
			}
		}
	}
	return TRUE;
}

int MolEditor::ClearAtomFFParams(MolSet* pmset)
{
	PrintLog(" MolEditor::ClearAtomFFParams() pt 1 \n");
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		aptr->SetCharge(0.0);
		aptr->SetFFSymbol("");
	}
	return TRUE;
}

int MolEditor::BondIfClose(HaAtom* sptr, HaAtom* dptr)
{
    HaBond  *bptr;
    double dx, dy, dz;
    double max, dist;
	
	if( sptr == NULL || dptr == NULL) return FALSE;

	MolSet* pmset = sptr->GetHostMolSet();

	if( pmset == NULL) return FALSE;
	if( pmset != dptr->GetHostMolSet()) return FALSE;

	dist = ElementArr[sptr->GetElemNo()].covalrad + 
		ElementArr[dptr->GetElemNo()].covalrad + 0.56;
	max = dist*dist;  
	
    dx = sptr->GetX()-dptr->GetX();   if( (dist=dx*dx)>max ) return FALSE;
    dy = sptr->GetY()-dptr->GetY();   if( (dist+=dy*dy)>max ) return FALSE;
    dz = sptr->GetZ()-dptr->GetZ();   if( (dist+=dz*dz)>max ) return FALSE;
	
    if( dist > MinBondDist )
    {   /* Reset Non-bonded flags! */
        sptr->flag &= ~NonBondFlag;
        dptr->flag &= ~NonBondFlag;
		
        bptr = pmset->AddBond(sptr,dptr );
		return TRUE;
    }
	return FALSE;
}

int MolEditor::CreateCovBonds(AtomContainer* pat_cont)
{
	if( pat_cont == NULL) return FALSE;
    double mx, my, mz; 
    int lx, ly, lz, ux, uy, uz;
    HaAtom  *aptr;

	BoxPartition part_table; // Table of Distributions of atoms
	                         // into quadrants between minimal and maximal 
	                         // coordinates of atoms of the molecule

	double MinX, MinY, MinZ, MaxX, MaxY, MaxZ;

	pat_cont->GetMinMaxCrd(MinX, MinY, MinZ, MaxX, MaxY, MaxZ);
	part_table.SetBoundaries(MinX - 0.05, MinY - 0.05, MinZ - 0.05, 
		                    MaxX + 0.05, MaxY + 0.05, MaxZ + 0.05);

	part_table.DistributePointsToCells(*pat_cont);
	part_table.SetRegionRad(max_bond_length);

    int ix, iy, iz, i;

	AtomIteratorGen aitr(pat_cont);

	AtomGroup close_atoms;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->flag |= NonBondFlag;
		
		part_table.GetNeighbors(*aptr,close_atoms);

		AtomIteratorAtomGroup aitr(&close_atoms);

		HaAtom* aptr2;
		MolSet* pmset = aptr->GetHostMolSet(); 
		for(aptr2 = aitr.GetFirstAtom(); aptr2; aptr2 = aitr.GetNextAtom())
		{			
			if( aptr2 > aptr && aptr2->GetHostMolSet() == pmset )
			{
				MolEditor::BondIfClose(aptr,aptr2);
			}
		}
	}
	return TRUE;
}

int MolEditor::FindAlphaHelix( HaMolecule* pmol, int pitch, int flag )
{
	HaChain  *chain;
    HaResidue  *group;
    HaResidue  *first;
    HaResidue  *ptr;
    int res,dist,prev;
	SecStructElement* fptr=NULL;

    /* Protein chains only! */
    ChainIteratorMolecule ch_itr(pmol);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{ 
		ResidueIteratorChain ritr_ch(chain);
		if( (first = ritr_ch.GetFirstRes()) && first->IsProtein() )
		{   
			prev = False; dist = 0;
			for( group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes() )
			{   
				if( group->IsAmino() )
				{   
					if( dist==pitch )
					{   
						res = first->HasBackBHBond(group);
						
						// find if there a H-bond between residue first and the residue 
						// which is pitch away along the chain, then set res=True

						if( res )
						{   
							if( prev )
							{   
								if( !(first->struc & HelixFlag) ) 
								{
									fptr=pmol->AddFeature();
									fptr->type=FeatHelix;
									fptr->init=first->serno;
									fptr->chain=chain->ident;
								}
								
								ptr = first;
								do {
									ptr->struc |= flag;
									fptr->term= ptr->serno;
									ptr = ptr->GetNextResInChain();
								} while( ptr != group );
							} 
							else 
								prev = True;
						} 
						else 
							prev = False;
					} 
				} 
				else 
					prev = False;
				
				if( group->struc & HelixFlag )
				{   
					first = group; prev = False; dist = 1;
				} 
				else if( dist==pitch )
				{   
					first = first->GetNextResInChain();
				} 
				else 
					dist++;
			}
		} 
	}
	return TRUE;
}

static void TestLadder(list<HaChain>::iterator chain_ref, list<HaChain>& Chains, 
					   HaResidue* prevri,HaResidue* curri, HaResidue* nexti)
{
	list<HaChain>::iterator chain= chain_ref;
   
	HaResidue* prevrj = NULL;
	HaResidue* currj = NULL;
	HaResidue* nextj = NULL;

    int count, result, found;

    // Already part of at least one ladder 
    found = curri->flag & SheetFlag;
    nextj = nexti->GetNextResInChain();

    while( True )
    {   
		if( nextj )
		{
			HaResidue* pres=(*chain).GetFirstRes();
            if( pres->IsProtein())
            {   
				count = 1;
                do {
                    if( count == 3 ) // skip first to residues in the chain to fill cprevj, currj,nextj
                    {   
						if( nexti->HasBackBHBond(currj) && currj->HasBackBHBond(prevri) )
                        {   
							result = ParaLadder;
                        } 
						else if( nextj->HasBackBHBond(curri) && curri->HasBackBHBond(prevrj) )
                        {   
							result = ParaLadder;
                        } 
						else if( nexti->HasBackBHBond(prevrj) &&  nextj->HasBackBHBond(prevri) )
                        {   
							result = AntiLadder;
                        } 
						else if( curri->HasBackBHBond(currj))
                        {   
							result = AntiLadder;
                        } 
						else result = NoLadder;
						
                        if( result )
                        {   
							curri->struc |= SheetFlag;
                            currj->struc |= SheetFlag;
                            if( found ) return;
                            found = True;
                        }
                    } 
					else 
						count++;
					
					prevrj = currj; 
					currj = nextj; 
					
				} while( (nextj = nextj->GetNextResInChain()) );
				
			} 
		}	
		chain++;
		if( chain != Chains.end() ) 
		{   
			nextj = (*chain).GetFirstRes();
		} 
		else 
			return;
    }
}

int MolEditor::FindBetaSheets(HaMolecule* pmol)
{
	ResidueIteratorMolecule ritr(pmol);
	ResidueIteratorMolecule ritr_next(pmol);

	MolSet* pmset = pmol->GetHostMolSet();

	HaChain* chain_i = NULL;

	HaResidue* pr_prev_i = NULL;                     // Previous residue 
	HaResidue* pr_i = NULL;                          // Current Residue
	HaResidue* pr_next_i = ritr_next.GetFirstRes();  // Next Residue in chain

	AtomGroup hb_acc_i;
	AtomGroup hb_acc_next_i;

	HaAtom* nh_i      = NULL;  // N of amino group of the current residue 
	HaAtom* nh_next_i = NULL;  // N of amino group of the next residue

	HaAtom* oc_i      = NULL;  // O of C=O group of the current residue 
	HaAtom* oc_prev_i = NULL;  // O of C=O group of the previous residue

	for( pr_i = ritr.GetFirstRes(); pr_i; pr_i = ritr.GetNextRes())
	{
		if(chain_i != pr_i->GetHostChain()) // new chain
		{
			chain_i = pr_i->GetHostChain();
			pr_prev_i = NULL;
			oc_prev_i = NULL;
		}
		else
		{
			oc_prev_i = pr_prev_i->GetAtomByName("O");  
		}

		pr_next_i = ritr_next.GetNextRes();
		if( pr_next_i != NULL && pr_next_i->GetHostChain() != chain_i)
		{
			pr_next_i = NULL;
		}

		nh_i = pr_i->GetAtomByName("N");
		oc_i = pr_i->GetAtomByName("O");

		hb_acc_i.clear();
		if(nh_i) nh_i->GetHBondAcc(hb_acc_i);
		
		hb_acc_next_i.clear();
		if( pr_next_i == NULL )
		{
			nh_next_i = NULL;
		}
		else
		{
			nh_next_i = pr_next_i->GetAtomByName("N");
			if( nh_next_i != NULL) nh_next_i->GetHBondAcc(hb_acc_next_i);
		}
		
		int k;
		int nb_i = hb_acc_i.size();

// Check for antiparallel beta sheet pattern ( i -> j && j->i )
		for(k = 0; k < nb_i; k++)  
		{
			if(oc_i == NULL) continue;
			HaAtom* oc_j = hb_acc_i[k];
			HaResidue* pr_j = oc_j->GetHostRes();
			HaAtom* nh_j = pr_j->GetAtomByName("N");
			if( !pr_j->IsAmino() || nh_j == NULL) continue;
	
			if( pmset->AreHBonded(nh_j,oc_i) ) 
			{
				pr_i->struc |= SheetFlag;
                pr_j->struc |= SheetFlag;
			}
		}

		int nb_next_i = hb_acc_next_i.size();
		for(k = 0; k < nb_next_i; k++)  
		{
			if(oc_prev_i == NULL) continue;
				
			HaAtom* oc_prev_j = hb_acc_next_i[k];
			HaResidue* pr_prev_j = oc_prev_j->GetHostRes();
			HaAtom* nh_prev_j = pr_prev_j->GetAtomByName("N");

			if(nh_prev_j != NULL)
			{
// Check for parallel beta sheet pattern ( (i+1) -> (j-1) && (j-1) ->(i-1) )
				if( pmset->AreHBonded(nh_prev_j,oc_prev_i) )
				{
					pr_i->struc |= SheetFlag;
					pr_prev_j->struc |= SheetFlag;
				}
			}

			HaResidue* pr_j      = pr_prev_j->GetNextResInChain();
			if( pr_j == NULL) continue;
			HaResidue* pr_next_j = pr_j->GetNextResInChain();
			if( pr_next_j == NULL) continue;
			HaAtom* nh_next_j = pr_next_j->GetAtomByName("N");
			if( nh_next_j == NULL) continue;
	
// Check for antiparallel beta sheet pattern ( (i+1) -> (j-1) && (j+1) ->(i-1) )

			if( pmset->AreHBonded(nh_next_j,oc_prev_i) ) 
			{
				pr_i->struc |= SheetFlag;
                pr_j->struc |= SheetFlag;
			}
		}
		pr_prev_i = pr_i;
	}

// Fill FeatureList

	chain_i = NULL;

//	HaMolecule* pMol= NULL;
	SecStructElement* feat_curr = NULL;
	for( pr_i = ritr.GetFirstRes(); pr_i; pr_i = ritr.GetNextRes())
	{
		if(chain_i != pr_i->GetHostChain()) // new chain
		{
			chain_i = pr_i->GetHostChain();
			feat_curr = NULL;
//			if( pMol != chain_i->GetHostMol() ) // new molecule
//			{
//				pMol->DeleteFeatures(FeatSheet);
//			}
		}
		if( pr_i->struc & SheetFlag)
		{
			if( feat_curr != NULL) 
			{
				feat_curr->term = pr_i->GetSerNo();
			}
			else
			{
				feat_curr = pmol->AddFeature();
				feat_curr->type = FeatSheet;
				feat_curr->chain = (*chain_i).ident;
				feat_curr->init = pr_i->GetSerNo();
				feat_curr->term = pr_i->GetSerNo();
			}
		}
		else
		{
			feat_curr= NULL;
		}
	}
	return TRUE;
}

#if 0 

void MolEditor::FindBetaSheets_old(HaMolecule* pmol)
{
    int ladder;
    int count;
	SecStructElement* fptr=NULL;

	HaResidue* prevri = NULL;
	HaResidue* curri = NULL;
	HaResidue* nexti = NULL;

	list<HaChain>::iterator chain;
    for(chain=Chains.begin(); chain != Chains.end(); chain++)
	{
		if( (nexti = (*chain).GetFirstRes()) )
		{
			if( nexti->IsProtein() )
			{   
				count = 1;
				ladder = False;
				do // cycle on residues in the chain (ptr current residue nexti)
				{
					
					if( count == 3 )
					{   
						TestLadder( chain,Chains, prevri,curri,nexti );
						if( curri->struc & SheetFlag )
						{   
							if( !ladder )
							{   
								fptr=AddFeature();
								fptr->type=FeatSheet;
								fptr->chain=(*chain).ident;
								ladder = True;
							}
						} 
						else 
						{
							ladder = False;
						}
					} 
					else
					{
						count++;
					}
					
					prevri = curri; 
					curri = nexti;  
				} 
				while( (nexti = nexti->GetNextResInChain()) );
			} 
		}
	} // End of a cycle on chains
}

#endif

int MolEditor::FindTurnStructure(HaMolecule* pmol)
{
    static HaAtom  *aptr[5];
    HaChain  *chain;
    HaResidue  *group;
    HaResidue  *prev;
    HaAtom  *ptr;
    double ux,uy,uz,mu;
    double vx,vy,vz,mv;
    int i,found,len;
    double CosKappa;
	SecStructElement* fptr;

    ChainIteratorMolecule ch_itr(pmol);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
		if( (!chain->res_map.empty()) && chain->GetFirstRes()->IsProtein() )
		{   
			len = 0;  found = False;
			ResidueIteratorChain ritr_ch(chain);
			for( group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes() )
			{    
				ptr = group->GetAtomByName("CA");
				if( ptr && (ptr->flag&BreakFlag) )
				{   
					found = False;
					len = 0;
				} 
				else if( len==5 )
				{   
					for( i=0; i<4; i++ )
						aptr[i] = aptr[i+1];
					len = 4;
				} 
				else if( len==2 )
					prev = group;
				
				aptr[len++] = ptr;
				if( len==5 ) 
				{   
					if( !(prev->struc&(HelixFlag|SheetFlag)) &&
						aptr[0] && aptr[2] && aptr[4] )
					{   
						ux = aptr[2]->GetX() - aptr[0]->GetX();
						uy = aptr[2]->GetY() - aptr[0]->GetY();
						uz = aptr[2]->GetZ() - aptr[0]->GetZ();
						
						vx = aptr[4]->GetX() - aptr[2]->GetX();
						vy = aptr[4]->GetY() - aptr[2]->GetY();
						vz = aptr[4]->GetZ() - aptr[2]->GetZ();
						
						mu = ux*ux + uy*uy + uz*uz;
						mv = vx*vx + vz*vz + vy*vy;
						if( mu > 0.00001 && mv > 0.00001)
						{   
							CosKappa = ux*vx + uy*vy + uz*vz;
							CosKappa /= sqrt( mu*mv );
							if( CosKappa < Cos70Deg )
							{   
								if( !found )
								{
									fptr=pmol->AddFeature();
									fptr->type=FeatTurn;
									fptr->chain = chain->ident;
									fptr->init=prev->serno;
								}
								prev->struc |= TurnFlag;
								fptr->term=prev->serno;
							}
						}
					}
					found = prev->struc&TurnFlag;
					prev = prev->GetNextResInChain();
				} /* len==5 */
			}
		}
	}
	return TRUE;
}

int MolEditor::FindBetaTurns(HaMolecule* pmol)
{
    static HaAtom  *aptr[4];
    HaChain  *chain;
    HaResidue  *group;
    HaResidue  *prev;
    HaResidue  *next;
    HaAtom  *ptr;
    double dx,dy,dz;
    int found,len;
    int flag;
	SecStructElement* fptr;

    ChainIteratorMolecule ch_itr(pmol);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
        if( (!chain->res_map.empty()) && chain->GetFirstRes()->IsProtein() )
        {   
			prev = chain->GetFirstRes();  
			len = 0;  found = False;
			ResidueIteratorChain ritr_ch(chain);
			for( next = ritr_ch.GetFirstRes(); next; next = ritr_ch.GetNextRes() )
			{   
				ptr = next->GetAtomByName("CA");
				if( ptr && (ptr->flag&BreakFlag) )
				{   
					found = False;
					prev = next;
					len = 0;
				} 
				else if( len==4 )
				{   
					aptr[0] = aptr[1];
					aptr[1] = aptr[2];
					aptr[2] = aptr[3];
					aptr[3] = ptr;
					
				} 
				else 
					aptr[len++] = ptr;
				if( len==4 ) 
				{   
					flag = False;
					if( aptr[0] && aptr[3] )
					{   
						dx = aptr[3]->GetX() - aptr[0]->GetX();
						dy = aptr[3]->GetY() - aptr[0]->GetY();
						dz = aptr[3]->GetZ() - aptr[0]->GetZ();
						if( ( dx*dx + dy*dy + dz*dz ) < (7.0*7.0) )
						{   
							group = prev;
							while( group!=next->GetNextResInChain() )
							{   
								if( !(group->struc&(HelixFlag|SheetFlag)) )
								{   
									group->struc |= TurnFlag;
									flag = True;
								}
								group = group->GetNextResInChain();
							}
							if( !found && flag ) 
							{	
								fptr=pmol->AddFeature();
								fptr->type= FeatTurn;
								fptr->chain = chain->ident;
								fptr->init=group->serno;
								fptr->term=group->serno;
							}
						}
					}
					prev = prev->GetNextResInChain();   
					found = flag;
				} /* len==4 */
			}
        }
	}
	return TRUE;
}

int
MolEditor::DetermineSecStructure(HaMolecule* pmol, int flag )
{
    HaChain  *chain;
    HaResidue  *group;
	
	MolSet*  pmset = pmol->GetHostMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

    if( !pmset->HBonds_found ) p_mol_editor->CalcHBonds(pmset);
	
    if( pmol->IsSecStructFound() )
	{
		ChainIteratorMolecule ch_itr(pmol);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for( group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes() )
				group->struc = 0;
		}	
	}
	pmol->structsource = SourceCalc;
	pmol->Features.clear();
	
	if( pmol->GetNHBonds() )
	{   
		FindAlphaHelix(pmol,4,Helix4Flag);
		FindBetaSheets(pmol);
		FindAlphaHelix(pmol,3,Helix3Flag);
		FindAlphaHelix(pmol,5,Helix5Flag);
		
		if( !flag )
		{   
			FindTurnStructure(pmol);
		} 
		else 
			FindBetaTurns(pmol);
	}
	
	pmol->sec_struct_found = TRUE;

	PrintLog("Number of Helices ... %d\n", pmol->GetNumFeatures(FeatHelix));
	PrintLog("Number of Strands ... %d\n", pmol->GetNumFeatures(FeatSheet));
	PrintLog("Number of Turns ..... %d\n", pmol->GetNumFeatures(FeatTurn));

	return TRUE;
}

void MolEditor::SetAlphaHelix(MolSet* pmset)
{
	ResidueIteratorMolSet ritr(pmset);
	
	HaResidue* pres1 = ritr.GetFirstRes();
	if( pres1 == NULL) return;
	HaResidue* pres2;
	for( pres2 = ritr.GetNextRes(); pres2; pres1 = pres2, pres2 = ritr.GetNextRes())
	{
		int sel_1 = FALSE;
		int sel_2 = FALSE;
		if(pres1->size() == 0) continue;
		HaAtom* aptr1 = (*pres1)[0];
		if(!aptr1->Selected() ) continue;
		if(pres2->size() == 0) continue;
		HaAtom* aptr2 = (*pres2)[0];
		if(!aptr2->Selected() ) continue;
		if(!pres1->IsAmino()) continue;
		if(!pres2->IsAmino()) continue;
		
		HaAtom* aptr_n1  = pres1->GetAtomByName("N");  if(aptr_n1 == NULL) continue;
		HaAtom* aptr_ca1 = pres1->GetAtomByName("CA"); if(aptr_ca1 == NULL) continue;
		HaAtom* aptr_c1  = pres1->GetAtomByName("C");  if(aptr_c1 == NULL) continue;
		HaAtom* aptr_n2  = pres2->GetAtomByName("N");  if(aptr_n2 == NULL) continue;
		HaAtom* aptr_ca2 = pres2->GetAtomByName("CA"); if(aptr_ca1 == NULL) continue;
		HaAtom* aptr_c2  = pres2->GetAtomByName("C");  if(aptr_c2 == NULL) continue;
		
		if( !aptr_c1->IsBonded(*aptr_n2) ) continue;

		double psi = -45.0*DEG_TO_RAD;
		double phi = -60.0*DEG_TO_RAD;

		SetTorsion(aptr_n1,aptr_ca1,aptr_c1,aptr_n2, psi);
		SetTorsion(aptr_c1,aptr_n2,aptr_ca2,aptr_c2, phi);		
	}
}

int MolEditor::CenterAtOrigin(AtomContainer* pat_cont)
{
	if( pat_cont == NULL) return FALSE;

	HaAtom* aptr;
	double xmin, ymin, zmin, xmax, ymax, zmax;
	pat_cont->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);
	double shiftx = - (xmax + xmin)/2.0;
	double shifty = - (ymax + ymin)/2.0;
	double shiftz = - (zmax + zmin)/2.0;

    AtomIteratorGen aitr(pat_cont);

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		double xx = aptr->GetX() + shiftx;
		double yy = aptr->GetY() + shifty;
		double zz = aptr->GetZ() + shiftz;
		aptr->SetX(xx);
		aptr->SetY(yy);
		aptr->SetZ(zz);
	}
	return TRUE;
}

int MolEditor::CenterAtOriginWithRad(AtomContainer* pat_cont)
{
	if( pat_cont == NULL) return FALSE;

	int na = pat_cont->GetNAtoms();
	if( na == NULL) return FALSE;

	double xmin = 1.0e10;
	double ymin = 1.0e10;
	double zmin = 1.0e10;
	double xmax = -1.0e10;
	double ymax = -1.0e10;
	double zmax = -1.0e10;

	HaAtom* aptr;
	AtomIteratorGen aitr(pat_cont);
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		int elemno = aptr->GetElemNo();
		double rad = aptr->vdw_rad;
		if( rad < 0.01 )
		{
			rad = aptr->HaAtom::ElemVDWRadius(elemno);
		}
		double xmin_at = aptr->GetX() - rad;
		double ymin_at = aptr->GetY() - rad;
		double zmin_at = aptr->GetZ() - rad;
		double xmax_at = aptr->GetX() + rad;
		double ymax_at = aptr->GetY() + rad;
		double zmax_at = aptr->GetZ() + rad;
		if( xmin_at < xmin) xmin = xmin_at;
		if( ymin_at < ymin) ymin = ymin_at;
		if( zmin_at < zmin) zmin = zmin_at;
		if( xmax_at > xmax) xmax = xmax_at;
		if( ymax_at > ymax) ymax = ymax_at;
		if( zmax_at > zmax) zmax = zmax_at;
	}

	double shiftx = - (xmax + xmin)/2.0;
	double shifty = - (ymax + ymin)/2.0;
	double shiftz = - (zmax + zmin)/2.0;

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		double xx = aptr->GetX() + shiftx;
		double yy = aptr->GetY() + shifty;
		double zz = aptr->GetZ() + shiftz;
		aptr->SetX(xx);
		aptr->SetY(yy);
		aptr->SetZ(zz);
	}
	return TRUE;
}

int MolEditor::Solvate(MolSet* pmset)
{
	MolSet* cur_mol_set = GetCurMolSet();

	HaResDB* p_res_db = HaResDB::GetDefaultResDB();
	std::string solv_fname = pApp->res_db_dir + solv_name + ".hlm";
	FILE* solv_file = fopen(solv_fname.c_str(),"r");
	if(solv_file == NULL)
	{
		PrintLog("Error In MolEditor::Solvate() \n");
	    PrintLog(" No solvent file %s in the residues_db directory", solv_fname.c_str() );
		return FALSE;
	}
	fclose(solv_file);

	this->CenterAtOriginWithRad(pmset);

	double xmin, ymin, zmin, xmax, ymax, zmax;
	pmset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

	xmax += solv_buffer_dist;
	ymax += solv_buffer_dist;
	zmax += solv_buffer_dist;

	MolSet* solvent = new MolSet();
	solvent->LoadHarlemFile(solv_fname.c_str());


    if( !solvent->per_bc->IsSet() )
	{
		PrintLog("MolEditor::Solvate() \n");
	    PrintLog("Solvent file %s Does not have periodic box information \n", solv_fname.c_str());
	      
//		delete solvent; // some errors when deleting solvent  TO FIX ?

		return False;		
	}

	int nx = (int)((2.0*xmax)/solvent->per_bc->GetA()); nx++;
	int ny = (int)((2.0*ymax)/solvent->per_bc->GetB()); ny++;
	int nz = (int)((2.0*zmax)/solvent->per_bc->GetC()); nz++;

	if( nx > 1 || ny > 1 || nz > 1)
 	{
		this->ReplicatePeriodBox(solvent,nx,ny,nz);
	}

	this->WrapToUnitCell(solvent,solvent->per_bc);

	AtomIteratorMolSet aitr(pmset);
	AtomGroup old_atoms;
	HaAtom* aptr;

	Vec3D tr;

	double a = solvent->per_bc->GetA();
	double b = solvent->per_bc->GetB();
	double c = solvent->per_bc->GetC();

	int i;
	for( i = 0; i < 3; i++)
	{
		tr[i] =   0.5*solvent->per_bc->ucell[0][i]; 
		tr[i] +=  0.5*solvent->per_bc->ucell[1][i];
		tr[i] +=  0.5*solvent->per_bc->ucell[2][i];
	}

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->SetCoordFrom((Vec3D)(*aptr) + tr );
		old_atoms.InsertAtom(aptr);
	}
	AtomIteratorAtomGroup aitr_old_atoms(&old_atoms);

	HaMolecule* pMol = pmset->AddNewMolecule();
	pMol->SetObjName("Solvent");
	MoleculesType::iterator mol_itr;
	HaChain* chain;
	HaResidue* group;

	HaChain* chain_cur = pMol->AddChain(' ');

	int nres = 0;
	for(mol_itr = solvent->HostMolecules.begin(); mol_itr != solvent->HostMolecules.end(); mol_itr++)
	{
		AtomAtomMap at_map;
		AtomAtomMap::iterator mitr;
			
		ChainIteratorMolecule chitr(*mol_itr);
		for( chain= chitr.GetFirstChain();chain; chain= chitr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes();group; group = ritr_ch.GetNextRes())
			{
				int overlap = False;
		        for(aptr = aitr_old_atoms.GetFirstAtom(); aptr; aptr = aitr_old_atoms.GetNextAtom())
				{
					if( group->IsWithinRadius(aptr, 25.0))
						overlap = True;
				}
				if(overlap)
					continue;

				nres++;
				HaResidue* new_res = chain_cur->AddResidue(nres);
				new_res->SetParamFrom(*group);
				new_res->serno = nres;

				HaAtom* aptr_new;
				AtomIteratorAtomGroup aitr_group(group);

				for(aptr = aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					aptr_new = new_res->AddNewAtom();
					aptr_new->SetParamFrom(*aptr);
					at_map[aptr] = aptr_new ;
				}
			}
		}

		HaBond* bptr;
		HaBond* bptr_new;

		AtomIteratorMolecule aitr(*mol_itr);
		for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
		{
			HaAtom::BondIterator bitr = aptr->Bonds_begin();
			for( ; bitr != aptr->Bonds_end(); ++bitr )
			{
				bptr = (*bitr);

				mitr = at_map.find(bptr->srcatom);
				if(mitr == at_map.end()) continue;
				HaAtom* aptr1= (*mitr).second;
				mitr = at_map.find(bptr->dstatom);
				if(mitr == at_map.end()) continue;
				HaAtom* aptr2=(*mitr).second;
//				if( aptr2 < aptr1 ) continue;
				bptr_new = pmset->AddBond( aptr1,aptr2 );
				bptr_new->SetParamFrom(*bptr);
			}
		}
	}

	pmset->per_bc->SetBox(solvent->per_bc->GetA(),solvent->per_bc->GetB(),solvent->per_bc->GetC());
    delete solvent;

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->InitialTransform();
	    pview->DefaultRepresentation();
	}
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );	

	SetCurMolSet(cur_mol_set);

	return TRUE;
}

int MolEditor::ReplicatePeriodBox(MolSet* pmset, int nx, int ny, int nz)
{
	 if( pmset == NULL) return FALSE;
	 if( nx <= 0 || ny <= 0 || nz <= 0)
	 {
		 ErrorInMod("MolEditor::ReplicatePeriodBox",  
			        " One of the multiplies < 1 ");
		 return FALSE;
	 }
	 if( pmset->GetNMol() == 0 )
	 {
		 ErrorInMod("MolEditor::ReplicatePeriodBox",  
			        " No Molecules in the molset");
		 return FALSE;
	 }
     list<HaResidue*> res_list;
	 list<HaResidue*>::iterator resl_itr;
	 HaResidue* group;
	 HaChain* chain;
	 HaAtom* aptr;
	 
	 ResidueIteratorMolSet ritr_mset(pmset);
	 for(group = ritr_mset.GetFirstRes(); group; group = ritr_mset.GetNextRes())
	 {
		 res_list.push_back(group);
	 }
	 int i,j,k;
	 int res_num = 0;
	 HaMolecule* pMol = pmset->HostMolecules[0];
	 char ident_max = pMol->GetChainIdentMax();
	 ident_max++;
	 HaChain* chain_cur = pMol->AddChain(ident_max);

	 ResidueIteratorMolecule ritr_mol(pMol);
	 for(group = ritr_mol.GetFirstRes(); group;group =  ritr_mol.GetNextRes())
	 {
		 if(group->GetSerNo() > res_num) res_num = group->GetSerNo();
	 }
	 list<HaBond*> bond_list;
	 list<HaBond*>::iterator bondl_itr;

	 HaBond* bptr;
	 BondIteratorMolSet bitr(pmset);
	 for(bptr=bitr.GetFirstBond(); bptr; bptr= bitr.GetNextBond())
	 {	
		 bond_list.push_back(bptr);
	 }

	 for( k = 0; k < nz; k++)
	 {
		 for( j = 0; j < ny; j++)
		 {
			 for( i = 0; i < nx; i++)
			 {
				  if( i == 0 && j == 0 && k == 0 ) continue;
				  
				  AtomAtomMap frag_at_map;
				  AtomAtomMap::iterator mitr;
				  
				  for(resl_itr = res_list.begin(); resl_itr != res_list.end(); resl_itr++)
				  {
					  res_num++;
					  HaResidue* res_ref = *resl_itr;
					  HaResidue* res_new = chain_cur->AddResidue(res_num);
					  res_new->SetParamFrom(*res_ref);
					  HaAtom* aptr_new;
					  AtomIteratorAtomGroup aitr_res_ref(res_ref);
					  for( aptr = aitr_res_ref.GetFirstAtom(); aptr; aptr = aitr_res_ref.GetNextAtom())
					  {
						   aptr_new = res_new->AddNewAtom();
						   frag_at_map[aptr] = aptr_new ;
						   aptr_new->SetParamFrom(*aptr);
						   double xx = aptr->GetX();
						   double yy = aptr->GetY();
						   double zz = aptr->GetZ();
						   aptr_new->SetX(xx + i * pmset->per_bc->GetA());
						   aptr_new->SetY(yy + j * pmset->per_bc->GetB());
						   aptr_new->SetZ(zz + k * pmset->per_bc->GetC());
					  }
				  }
				  
				  for(bondl_itr = bond_list.begin(); bondl_itr != bond_list.end(); bondl_itr++)
				  {	  
					  bptr = *bondl_itr;
					  mitr = frag_at_map.find(bptr->srcatom);
					  if(mitr == frag_at_map.end())
						  continue;
					  HaAtom* aptr1_new = (*mitr).second;
					  mitr = frag_at_map.find(bptr->dstatom);
					  if(mitr == frag_at_map.end())
						  continue;
					  HaAtom* aptr2_new = (*mitr).second;
					  HaBond* bptr_new = pmset->AddBond(aptr1_new, aptr2_new );
					  bptr_new->SetParamFrom(*bptr);
				  }
			 }
		 }
	 }

	 double shiftx = - pmset->per_bc->GetA() * (nx - 1)/2.0;
	 double shifty = - pmset->per_bc->GetB() * (ny - 1)/2.0;
	 double shiftz = - pmset->per_bc->GetC() * (nz - 1)/2.0;
 
     AtomIteratorMolSet aitr(pmset);

	 for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	 {
		 double xx = aptr->GetX();
		 double yy = aptr->GetY();
		 double zz = aptr->GetZ();
		 aptr->SetX(xx + shiftx);
		 aptr->SetY(yy + shifty);
		 aptr->SetZ(zz + shiftz);
	 }
		
	 pmset->per_bc->SetBox( pmset->per_bc->GetA()*nx, 
		                    pmset->per_bc->GetB()*ny,
							pmset->per_bc->GetC()*nz );

	 HaMolView* pview = pmset->GetActiveMolView();
	 if(pview)pview->CPKColourAttrib();
	 pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );	

	 return True;
}

int MolEditor::WrapToUnitCell(AtomContainer* pat_coll, PeriodicUnitInfo* per_info)
{
	if( pat_coll == NULL || per_info == NULL) return FALSE;
	if( !per_info->IsSet() )
	{
		PrintLog("Error in MolEditor::WrapToUnitCell() \n");
		PrintLog("Periodic Box Info is not set \n");
		return FALSE;
	}

	vector<AtomGroup> mols;
	SplitToMolecules(pat_coll,mols);

	int nmol = mols.size();

	int imol,i,j;
	HaAtom* aptr;

	for(imol = 0; imol < nmol; imol++)
	{
		Vec3D f,g;

		f.SetZeros();
		int na = mols[imol].size();

		if( na == 0 )
		{
			PrintLog("MolEditor::WrapToUnitCell() mol %d  na = %d\n",i,na);
		}

		for(j = 0; j < na; j++)
		{
			aptr = mols[imol][j];
			f[0] += Vec3D::DotProduct(*aptr,per_info->recip_ucell[0]);
			f[1] += Vec3D::DotProduct(*aptr,per_info->recip_ucell[1]);
			f[2] += Vec3D::DotProduct(*aptr,per_info->recip_ucell[2]);
		}	

		for(i = 0; i < 3; i++)
		{
			f[i] = f[i]/na - 0.5;
		}

		for(i = 0; i < 3; i++)
		{
			int n = (int)f[i];
			double dn = (double) n; 
			double nearest_n;
			nearest_n = dn;
			if( (dn < f[i]) && ((f[i] - dn) > 0.5) ) nearest_n = dn + 1.0; 
			if( (dn > f[i]) && ((dn - f[i]) > 0.5) ) nearest_n = dn - 1.0;
			g[i] = f[i] - nearest_n;
		}

		if ( fabs(f[0] - g[0]) > 0.00001 || fabs(f[1] - g[1]) > 0.00001 || fabs(f[2] - g[2]) > 0.00001 )
		{
			Vec3D tr;
			Vec3D gmf = g - f;

			for( i = 0; i < 3; i++)
			{
				tr[i] = gmf[0]* per_info->ucell[0][i] + gmf[1]* per_info->ucell[1][i] + gmf[2]* per_info->ucell[2][i];
			}

			for( j = 0; j < na; j++)
			{
				aptr = mols[imol][j];
				aptr->SetCoordFrom((Vec3D)(*aptr) + tr);
			}
		}
	}
	return TRUE;
}

int MolEditor::DeleteOverlapMols(MolSet* pmset, AtomGroup& at_coll)
{
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;

	std::set<HaAtom*> atset_sel;
	std::set<HaAtom*> atset_unsel;
	std::set<HaAtom*> atset_ovlp;

	AtomIteratorAtomGroup atc_itr(&at_coll);

	for(aptr = atc_itr.GetFirstAtom();aptr; aptr = atc_itr.GetNextAtom())
	{
		if(aptr->Selected()) atset_sel.insert(aptr);
	}

	for( aptr = aitr.GetFirstAtom();aptr; aptr = aitr.GetNextAtom())
	{
		if( atset_sel.find(aptr) == atset_sel.end()) atset_unsel.insert(aptr);
	}

	if( atset_sel.size() == 0 || atset_unsel.size() == 0) return FALSE;
	
	std::set<HaAtom*>::iterator aitr1;
	std::set<HaAtom*>::iterator aitr2;

	for( aitr1 = atset_sel.begin(); aitr1 != atset_sel.end(); aitr1++)
	{
		aptr1 = (*aitr1);
		double rad1 = HaAtom::ElemVDWRadius(aptr1->GetElemNo());
		for( aitr2 = atset_unsel.begin(); aitr2 != atset_unsel.end(); aitr2++)
		{
			aptr2 = (*aitr2);
			double rad2 = HaAtom::ElemVDWRadius(aptr2->GetElemNo());
			double dist = Vec3D::CalcDistance(aptr1,aptr2);
			if( dist < (rad1+rad2)*0.9 ) atset_ovlp.insert(aptr2);
		}
	}

	int i,na;

	if( atset_ovlp.empty()) return TRUE;

	AtomGroup bonded_atoms;
	std::set<HaAtom*> atset_ovlp_2 = atset_ovlp;
	

	for(;;)
	{
		for( aitr1 = atset_ovlp.begin(); aitr1 != atset_ovlp.end(); aitr1++)
		{
			aptr = (*aitr1);
			aptr->GetBondedAtoms(bonded_atoms);
			na = bonded_atoms.size();
			for(i = 0; i < na; i++)
			{
				aptr2 = bonded_atoms[i];
				if( atset_sel.find(aptr2) != atset_sel.end())
				{
					PrintLog("Error in MolEditor::RemoveOverlapMols() \n");
					PrintLog("Excluded molecules are connected to selection \n");
					return FALSE;
				}
				atset_ovlp_2.insert(aptr2);
			}
		}
		if( atset_ovlp_2.size() == atset_ovlp.size()) break;
			atset_ovlp = atset_ovlp_2;
	}

	AtomGroup atoms_del;
	na = atset_ovlp.size();
	atoms_del.reserve(na);

	for( aitr1 = atset_ovlp.begin(); aitr1 != atset_ovlp.end(); aitr1++)
	{
		aptr = (*aitr1);
		atoms_del.push_back(aptr);
	}
	pmset->DeleteAtoms(atoms_del);
	return TRUE;
}

int MolEditor::SplitToMolecules(AtomContainer* p_at_coll, vector<AtomGroup>& mols)
{
	mols.clear();
	if( p_at_coll == NULL ) return FALSE;

	int na = p_at_coll->GetNAtoms();

	set<HaAtom*, less<HaAtom*> > atoms_left;
	set<HaAtom*, less<HaAtom*> > all_atoms;

	vector< set<HaAtom*, less<HaAtom*> > > clusters;
	set<HaAtom*> empty_atset;
	set<HaAtom*>* p_last_atset = NULL;

	AtomIteratorGen aitr(p_at_coll);
	HaAtom* aptr;

	std::queue< HaAtom* > conn_atoms;

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		atoms_left.insert(aptr);
		all_atoms.insert(aptr);
	}

	AtomGroup bonded_atoms;
	int i,j;

	for(;;)
	{
		if( conn_atoms.empty() )
		{
			if( atoms_left.empty() ) break;

			set<HaAtom*, less<HaAtom*> >::iterator atoms_left_itr = atoms_left.begin();
			aptr = (*atoms_left_itr);
			conn_atoms.push(aptr);
			atoms_left.erase(atoms_left_itr);
			clusters.push_back(empty_atset);
			p_last_atset = &(clusters.back());
		}
		aptr = conn_atoms.front();
		aptr->GetBondedAtoms(bonded_atoms);
		int na = bonded_atoms.GetNAtoms();
		for(i = 0; i < na; i++)
		{
			HaAtom* aptr_b = bonded_atoms[i];
			if( (all_atoms.count(aptr_b) > 0) && (atoms_left.count(aptr_b) > 0) )
			{
				conn_atoms.push(aptr_b);
				atoms_left.erase(aptr_b);
			}
		}
		p_last_atset->insert(aptr);
		conn_atoms.pop();
	}

	int nc = clusters.size();
//	PrintLog("MolEditor::SplitToMolecules() : \n");
//	PrintLog("Number of clusters: %d \n", nc );

	mols.resize(nc);
	for(i = 0; i < nc; i++)
	{
		na = clusters[i].size();
		mols[i].resize(na);
		j = 0;
		set<HaAtom*, less<HaAtom*> >::iterator asitr;
		for( asitr = clusters[i].begin(); asitr != clusters[i].end(); asitr++)
		{
			mols[i][j] = (*asitr);
			j++;
		}
	}
	return TRUE;
}

int MolEditor::CenterSoluteInSolvent(MolSet* pmset)
{
	if(!pmset->per_bc->IsSet())
	{
		PrintLog("Error in MolEditor::CenterSoluteInSolvent() \n");
		PrintLog("Periodicity information is not set \n");
		return FALSE;
	}
	
	AtomIteratorGen aitr(pmset);
	HaAtom* aptr;

	int n_solute = 0;
	double xc = 0.0;
	double yc = 0.0;
	double zc = 0.0;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		HaResidue* pres = aptr->GetHostRes();
		HaMolecule* pmol = aptr->GetHostMol();
		if( pres->IsSolvent() ) continue;
		if( stricmp_loc(pmol->GetObjName(),"SOLVENT") == 0 ) continue;
		n_solute++;
		xc += aptr->GetX();
		yc += aptr->GetY();
		zc += aptr->GetZ();
	}

	if( n_solute == 0)
	{
		PrintLog("Error in MolEditor::CenterSoluteInSolvent() \n");
		PrintLog("No Solute Atoms found \n");
		return FALSE;
	}

	xc = xc/n_solute;
	yc = yc/n_solute;
	zc = zc/n_solute;

	double a = pmset->per_bc->GetA();
	double b = pmset->per_bc->GetB();
	double c = pmset->per_bc->GetC();

	Vec3D tr;

	int i;
	for( i = 0; i < 3; i++)
	{
		tr[i] =   0.5*pmset->per_bc->ucell[0][i]; 
		tr[i] +=  0.5*pmset->per_bc->ucell[1][i];
		tr[i] +=  0.5*pmset->per_bc->ucell[2][i];
	}

	tr[0] = tr[0] - xc;
	tr[1] = tr[1] - xc;
	tr[2] = tr[2] - xc;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->SetCoordFrom((Vec3D)(*aptr) + tr );
	}

	this->WrapToUnitCell(pmset,pmset->per_bc);

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->InitialTransform();
	    pview->DefaultRepresentation();
	}
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );		

	return TRUE;
}

int MolEditor::CenterMolInPBox(MolSet* pmset)
{
	if(!pmset->per_bc->IsSet())
	{
		PrintLog("Error in MolEditor::CenterMolInPBox() \n");
		PrintLog("Periodicity information is not set \n");
		return FALSE;
	}

	int na = pmset->GetNAtoms();
	if( na == 0 ) return FALSE;
	
	AtomIteratorGen aitr(pmset);
	HaAtom* aptr;

	Vec3D cnt(0.0,0.0,0.0);

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		cnt += (*aptr);
	}

	int i;
	for(i=0; i < 3; i++)
		cnt[i] = cnt[i]/na;

	HaVec_double trans(3);
	HaMat_double tr_mat(3,3,0.0);
	tr_mat.r0(0,0) = 1.0;
	tr_mat.r0(1,1) = 1.0;
	tr_mat.r0(2,2) = 1.0;

	trans[0] = pmset->per_bc->GetA()/2 - cnt[0];
	trans[1] = pmset->per_bc->GetB()/2 - cnt[1];
	trans[2] = pmset->per_bc->GetC()/2 - cnt[2];

	pmset->Transform(tr_mat,trans);

	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->InitialTransform();
	    pview->DefaultRepresentation();
	}
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );		
	return TRUE;
}

int MolEditor::MergeMolecules(HaMolecule* pMol1, HaMolecule* pMol2)
{
	if(pMol1 == pMol2) return TRUE;
	if(pMol1 == NULL || pMol2 == NULL)
	{
		PrintLog(" Error in MolEditor::MergeMolecules() \n");
		PrintLog(" One of the molecule pointer is NULL \n");
		return FALSE;
	}

	MolSet* pmset = pMol1->GetHostMolSet();
	
	pMol1->AddMolCopy(*pMol2,false,NULL);
	if(pMol2->GetHostMolSet() == pmset )
	{
		pmset->DeleteMol(pMol2);
	}
	return TRUE;
}

HaMolecule* MolEditor::CreateTransAlk(MolSet* pmset, const int nunit, const std::string& name)
{
	HaMolecule* pMol= pmset->AddNewMolecule();
	if(pMol == NULL) return NULL;
	
	pMol->SetObjName(name.c_str());
	HaChain* chain = pMol->AddChain(' ');
	
	int nres = 0; int atid = 0;
	
	int    lab_metyl[] ={ 6,1,1};
	int    lab_h = 1;
	
	double crd_methyl[3][3] = {
		{ -0.444560 , 1.886107 ,  0.000000 },
		{ -1.062324 , 1.886107 ,  0.873651 },
		{ -1.062324 , 1.886107 , -0.873651 },
	};
	
	double crd_h [3] = { -0.173206 , 2.759758 , 0.000000 } ;
	
	double t_ch = { 0.873651 };
	double t_cc = { 1.257406 };
	
	double y_h, y_c;
	int i,j;
	
	y_h = (t_cc * (nunit-1) + 2 * t_ch)/2.0;
	
	crd_h[1]=y_h ;

	HaAtom* aptr;
	HaAtom* aptr1;
	
	nres++;
	HaResidue* pres = chain->AddResidue(nres);
	pres->SetName("HTM");
	
	atid++;
	aptr=pres->AddNewAtom();
	aptr->SetElemNo(1);
	aptr->SetName("HT");
	aptr->SetX(crd_h[0] );
	aptr->SetY(crd_h[1] );
	aptr->SetZ(crd_h[2] );
	aptr1= aptr;

	y_c = y_h - t_ch;
	
	for (i = 1; i<= nunit; i++) {
		
		for(j = 0; j < 3; j++) {
			crd_methyl[j][1] = y_c;
			crd_methyl[j][0] = -crd_methyl[j][0] ; }
		
		crd_methyl[1][2] = -crd_methyl[1][2] ;
		crd_methyl[2][2] = -crd_methyl[2][2] ;
		
		nres++;
		pres = chain->AddResidue(nres);
		pres->SetName("MTL");
		
		atid++;
		aptr= pres->AddNewAtom();
		aptr->SetElemNo(6);
		aptr->SetName("CT");
		aptr->SetX(crd_methyl[0][0] );
		aptr->SetY(crd_methyl[0][1] );
		aptr->SetZ(crd_methyl[0][2] );
		pmset->AddBond(aptr,aptr1);
			
		aptr1= aptr;
		

		atid++;
		aptr=pres->AddNewAtom();
		aptr->SetElemNo(1);
		aptr->SetName("HT1");
		aptr->SetX(crd_methyl[1][0] );
		aptr->SetY(crd_methyl[1][1] );
		aptr->SetZ(crd_methyl[1][2] );
		pmset->AddBond(aptr,aptr1);
		
		atid++;
		aptr=pres->AddNewAtom();
		aptr->SetElemNo(1);
		aptr->SetName("HT2");
		aptr->SetX(crd_methyl[2][0] );
		aptr->SetY(crd_methyl[2][1] );
		aptr->SetZ(crd_methyl[2][2] );
		pmset->AddBond(aptr,aptr1);
		
		if( i != nunit) y_c=y_c - t_cc;
	}
	
    y_h=y_c - t_ch;

    crd_h[1] = y_h;
	
    if( nunit % 2 == 0)
		crd_h[0]=-crd_h[0];
	
	nres++;
	pres = chain->AddResidue(nres);
	pres->SetName("HTM");
	
	atid++;
	aptr=pres->AddNewAtom();

	aptr->SetElemNo(1);
	aptr->SetName("HT");
	aptr->SetX(crd_h[0] );
	aptr->SetY(crd_h[1] );
	aptr->SetZ(crd_h[2] );
	pmset->AddBond(aptr,aptr1);
	
	return pMol;
}

bool MolEditor::Create2DMolArray( MolSet* pmset, HaMolecule* pMol_ref,
	        					 const double deltx, const double delty, const int nx, const int ny,
						         const double alpha, const double tilt)
// Create 2D Array of Molecular Chains
// steps along symmetry vectors deltx and delty assumed to be in Ang
// alpha- angle between them is in Radians
// tilt - tilt angle of the molecules is in Radians
{
	int i,j;

	AtomIteratorMolecule aitr_ref(pMol_ref);

	HaAtom* aptr= aitr_ref.GetFirstAtom();

	if(aptr == NULL)
		return false;

	double xrot = aptr->GetX();
	double yrot = aptr->GetY();
	double zrot = aptr->GetZ();

	for(aptr= aitr_ref.GetFirstAtom();aptr; aptr = aitr_ref.GetNextAtom() )
	{
		if(aptr->GetY() < yrot)
		{
			xrot = aptr->GetX();
			yrot = aptr->GetY();
			zrot = aptr->GetZ();
		}
	}

	for(aptr= aitr_ref.GetFirstAtom();aptr; aptr = aitr_ref.GetNextAtom() )
	{
		double x = aptr->GetX();
		double y = aptr->GetY();
		
		aptr->SetX( x*cos(tilt) - y*sin(tilt));
		aptr->SetY( y*cos(tilt) + x*sin(tilt));
	}

	double cos_tilt = cos(tilt);
	if( cos_tilt == 0.0 ) cos_tilt=0.01;

	double dx1= ( deltx*sin(alpha)*cos(tilt) + delty*sin(tilt) );
	double dz1 = deltx* cos(alpha);
	double dz2 = delty;
//	double dy1 = dz2*sin(tilt);

	double dy1 = 0.0;

	for( i= 0; i < nx; i++)
	{
		for( j = 0; j < ny; j++)
		{	
			if(i == 0 && j == 0)
				continue;
			std::string mol_name= pmset->GetUniqueMolName(pMol_ref->GetObjName());
			HaMolecule* pMol= pmset->AddNewMolecule();
			pMol->AddMolCopy(*pMol_ref);
			pMol->SetObjName(mol_name.c_str());

			AtomIteratorMolecule aitr(pMol);
			
			for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				aptr->SetX( aptr->GetX() + dx1 * i );
				aptr->SetY( aptr->GetY() + dy1 * i );
				aptr->SetZ( aptr->GetZ() + dz1 * i + dz2 * j );
			}
		}
	}
	return true;
}


bool MolEditor::AddElectrSurf( MolSet* pmset,int add_surf_below_flag, int add_surf_top_flag, int add_atom_top_flag,
					    	 int add_atom_below_flag)
{
	double MinX, MinY, MinZ;
	double MaxX, MaxY, MaxZ;

	pmset->GetMinMaxCrd(MinX, MinY, MinZ, MaxX, MaxY, MaxZ);

	MinY-= 3.0;
	MaxY+= 3.0;

	double x0 = (MinX + MaxX)/2.0;
	double z0 = (MinZ + MaxZ)/2.0;

	double dlat = 1.5;

	double dx1 = dlat;
	double dz1 = dlat*0.5;
	double dz2 = dlat*sin(PI/3.0);

	int imax= (int) (( MaxX - MinX + 5.0)/dlat);
	int jmax= imax;

	int i,j;

	int nres=0;
	int atid=0;
	char atn[5];

	HaMolecule* pMol;
	HaChain* chain;
	HaResidue* pres;

	while(add_surf_below_flag || add_surf_top_flag )
	{
		pMol= pmset->AddNewMolecule();
		if(pMol == NULL) return false;

		if(add_surf_below_flag)
			pMol->SetObjName("SURF_LOW");
		else if( add_surf_top_flag)
			pMol->SetObjName("SURF_TOP");

		chain = pMol->AddChain(' ');

		int ic=-1;
		for(i= - imax; i <= imax; i++)
		{
			for(j = -jmax; j <= jmax; j++)
			{
				ic++;
				if(ic%10 == 0)
				{
					ic=0;
					nres++;
					pres = chain->AddResidue(nres);
					pres->SetName("SBE");		
				}
				HaAtom* aptr= pres->AddNewAtom();
				atn[0]='B';
				atn[1]='e';
				atn[2]='0'+ic;
				atn[3]=' ';
				atn[4] = 0;
				aptr->SetElemNo(4);
				aptr->SetName(atn);
				aptr->SetX(x0 + dx1 * i );

				if(add_surf_below_flag)
					aptr->SetY( MinY );
				else if(add_surf_top_flag)
					aptr->SetY( MaxY);
				aptr->SetZ(z0 + dz1 * i + dz2 * j);
			}
		}

		if(add_surf_below_flag)
			add_surf_below_flag = 0;
		else if(add_surf_top_flag)
			add_surf_top_flag = 0;
	}

	if(pMol != NULL) CreateCovBonds(pMol);

	if( add_atom_top_flag )
	{
		pMol= pmset->AddNewMolecule();
		if(pMol == NULL) return false;

		pMol->SetObjName("ATOM_TOP");
		chain = pMol->AddChain(' ');
		pres = chain->AddResidue(1);
		pres->SetName("DBE");		

		HaAtom* aptr= pres->AddNewAtom();

		atn[0]='B';
		atn[1]='e';
		atn[2]='1';
		atn[3]=' ';
		atn[4] =0;

		aptr->SetElemNo(4);
		aptr->SetName(atn);
		aptr->SetX( x0 );
		aptr->SetY( MaxY );
		aptr->SetZ( z0 );
	}

	if( add_atom_below_flag )
	{
		pMol= pmset->AddNewMolecule();
		if(pMol == NULL) return false;

		pMol->SetObjName("ATOM_BOTTOM");
		chain = pMol->AddChain(' ');
		pres = chain->AddResidue(1);
		pres->SetName("ABE");		

	    HaAtom* aptr= pres->AddNewAtom();

		atn[0]='B';
		atn[1]='e';
		atn[2]='1';
		atn[3]=' ';
		atn[4] =0;

		aptr->SetElemNo(4);
		aptr->SetName(atn);
		aptr->SetX( x0 );
		aptr->SetY( MinY );
		aptr->SetZ( z0 );
	}
	return true;
}

HaMolecule* MolEditor::CreateSurf(MolSet* pmset, const int num_layers, const std::string& name)
{
	if( pmset == NULL) return NULL;
	HaMolecule* pMol= pmset->AddNewMolecule();
	
	
	if(pMol == NULL) return NULL;
	
	pMol->SetObjName(name.c_str());
	HaChain* chain_cur = pMol->AddChain(' ');
	HaResidue* pres_cur = NULL;
	
	int nres = 0; int atid = 0;

	HaAtom* aptr;
	
	int ilayer;
	double dist = 5.44958;

	for(ilayer=0; ilayer < num_layers; ilayer++)
	{	
		int k,l;
		for( k = -10; k <= 10; k++)
		{
			for( l = -10; l <= 10; l++)
			{
				nres++;
				pres_cur = chain_cur->AddResidue(nres);
				pres_cur->SetName("AU");
				
				atid++;
				aptr=pres_cur->AddNewAtom();
				aptr->SetElemNo(79);
				aptr->SetName("AU");
			
				//double xx = dist*k;
				//double zz = dist*l;
				double xx = dist*k+dist*l/2;
				double zz = dist*l*1.73205/2;
				double yy = dist*ilayer;

				aptr->SetX(xx);
				aptr->SetY(yy);
				aptr->SetZ(zz);
			}
		}
	}
	
	return pMol;
}
