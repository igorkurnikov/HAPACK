//  hascattermod.cpp 
//
//  HARLEM
//
//  Classes 
//  to model electron scattering 
//
//  implementation
//
//  Igor Kurnikov , University of Pittsburgh 
//
//  Created:   August 28 2000
//

#define HASCATTERMOD_CPP

#include <math.h>
#include <mpi.h>

#include "tokens.h"
#include "hasurface.h"
#include "hascattermod.h"
#include "hamolset.h"
#include "haqchem.h"
#include "hamultipole.h"

#ifdef USE_IPACK  // IPACK headers

#include "errorip.h"
#include "const.h"
#include "basis.h"
#include "operators.h"
#include "parallel.h"
#include "f.h"
#endif 


HaScatterMod::HaScatterMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_SCATTER,new_phost_mset)
{
	pot_grid = new HaNonLocField3D;
	SetStdParams();
}

HaScatterMod::~HaScatterMod()
{
	if(pot_grid) delete pot_grid;
}

void
HaScatterMod::SetStdParams()
{
	if(pot_grid) pot_grid->SetDimensions(10,10,10);

	MolSet* pmset= GetMolSet();
	if(pmset == NULL)
		return;

	double xmin, ymin, zmin;
	double xmax, ymax, zmax;

	pmset->GetMinMaxCrd(xmin,ymin,zmin,xmax,ymax,zmax);
	
	xmin = xmin - 1.0;
	ymin = ymin - 1.0;
	zmin = zmin - 1.0;

	xmax = xmax + 1.0;
	ymax = ymax + 1.0;
	zmax = zmax + 1.0;

	if(pot_grid != NULL)
	{
		pot_grid->SetGridCornersCoord(xmin,ymin,zmin,xmax,ymax,zmax);
	}

	psp_thresh_val = 1.0e-2;
	pt_exp = 1.0;
	psp_fname = "psp_grid.dat";
}

int
HaScatterMod::CalcPseudoPotGrid()
{
	int ires = CalcPseudoPotAO();

	if( ires == FALSE)
		return FALSE;

	MolSet* pmset = this->GetMolSet();
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(false);
	
	if(pot_grid == NULL)
	{
		ErrorInMod("HaScatterMod::CalcPseudoPotGrid()",
			       " Grid for Pseudopotential is not set ");
		return FALSE;
	}

	int nb = ptr_qc_mod->GetNBfunc();

	int nx = pot_grid->GetNx();
	int ny = pot_grid->GetNy();
	int nz = pot_grid->GetNz();

	if( nx < 2 || ny < 2 || nz < 2)
	{
		ErrorInMod("HaScatterMod::CalcPseudoPotGrid()",
			       " Some of the grid dimensions are less than 2");
		return FALSE;
	}

	double xstep = pot_grid->GetXstep();
	double ystep = pot_grid->GetYstep();
	double zstep = pot_grid->GetZstep();

	int ix,iy,iz;

	int ngrid_size = nx * ny * nz;

	pot_grid->fvals.resize(ngrid_size);
	
	HaMat_double bf_grid_val; // basis functions values in the grid points
	
	bf_grid_val.newsize(nb,ngrid_size);
	bf_grid_val = 0.0;

	PrintLog("Building an array of values of basis functions on the grid \n");

	int bf_idx_glob = 0;
	AtBasisType::iterator atb_itr;
	for(atb_itr= ptr_qc_mod->AtBasis.at_bas_vec.begin(); 
	    atb_itr != ptr_qc_mod->AtBasis.at_bas_vec.end(); atb_itr++)
	{
		ShellsType::iterator shitr;
		const HaAtom* aptr= (*atb_itr).GetAtHost();

		for(shitr = (*atb_itr).Shells.begin(); 
		    shitr != (*atb_itr).Shells.end(); shitr++)
		{
			int nf= (*shitr).GetNBfunc();
			int ifun;

			for( ifun = 1; ifun <= nf; ifun++)
			{
				bf_idx_glob++;
				HaVec_double cf_axx(nf);
				cf_axx = 0.0;
				cf_axx(ifun) = 1.0;

				int grid_lin_idx = 0;
				for(iz = 0; iz < nz; iz++)
				{
					for(iy = 0; iy < ny; iy++)
					{
						for(ix = 0; ix < nx; ix++)
						{
							grid_lin_idx++;
							float x,y,z;
							pot_grid->GetXYZ(x,y,z,ix,iy,iz);
							double dx, dy,dz;
							dx = x - aptr->GetX(); dy = y - aptr->GetY(); dz = z - aptr->GetZ();
							double ff = 0.0;
							ff+= 0.25* (*shitr).EvalLinCombInPoint( dx, dy, dz, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx + xstep/2.0, dy, dz, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx - xstep/2.0, dy, dz, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx, dy + ystep/2.0, dz, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx, dy - ystep/2.0, dz, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx, dy, dz + zstep/2.0, cf_axx.begin() );
							ff+= 0.125* (*shitr).EvalLinCombInPoint( dx, dy, dz - zstep/2.0, cf_axx.begin() );
							
							bf_grid_val(bf_idx_glob,grid_lin_idx) = ff;
						}
					}
				}
			}
		}
	}

//	PrintLog("Find  interactions between grid points \n");


	int ir2,ic2,il2,lidx2;

	PrintLog("Find  interactions between grid points \n");


	int lidx1 = 0;
	for(iz=1; iz <= nz; iz++)
	  for(iy =1 ; iy <= ny; iy++)
	    for(ix = 1; ix <= nx; ix++)
		{
		   lidx1++;
		   lidx2 = 0;
		   HaMat_double c1;
		   c1.set_ext_alloc(&bf_grid_val(1,lidx1),nb,1);
		   for(il2 = 1; il2 <= nz; il2++)
			 for(ic2 = 1; ic2 <= ny; ic2++)
				for(ir2 = 1; ir2 <= nx; ir2++ )
				{
					lidx2++;
					HaMat_double c2;
					c2.set_ext_alloc(&bf_grid_val(1,lidx2),nb,1);
					HaMat_double scr, val;

					matmult_T1(scr,c1,pseudo_pot_ao);
					matmult(val,scr,c2);

					if(fabs(val(1,1)) > psp_thresh_val)
					{
						pot_grid->fvals[lidx1-1].push_back(ValAtPoint(ir2,ic2,il2,val(1,1)) );
					}
				}

		}
	return True;
}

int
HaScatterMod::CalcPseudoPotAO()
{
	MolSet* pmset = this->GetMolSet();
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(false);

	if(ptr_qc_mod == NULL)
	{
		ErrorInMod("HaScatterMod::CalcPseudoPotAO()",
			       " QChem Module is not set ");
		return FALSE;
	}

	int nb = ptr_qc_mod->GetNBfunc();

	HaMat_double fock_mat, t_kin_mat, ss_m1;
	double ene_cut = -5.0; // cutoff ennergy 
	ptr_qc_mod->BuildFockMatFromMOs(fock_mat,ene_cut);

	HaOperKinEner T_oper;
	
	ss_m1 = ptr_qc_mod->GetOvlpMat();
	HaMat_double::mat_inverse(ss_m1);

	T_oper.EvalGauBasisSet(&(ptr_qc_mod->AtBasis),t_kin_mat);

	pseudo_pot_ao.newsize(nb,nb);

	mat_diff(pseudo_pot_ao,fock_mat,t_kin_mat);

	cout << " Pseudo-potential matrix on AOs before multiplying by S-1 " << endl; 
	pseudo_pot_ao.Print_format(cout,"%7.4f ");

	t_kin_mat.newsize(0,0);
	
	HaMat_double scr;

	matmult(scr,ss_m1,pseudo_pot_ao);
	matmult(pseudo_pot_ao,scr,ss_m1);

	return True;
}

double
HaScatterMod::GetPsPGridElem(  int nx_1, int ny_1, int nz_1, 
								int nx_2, int ny_2, int nz_2)
{
	if( pot_grid == NULL)
		return 0.0;

	if( nx_1 < 1 || ny_1 < 1 || nz_1 < 1 || 
		nx_2 < 1 || ny_2 < 1 || nz_2 < 1 )
	{
		return 0.0;
	}

	int nx_max = pot_grid->GetNx();
	int ny_max = pot_grid->GetNy();
	int nz_max = pot_grid->GetNz();

	if( nx_1 > nx_max ||  nx_2 > nx_max ||
		ny_1 > ny_max ||  ny_2 > ny_max ||
		nz_1 > nx_max ||  nz_2 > nz_max )
	{
		return 0.0;
	}
	
	float xf,yf,zf;

	double x1, y1, z1, x2, y2, z2;

	bool bres = pot_grid->GetXYZ(xf,yf,zf,nx_1,ny_1,nz_1);

	x1 = xf; y1 = yf ; z1 = zf;

	bres = pot_grid->GetXYZ(xf,yf,zf,nx_1,ny_1,nz_1);

	x2 = xf; y2 = yf ; z2 = zf;

	return 0.0;
}

double
HaScatterMod::GetPsP_xyz(double x1, double y1, double z1,
						 double x2, double y2, double z2)
{
	MolSet* pmset = GetMolSet();
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(false);

	if(ptr_qc_mod == NULL)
		return 0.0;

	int nb = ptr_qc_mod->GetNBfunc();

	HaMat_double c1(nb,1,0.0), c2(nb,1,0.0);
	int idx_glob = 0;

	AtBasisType::iterator atb_itr;
	for(atb_itr= ptr_qc_mod->AtBasis.at_bas_vec.begin(); 
	    atb_itr != ptr_qc_mod->AtBasis.at_bas_vec.end(); atb_itr++)
	{
		ShellsType::iterator shitr;
		const HaAtom* aptr= (*atb_itr).GetAtHost();

		for(shitr = (*atb_itr).Shells.begin(); 
		    shitr != (*atb_itr).Shells.end(); shitr++)
		{
			int nf= (*shitr).GetNBfunc();
			int ifun;

			for( ifun = 1; ifun <= nf; ifun++)
			{
				idx_glob++;
				HaVec_double cf_axx(nf);
				cf_axx = 0.0;
				cf_axx(ifun) = 1.0;
				double dx, dy,dz;
				dx = x1 - aptr->GetX(); dy = y1 - aptr->GetY(); dz = z1 - aptr->GetZ();
				c1(idx_glob,1) = (*shitr).EvalLinCombInPoint( dx, dy, dz, cf_axx.begin() );
				dx = x2 - aptr->GetX(); dy = y2 - aptr->GetY(); dz = z2 - aptr->GetZ();
				c2(idx_glob,1) = (*shitr).EvalLinCombInPoint( dx, dy, dz, cf_axx.begin() );
			}
		}
	}
	
	if(pseudo_pot_ao.num_rows() != nb || pseudo_pot_ao.num_cols() != nb)
		CalcPseudoPotAO();

	if(pseudo_pot_ao.num_rows() != nb || pseudo_pot_ao.num_cols() != nb)
		return 0.0;

	
	HaMat_double scr, val;

	matmult_T1(scr,c1,pseudo_pot_ao);
	matmult(val,scr,c2);

	return (val(1,1));
}


double
HaScatterMod::GetKinEne_Gauss_xyz(double x1, double y1, double z1,
	                           double x2, double y2, double z2)
{	
	HaQCMod::InitIPack();

	PureBasisSet grid_pure_bf_set;
	PureBasisFunction* pure_grid_bfun =  new PureBasisFunction(0);
	pure_grid_bfun->add(pt_exp,1.0);
	grid_pure_bf_set.add(pure_grid_bfun);

	BasisSet grid_bset;
	String str;
	Location pt1(x1,y1,z1);
	Location pt2(x2,y2,z2);
	
	grid_bset.add(&grid_pure_bf_set,str,pt1);
	grid_bset.add(&grid_pure_bf_set,str,pt2);

	Location dum_loc;
	OrbInfo  orbinfo;

	InternalBasis* grid_int_bas = new InternalBasis(grid_bset, dum_loc, ALL, orbinfo);

	ARRAY<IntegMat> ipack_kin_ene_matr;
	ipack_kin_ene_matr.reset(1);
	ipack_kin_ene_matr[0].reset(2,2);
	ipack_kin_ene_matr[0].set(0);

	compute_kinetic(*grid_int_bas,*grid_int_bas,ipack_kin_ene_matr,TRUE);

	return( ipack_kin_ene_matr[0](0,1) );

}

double
HaScatterMod::GetPsP_Gauss_xyz(double x1, double y1, double z1,
	                           double x2, double y2, double z2)
{
#ifdef USE_IPACK	
	MolSet* pmset = GetMolSet();
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(false);

	if(ptr_qc_mod == NULL)
		return 0.0;
	
	int nb = ptr_qc_mod->GetNBfunc();

	if(pseudo_pot_ao.num_rows() != nb || pseudo_pot_ao.num_cols() != nb)
		CalcPseudoPotAO();

	HaQCMod::InitIPack();

	InternalBasis* mol_bas = ptr_qc_mod->AtBasis.CreateIPackBas();

	PureBasisSet grid_pure_bf_set;
	PureBasisFunction* pure_grid_bfun =  new PureBasisFunction(0);
	pure_grid_bfun->add(pt_exp,1.0);
	grid_pure_bf_set.add(pure_grid_bfun);

	BasisSet grid_bset;
	String str;
	Location pt1(x1,y1,z1);
	Location pt2(x2,y2,z2);
	
	grid_bset.add(&grid_pure_bf_set,str,pt1);
	grid_bset.add(&grid_pure_bf_set,str,pt2);

	Location dum_loc;
	OrbInfo  orbinfo;

	InternalBasis* grid_int_bas = new InternalBasis(grid_bset, dum_loc, ALL, orbinfo);

	ARRAY<IntegMat> ipack_ovlp_matr;
	ipack_ovlp_matr.reset(1);
	ipack_ovlp_matr[0].reset(2,nb);
	ipack_ovlp_matr[0].set(0);

	compute_overlap(*grid_int_bas,*mol_bas,ipack_ovlp_matr,TRUE);

	HaMat_double s1(nb,1), s2(nb,1);

	cout << " dimensions of the overlap matrix "<< ipack_ovlp_matr[0].size1() << "  " << ipack_ovlp_matr[0].size2() << endl;

	int i;
	for(i = 0; i < nb; i++) 
	{

		s1(i+1,1) = ipack_ovlp_matr[0](0,i);
		s2(i+1,1) = ipack_ovlp_matr[0](1,i);
	}

	HaMat_double scr, psp_grid;
	
	matmult_T1(scr,s1,pseudo_pot_ao);
	matmult(psp_grid,scr,s2);

	return psp_grid(1,1);

#endif
	return 0.0;
}
int
HaScatterMod::FindGridHamEigVec(int num_vec)
{
	if(pot_grid == NULL)
		return FALSE;

	int nx = pot_grid->GetNx();
	int ny = pot_grid->GetNy();
	int nz = pot_grid->GetNz();

	if( nx < 2 || nz < 2 || ny < 2 )
	{
		ErrorInMod(" HaScatterMod::FindGridHamEigVec()",
			       " Some grid diemensions less than 2");
		return FALSE;
	}
	
	int ngrid_size = nx*ny*nz;

	double dx = pot_grid->GetXstep();
	double dy = pot_grid->GetYstep();
	double dz = pot_grid->GetZstep();

	if( dx < 0.0001 || dy < 0.0001 || dz < 0.0001)
	{
		ErrorInMod(" HaScatterMod::FindGridHamEigVec()",
			       " Some grid steps are less than 0.0001 Bohr ");
		return FALSE;
	}

	double dV= dx*dy*dz;

	if( ngrid_size == 0)
	{
		ErrorInMod(" HaScatterMod::FindGridHamEigVec()",
			       " Grid Size in zero ");
		return FALSE;
	}

	if( ngrid_size > 1500 )
	{
		ErrorInMod(" HaScatterMod::FindGridHamEigVec()",
			       " Grid Size is > 1500 , too large for this function");
		return FALSE;
	}

	HaMat_double ham(ngrid_size,ngrid_size,0.0);

	int fval_size = pot_grid->fvals.size();

	if(fval_size != ngrid_size)
	{
		ErrorInMod(" HaScatterMod::FindGridHamEigVec()",
			       " Grid Size is not equal to the dimension of pseudopot values array");
		return FALSE;
	}

	list<ValAtPoint>::iterator vitr;

	int i;
	for( i = 1; i <= ngrid_size ; i++)
	{
		list<ValAtPoint>& vlist = pot_grid->fvals[i-1];
		for( vitr = vlist.begin(); vitr != vlist.end(); vitr++)
		{
			int lidx = pot_grid->GetLinIdx((*vitr).ix,(*vitr).iy,(*vitr).iz);
			ham(i,lidx) = (*vitr).val;

//			cout << " i = " << i << ";   lidx = " <<  lidx <<  endl;
//			ham(i,lidx) = -1.0;
		}
	}

    int ix, iy, iz;
		
	double rdx2 = 1.0/(dx*dx);
	double rdy2 = 1.0/(dy*dy);
	double rdz2 = 1.0/(dz*dz);

//	ham = 0.0;

	for( iz = 0; iz < nz; iz++ )
	{
		for( iy = 0; iy < ny; iy++ )
		{
			for( ix = 0; ix < nx; ix++ )
			{
				int lidx_0 = pot_grid->GetLinIdx(ix,iy,iz); 
				ham(lidx_0,lidx_0)+= rdx2 + rdy2 + rdz2;
				int lidx_1,lidx_2;
				if( ix != 0)
				{
					lidx_1 = pot_grid->GetLinIdx(ix-1,iy,iz);
					ham(lidx_1,lidx_0)-= 0.25* rdx2;
					ham(lidx_0,lidx_1)-= 0.25* rdx2;
				}
				if( ix != nx-1 )
				{
					lidx_2 = pot_grid->GetLinIdx(ix+1,iy,iz);
					ham(lidx_2,lidx_0)-= 0.25* rdx2;
					ham(lidx_0,lidx_2)-= 0.25* rdx2;
				}
				if( iy != 0 )
				{
					lidx_1 = pot_grid->GetLinIdx(ix,iy-1,iz);
					ham(lidx_1,lidx_0)-= 0.25* rdy2;
					ham(lidx_0,lidx_1)-= 0.25* rdy2;
				}
				if( iy != ny-1 )
				{
					lidx_2 = pot_grid->GetLinIdx(ix,iy+1,iz);
					ham(lidx_2,lidx_0)-= 0.25* rdy2;
					ham(lidx_0,lidx_2)-= 0.25* rdy2;
				}
				if( iz != 0 )
				{
					lidx_1 = pot_grid->GetLinIdx(ix,iy,iz-1);
					ham(lidx_1,lidx_0)-= 0.25* rdz2;
					ham(lidx_0,lidx_1)-= 0.25* rdz2;
				}
				if( iz != nz-1 )
				{
					lidx_2 = pot_grid->GetLinIdx(ix,iy,iz+1);
					ham(lidx_2,lidx_0)-= 0.25* rdz2;
					ham(lidx_0,lidx_2)-= 0.25* rdz2;
				}
			}
		}
	}
	
	HaVec_double eig_val;
	HaMat_double cc(ngrid_size,ngrid_size);

	HaMat_double::mat_sdiag(ham,cc,eig_val);
	PrintLog("Lowest Eigenvalues of the grid hamiiltonian: \n");

	for(int j = 1; j <= num_vec; j++)
	{
		PrintLog(" E(%3d)= %12.6f \n",j,eig_val(j) );
	}
	
	return TRUE;
}


int 
HaScatterMod::SetCoreHamOnGrid()
{
	MolSet* pmset = GetMolSet();	
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(true);

	int cntr_num = ptr_qc_mod->GetNumCnt();
	HaMat_double cntr_coord(3,cntr_num);
	HaVec_double cntr_charges(cntr_num);

	ptr_qc_mod->GetCntCharges(cntr_charges);
	ptr_qc_mod->AtBasis.GetCntCoord(cntr_coord);

	if( pot_grid == NULL)
		return FALSE;

	int nx = pot_grid->GetNx();
	int ny = pot_grid->GetNy();
	int nz = pot_grid->GetNz();

	if(nx < 2)
		return FALSE;

	double xstep = pot_grid->GetXstep();

	if(xstep < 0.001)
		return FALSE;

	pot_grid->fvals.resize(nx*ny*nz);

	int ix,iy,iz;
	int lidx1 = 0;
	for(iz=0; iz < nz; iz++)
	{
	  for(iy =0 ; iy < ny; iy++)
	  {
	    for(ix = 0; ix < nx; ix++)
		{
		   lidx1++;
		   pot_grid->fvals[lidx1-1].clear();
		   float xg,yg,zg;
		   pot_grid->GetXYZ(xg,yg,zg,ix,iy,iz);

		   double fv = 0.0;
		   int i;
		   for(i = 1; i <= cntr_num; i++)
		   {
			   double dd =0.0;
			   dd += (xg - cntr_coord(1,i))* (xg - cntr_coord(1,i));
			   dd += (yg - cntr_coord(2,i))* (yg - cntr_coord(2,i));
			   dd += (zg - cntr_coord(3,i))* (zg - cntr_coord(3,i));
			   dd = sqrt(dd);

			   if( dd < 0.01)
			   {
				   fv+=  -cntr_charges(i)* (3.0/xstep);
			   }
			   else
				   fv+= -cntr_charges(i)/dd;
		   }
		   pot_grid->fvals[lidx1-1].push_back(ValAtPoint(ix,iy,iz,fv) );
		}
	  }
	}
	return TRUE;
}
