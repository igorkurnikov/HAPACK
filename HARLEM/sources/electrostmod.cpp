/*! electrostmod.cpp

    Classes to perform Continuum Elctrostatics Calculations in HARLEM  

    \author Igor Kurnikov 
    \date   1999-2002
*/

#define ELECTROSTMOD_CPP

#include <mpi.h>

#if defined(_MSC_VER)
#include <process.h>
#endif

#include "command.h"
#include "stdlib.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include "electrostmod.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "etcoupl.h"
#include "hasurface.h"
#include "hamolview.h"
#include "harlemapp.h"
#include "protonredox.h"

//#include "cworld.h"

int ElectrostMod::ActiveElectrMod = COMP_MOD_ELECTROST;

ElectrostMod::ElectrostMod(MolSet* new_phost_mset,const int new_mtype) :
HaCompMod(new_mtype,new_phost_mset)
{
	SetStdParam();
}

ElectrostMod::~ElectrostMod()
{
	clear();
}


bool
ElectrostMod::SetStdParam()
{
	param_file_title = "! DELPHI version 3.0 qdiffxs parameter file" ;
	nx=ny=nz=65;
	perfil= 60.0;
	offset[0]=offset[1]=offset[2]=0.0;
	epsi=4.0; epsout = 81.0;
	rionst = 0.0;
	exrad = 2.0; radprb= 1.4;
	boundary = 0;
	iper[0]= iper[1]= iper[2]= 0;
	nlit  = -1;
	nnit  = 0;
	iconc = 0; ibios = 1;
    isite = 0;
	iatout = 1;
	toplbl = "Sample Delphi calculations";
	isph = 0;
	ipdbwrt = 0;
	ifrcwrt = 0;
	enc = "GCSA";
	igraph=0; ipotent=0; icon1=10; icon2=1;
	imem = 0;
	phiwrt = 0;
	ihs=0,isen=0,ish=0;

	pot_isolevel = 2.0;
	dots_number = 100;
	elpot_low_val = -2.0;
	elpot_high_val = 2.0;

	rdx_shift = 0.0;
	rdx_shft_mode = RDX_SHFT_VS_VAC;

	fill_charge_mode = FILL_CHARGE_PROD_REV_DIST;

	elfield_fname = "elfield.dat";

	ClearBoundaryAtoms();

	if(phost_mset != NULL)
	{
		std::string mset_name = phost_mset->GetName();
		param_file_name  = mset_name + "_delph.prm";
		radius_file_name = mset_name + "_delph.siz";
		charge_file_name = mset_name + "_delph.crg";
		coord_file_name  = mset_name + "_delph_inp.pdb";
		log_file_name    = mset_name + "_delph.log";
	}
	else
	{	
		param_file_name =  "DELPHI_INP.prm";
		radius_file_name = "DELPHI_INP.siz";
		charge_file_name = "DELPHI_INP.crg";
		coord_file_name  = "DELPHI_INP.pdb";
		log_file_name    = "DELPHI.log";
	}

	return true;
}

void
ElectrostMod::clear()
{
	el_pot_map.clear();
}


#define BOOL_SYMB(x)  ( ((x) == 0)? "f" : "t" )  


bool ElectrostMod::SaveParamFile()
{
	char buf[256];
	std::string tail, tmp;
	std::ofstream fpar(param_file_name.c_str());
	
//	fpar << param_file_title << endl;

    tail="                 ! igrid - Grid size";
    sprintf(buf," %3d %3d %3d %s\n",nx,ny,nz,tail.c_str());
	fpar << buf;

    tail="                 ! perfil - dimensions of the box, or scale -? ";
	sprintf(buf,"  %6.3f %s %s %s \n",perfil,"j","j",tail.c_str());
	fpar << buf;

    tail="! offset - Molecule Center (x,y,z)";
    sprintf(buf,"%9.4f %9.4f %9.4f %s\n", offset[0], offset[1], offset[2], tail.c_str());
	fpar << buf;

    tail="! epsi, epsout - int. and ext. diel. constants";
    sprintf(buf,"%12.6f %12.6f     %s\n", epsi, epsout, tail.c_str());
	fpar << buf;

    tail="               ! rionst - ionic strength (M/l)";
    sprintf(buf,"%12.6f   %s\n", rionst, tail.c_str());
	fpar << buf;

    tail="    ! exrad, radprb - Ion excl. radius, probe rad.";
    sprintf(buf,"%12.6f %12.6f %s\n", exrad, radprb,tail.c_str());
	fpar << buf;

    tail="               ! Boundary: 1=Zero,2=Apprx Coul,3=Focus,4=Coul";
    sprintf(buf,"%2d             %s\n",(boundary+1),tail.c_str());
	fpar << buf;

    tail="                       ! x,y,z periodic boundary flags (t/f)";
    sprintf(buf," %s %s %s %s\n",
                     BOOL_SYMB(iper[0]),BOOL_SYMB(iper[1]),BOOL_SYMB(iper[2]),tail.c_str());
	fpar << buf;

    tail="                        ! nlit - lin. P-B iter. or a-auto,or gten";
	if( nlit < 0 )
	{
		sprintf(buf,"%s %s\n", "a    ", tail.c_str());
	}
	else
	{
		sprintf(buf,"%5d %s\n", nlit, tail.c_str());
	}
	fpar << buf;

    tail="                        ! nnit- # non-linear P-B iterations";
    sprintf(buf,"%5d %s\n", nnit, tail.c_str());
	fpar << buf;

    tail="                        ! iconc,ibios- pot. output control";
    sprintf(buf," %s %s  %s\n", BOOL_SYMB(iconc), BOOL_SYMB(ibios),tail.c_str() );
	fpar << buf;

    tail="                           ! isite - flag print el. fld at spec. points";
    sprintf(buf," %s %s\n",BOOL_SYMB(isite),tail.c_str());
	fpar << buf;

    tail="                           ! iatout - flag to to write modif. pdb ";
    sprintf(buf," %s %s\n",BOOL_SYMB(iatout),tail.c_str());
	fpar << buf;

	fpar << toplbl << endl;

    tail="                           ! isph - flag for spherical charge distrib";
    sprintf(buf," %s %s\n",BOOL_SYMB(isph),tail.c_str());
	fpar << buf;

    tail="                           ! ipdbwrt - flag for writing pf pdb file";
    sprintf(buf," %s %s\n",BOOL_SYMB(ipdbwrt),tail.c_str());
	fpar << buf;

    tail="                        ! ifrcwrt - flag for writing frc file";
    sprintf(buf," %s    %s\n",BOOL_SYMB(ifrcwrt),tail.c_str());
	fpar << buf;

    tail="                       ! enc -energy print control flags ";
    sprintf(buf," %s  %s\n", enc.c_str(),tail.c_str());
	fpar << buf;

    tail="                 ! igraph, ipotent, icon1, icon2 ";
    sprintf(buf," %s %s %3d %3d %s\n", BOOL_SYMB(igraph),BOOL_SYMB(ipotent),icon1,icon2,tail.c_str());
	fpar << buf;

    tail="                          ! imem - flag to model  membrane ";
    sprintf(buf," %s  %s\n", BOOL_SYMB(imem), tail.c_str());
	fpar << buf;

    tail="                          ! phiwrt - flag to write potential map ";
    sprintf(buf," %s  %s\n", BOOL_SYMB(phiwrt),tail.c_str());
	fpar << buf;

    tail="                     ! ihs,isen,isch - flags to write surface info ";
    sprintf(buf," %s %s %s   %s\n", BOOL_SYMB(ihs), BOOL_SYMB(isen),BOOL_SYMB(ish),tail.c_str());
	fpar << buf;

	fpar.close();

	return true;
}


bool ElectrostMod::SaveRadiusFile()
{
	map< std::string, double, less<std::string> > rad_map;
	map< std::string, double, less<std::string> >::iterator ritr;

	if(phost_mset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(phost_mset);
	    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			std::string label;
			std::string atname = aptr->GetName();
			label = atname;
			for(;;)
			{
				if(label.length() < 6 )
					label+= ' ';
				else
					break;
			}
			HaResidue* pres= aptr->GetHostRes(); 
			char res_name[5];
			strncpy(res_name,pres->GetName(),4);
			if(res_name[0]== ' ')
			{
				res_name[0] = res_name[1];
				res_name[1] = res_name[2];
				res_name[2] = res_name[3];
			}
			res_name[4]=0;
			label+= res_name;
            rad_map[label]= aptr->radius;
		}
	}

	FILE* frad = fopen(radius_file_name.c_str(),"w");
	fprintf(frad,"! Molecule atom radii \n");
	fprintf(frad,"atom__res_radius_ \n");

    for(ritr = rad_map.begin(); ritr != rad_map.end(); ritr++)
	{
		std::string label = (*ritr).first;
		double drad = (*ritr).second;
		fprintf(frad,"%s%8.3f \n",label.c_str(),drad);
	}
	fprintf(frad,"xx          0.00 \n");
	fclose(frad);

//	ofstream frad(radius_file_name.c_str());
//	frad << "! PARSE atom radii set" << endl;
//	frad << "atom__res_radius_" << endl;
//	frad << "c           1.70 " << endl;
//	frad << "o           1.60 " << endl;
//	frad << "n           1.50 " << endl;
//	frad << "h           1.00 " << endl;
//	frad << "s           1.90 " << endl;
//	frad << "p           2.00 " << endl;
//	frad << "ir1         2.00 " << endl;
//	frad << "ir2         2.00 " << endl;
//	frad << "fe          2.00 " << endl;
//	frad << "fe1         2.00 " << endl;
//	frad << "fe2         2.00 " << endl;
//	frad << "fe3         2.00 " << endl;
//	frad << "fe4         2.00 " << endl;
//	frad << "fe5         2.00 " << endl;
//	frad << "fe6         2.00 " << endl;
//	frad << "fe7         2.00 " << endl;
//	frad << "fe8         2.00 " << endl;
//	frad << "ru          2.00 " << endl;
//	frad << "mg          2.00 " << endl;
//	frad << "al          2.00 " << endl;
//	frad << "f           2.00 " << endl;
//	frad << "mo          2.00 " << endl;
//	frad << "be          1.50 " << endl;
//	frad << "xx          0.00 " << endl;
//
//	frad.close();
	return true;
}

bool
ElectrostMod::SaveChargeFile()
{
	char buf[256];
	ofstream fcrg(charge_file_name.c_str());
	fcrg << "! default charge file   " << endl;
	fcrg << "atom__resnumbc_charge_  " << endl;
	if(phost_mset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(phost_mset);
	    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			if( fabs( aptr->GetCharge() ) < DBL_EPSILON )
				continue;
			int j=0;
			j+= sprintf(buf+j,"%s",aptr->GetName());
			for(; j < 6 ;)
			{
				j+=sprintf(buf+j,"%s"," ");
			}
			HaResidue* pres= aptr->GetHostRes(); 
			char res_name[5];
			strncpy(res_name,pres->GetName(),4);
			if(res_name[0]== ' ')
			{
				res_name[0] = res_name[1];
				res_name[1] = res_name[2];
				res_name[2] = res_name[3];
			}
			res_name[4]=0;
			j+= sprintf(buf+j,"%s", res_name );
			for(; j <  9; )
			{
				j+=sprintf(buf+j,"%s"," ");
			}
			j+= sprintf(buf+j,"%4d", pres->GetSerNo() );
			HaChain* chain= aptr->GetHostChain();
			j+= sprintf(buf+j,"%c",chain->ident);
			j+= sprintf(buf+j,"%10.5f", aptr->charge);
			fcrg << buf << endl;
		}
	}
	fcrg.close();
	return true;
}

bool
ElectrostMod::SaveCoordFile()
{
	if(phost_mset == NULL)
		return false;
	phost_mset->SavePDBFile(coord_file_name.c_str());
    AddBoundaryAtoms();

	return true;
}


int
ElectrostMod::LoadElPotFromFile(int format)
{
	char fname_buf[256];
	
	int ires = FALSE;
	if(format == HLM_F3D_FORM) format = HLM_F3D_BIN;

	if( format == HLM_F3D_BIN )
	{
		ires = el_pot_map.RestoreFromFile(elfield_fname.c_str(),TRUE);
		if(ires)
		{
			nx = el_pot_map.GetNx();
			ny = el_pot_map.GetNy();
			nz = el_pot_map.GetNz();
		}
		return ires;
	}
	else if(format == 2)
	{
		int nx_f, ny_f, nz_f;
		freal   xmin_f, xmax_f, ymin_f, ymax_f, zmin_f, zmax_f;  
	
		strncpy(fname_buf,elfield_fname.c_str(),255);
		openpm_(&nx_f,&ny_f,&nz_f,&xmin_f,&xmax_f,&ymin_f,&ymax_f,&zmin_f,&zmax_f,
		          fname_buf, elfield_fname.size() );

		if (nx_f <= 0 || ny_f <= 0 || nz_f <= 0)
			return false;
	
		el_pot_map.SetGridCornersCoord(xmin_f,ymin_f,zmin_f,xmax_f,ymax_f,zmax_f); 

		el_pot_map.SetDimensions(nx_f, ny_f, nz_f);

		loadpm_(&nx_f, &ny_f, &nz_f, el_pot_map.GetFieldPtr() );
	}
	return TRUE;
}

bool ElectrostMod::BuildPotIsoSurface()
{
	if( el_pot_map.GetNx() <= 0 || el_pot_map.GetNy() <= 0 ||
		el_pot_map.GetNz() <= 0)
	{
		return false;
	}
	
	MolSet* pmset = GetMolSet();
	if( pmset == NULL)
		return false ;

	char buf[256];
	HaDisplayedSurface* sptr;

// Positive Potential Isosurface
	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return false;

	Surfaces.push_back(sptr);
	pmset->AddObject3D(sptr);
	
	bool result;
	result = sptr->calc_isosurf(&el_pot_map, pot_isolevel);
	sptr->ColourUniform(0,0,255);
	sprintf(buf,"isopot_surf_%f",pot_isolevel);
	sptr->SetObjName(buf);

// Negative Potential Isosurface

	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return false;

	Surfaces.push_back(sptr);
	pmset->AddObject3D(sptr);

	result = sptr->calc_isosurf(&el_pot_map, -pot_isolevel);
	sptr->ColourUniform(255,0,0);
	
	sprintf(buf,"isopot_surf_%f",-pot_isolevel);
	sptr->SetObjName(buf);


//	sptr->Print_info(cout,1);
	return result;
}
bool 
ElectrostMod::FillChargeMap()
{
	charge_map.SetGridCornersCoord( el_pot_map.GetXmin(), el_pot_map.GetYmin(),el_pot_map.GetZmin(),
		                            el_pot_map.GetXmax(), el_pot_map.GetYmax(),el_pot_map.GetZmax()); 

	charge_map.SetDimensions(el_pot_map.GetNx(), el_pot_map.GetNy(),el_pot_map.GetNz());
    if(phost_mset == NULL)
		return false;   
	
	HaAtom* aptr;
	double xstep,ystep,zstep;
	xstep=charge_map.GetXstep();
	ystep=charge_map.GetYstep();
	zstep=charge_map.GetZstep();

    if( fabs(xstep*ystep*zstep) < 0.000000001 ) { return false; } 

	double el_vol = xstep*ystep*zstep; // volume element

	double x_min=charge_map.GetXmin();
	double y_min=charge_map.GetYmin();
	double z_min=charge_map.GetZmin();
	AtomIteratorMolSet aitr(phost_mset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		double at_charge;
		
		at_charge=aptr->GetCharge();		
		at_charge = at_charge/el_vol; // find charge density in the grid corresponding to this charge
		
		double x_at=aptr->GetX();
		double y_at=aptr->GetY();
		double z_at=aptr->GetZ();
		int ix,iy,iz;
		charge_map.GetClosestGridPoint(x_at,y_at,z_at, ix, iy,iz);
		
		double grid_x=x_min+ix*xstep;
		double grid_y=y_min+iy*ystep;
		double grid_z=z_min+iz*zstep;
		
        int i;
		double q[8],w[8];
		
		double dx1 = fabs(grid_x - x_at);
		if(dx1 < 1.0E-8) dx1 = 1.0E-8;
        double dx2 = fabs(grid_x + xstep - x_at);
        if(dx2 < 1.0E-8) dx2 = 1.0E-8;
		double dy1 = fabs(grid_y - y_at);
		if(dy1 < 1.0E-8) dy1 = 1.0E-8;
        double dy2 = fabs(grid_y + ystep - y_at);
		if(dy2 < 1.0E-8) dy2 = 1.0E-8;
		double dz1 = fabs(grid_z - z_at);
        if(dz1 < 1.0E-8) dz1 = 1.0E-8;
        double dz2 = fabs(grid_z + zstep - z_at);
		if(dz2 < 1.0E-8) dz2 = 1.0E-8;
		
		dx1 = 1.0/dx1;
		dx2 = 1.0/dx2;
		dy1 = 1.0/dy1;
		dy2 = 1.0/dy2;
		dz1 = 1.0/dz1;
		dz2 = 1.0/dz2;
		
		double ww = 0.0;
		
		if( fill_charge_mode == FILL_CHARGE_SUM_REV_DIST) 
        {
			w[0] =  dx1 + dy1 + dz1;
			w[1] =  dx2 + dy1 + dz1;
			w[2] =  dx1 + dy2 + dz1;
			w[3] =  dx2 + dy2 + dz1;
			w[4] =  dx1 + dy1 + dz2;
			w[5] =  dx2 + dy1 + dz2;
			w[6] =  dx1 + dy2 + dz2;
			w[7] =  dx2 + dy2 + dz2;			
		}
		else if( fill_charge_mode == FILL_CHARGE_PROD_REV_DIST) 
        {
			w[0] = dx1*dy1*dz1;
			w[1] = dx2*dy1*dz1;
			w[2] = dx1*dy2*dz1;
			w[3] = dx2*dy2*dz1;
			w[4] = dx1*dy1*dz2;
			w[5] = dx2*dy1*dz2;
			w[6] = dx1*dy2*dz2;
			w[7] = dx2*dy2*dz2;
		}
		
		
		for(i=0; i < 8; i++)
		{
			ww += w[i];
		}
		
		int ix2,iy2,iz2;
		
		//			if( fabs(at_charge) > 0.0000001)
		//			{
		//				PrintLog(" dx1 = %12.6f dx2 = %12.6f \n",dx1,dx2);
		//				PrintLog(" dy1 = %12.6f dy2 = %12.6f \n",dy1,dy2);
		//				PrintLog(" dz1 = %12.6f dz2 = %12.6f \n",dz1,dz2);
		//			}
		
		for(i=0; i < 8; i++)
		{
			q[i] = at_charge*w[i]/ww;
			
			if(i == 0) { ix2 = ix; iy2 = iy; iz2 = iz;}
			if(i == 1) { ix2 = ix+1; iy2 = iy; iz2 = iz;}
			if(i == 2) { ix2 = ix; iy2 = iy+1; iz2 = iz;}
			if(i == 3) { ix2 = ix+1; iy2 = iy+1; iz2 = iz;}
			if(i == 4) { ix2 = ix; iy2 = iy; iz2 = iz +1;}
			if(i == 5) { ix2 = ix+1; iy2 = iy; iz2 = iz +1;}
			if(i == 6) { ix2 = ix; iy2 = iy+1; iz2 = iz +1;}
			if(i == 7) { ix2 = ix+1; iy2 = iy+1; iz2 = iz +1;}
			
			//				if( fabs(at_charge) > 0.0000001)
			//				{
			//					PrintLog(" Ch placed to pt(%3d,%3d,%3d)=%12.6f \n", ix2,iy2,iz2, q[i]*el_vol);
			//				}
		}
		
		
		float* valPtr;
		valPtr=charge_map.GetValPtr(ix,iy,iz);
		*valPtr=*valPtr+q[0];
		valPtr=charge_map.GetValPtr(ix+1,iy,iz);
		*valPtr=*valPtr+q[1];
		valPtr=charge_map.GetValPtr(ix,iy+1,iz);
		*valPtr=*valPtr+q[2];
		valPtr=charge_map.GetValPtr(ix+1,iy+1,iz);
		*valPtr=*valPtr+q[3];
		valPtr=charge_map.GetValPtr(ix,iy,iz+1);
		*valPtr=*valPtr+q[4];
		valPtr=charge_map.GetValPtr(ix+1,iy,iz+1);
		*valPtr=*valPtr+q[5];
		valPtr=charge_map.GetValPtr(ix,iy+1,iz+1);
		*valPtr=*valPtr+q[6];
		valPtr=charge_map.GetValPtr(ix+1,iy+1,iz+1);
		*valPtr=*valPtr+q[7];
	}
		
	return true;
}
bool 
ElectrostMod::CalcIndCharge()
{   
    FillChargeMap();
	double val1x,val2x,val3x,valx,val;
	double val1y,val2y,val3y,valy;
	double val1z,val2z,val3z,valz;
	int i,j,k;
    int Nx,Ny,Nz;
	double xstep,ystep,zstep;
    xstep=el_pot_map.GetXstep();
    ystep=el_pot_map.GetYstep();
    zstep=el_pot_map.GetZstep();
	
    double el_vol = xstep*ystep*zstep; // volume element

	Nx=el_pot_map.GetNx();
	Ny=el_pot_map.GetNy();
	Nz=el_pot_map.GetNz();
	std::string fname;
	fname=phost_mset->GetName();
	fname=fname + "_ch.dat";
	FILE* fp=fopen(fname.c_str(),"w");

	ind_charge_map.SetGridCornersCoord( el_pot_map.GetXmin(), el_pot_map.GetYmin(),el_pot_map.GetZmin(),
		el_pot_map.GetXmax(), el_pot_map.GetYmax(),el_pot_map.GetZmax()); 
	
	ind_charge_map.SetDimensions(el_pot_map.GetNx(), el_pot_map.GetNy(),el_pot_map.GetNz());

	double tot_ind_charge = 0.0;
	double tot_pos_ind_charge = 0.0;
	double tot_neg_ind_charge = 0.0;

	for (k=1;k < Nz-1;k++)
	{
		for (j=1;j < Ny-1;j++)
		{
			for (i=1;i < Nx-1;i++)
			{
				val1x=el_pot_map.GetValue(i-1,j,k);
				val2x=el_pot_map.GetValue(i,j,k);
				val3x=el_pot_map.GetValue(i+1,j,k);
				valx=(val1x-2*val2x+val3x)/(xstep*xstep);
				
				val1y=el_pot_map.GetValue(i,j-1,k);
				val2y=el_pot_map.GetValue(i,j,k);
				val3y=el_pot_map.GetValue(i,j+1,k);
				
				valy=(val1y-2*val2y+val3y)/(ystep*ystep);
				
				val1z=el_pot_map.GetValue(i,j,k-1);
				val2z=el_pot_map.GetValue(i,j,k);
				val3z=el_pot_map.GetValue(i,j,k+1);
				
				valz=(val1z-2*val2z+val3z)/(zstep*zstep);
				
				val=valx+valy+valz;	
				
				val = val/HARTREE_TO_KT; // convert to kT

				double ind_charge;
				ind_charge= - (val/(4*PI)) - charge_map.GetValue(i,j,k);

//              ind_charge= (val/(4*PI));

				float grid_x,grid_y,grid_z;
				charge_map.GetXYZ(grid_x,grid_y,grid_z,i,j,k);

				
				if ( fabs(ind_charge) > 0.0001)
				{
					fprintf(fp,"%12.6f %12.6f %12.6f %12.6f \n",grid_x,grid_y,grid_z,ind_charge*el_vol);

                    tot_ind_charge += ind_charge*el_vol;

//					PrintLog(" Ind ch in (%3d,%3d,%3d)=%12.6f \n", i,j,k,ind_charge*el_vol);
					if(ind_charge > 0 )
					{
                        tot_pos_ind_charge += ind_charge*el_vol;
					}
					else
					{
                        tot_neg_ind_charge += ind_charge*el_vol;
					}

					float* valPtr;
					valPtr= ind_charge_map.GetValPtr(i,j,k);
					*valPtr=*valPtr+ind_charge;
				}	
			}
		}
	}	

	fclose(fp);
	PrintLog("Total induced charge = %9.3f \n",tot_ind_charge); 
	PrintLog("Total positive induced charge = %9.3f \n",tot_pos_ind_charge);
	PrintLog("Total negative induced charge = %9.3f \n",tot_neg_ind_charge); 

	return true;
}

bool
ElectrostMod::PlotIndCharge()
{
	PrintLog("Ind_Charge_Map \n");
    if( ind_charge_map.GetNx() <= 0 || ind_charge_map.GetNy() <= 0 ||
		ind_charge_map.GetNz() <= 0)
	{
		return false;
	}
	
	MolSet* pmset = GetMolSet();
	if( pmset == NULL)
		return false ;

	char buf[256];
	HaDisplayedSurface* sptr;

// Positive Potential Isosurface
	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return false;

	Surfaces.push_back(sptr);
	pmset->AddObject3D(sptr);
	
	bool result;
	result = sptr->calc_isosurf(&ind_charge_map, pot_isolevel);
	sptr->ColourUniform(255,0,0);
	sprintf(buf,"isopot_surf_%f",pot_isolevel);
	sptr->SetObjName(buf);

// Negative Potential Isosurface

	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return false;

	Surfaces.push_back(sptr);
	pmset->AddObject3D(sptr);

	result = sptr->calc_isosurf(&ind_charge_map, -pot_isolevel);
	sptr->ColourUniform(0,0,255);
	
	sprintf(buf,"isopot_surf_%f",-pot_isolevel);
	sptr->SetObjName(buf);


//	sptr->Print_info(cout,1);
	return result;    
}	

bool ElectrostMod::ColorMolSurfElPot()
{
	if( el_pot_map.GetNx() <= 0 || el_pot_map.GetNy() <= 0 ||
		el_pot_map.GetNz() <= 0)
	{
		return false;
	}
	
	MolSet* pmset = GetMolSet();
	if( pmset == NULL)
		return false;

	HaDisplayedSurface* sptr = pmset->GetMolSurface(FALSE);
	
	if( sptr == NULL)
	{
		sptr = pmset->CalcMolSurface();
	}
	
	if(sptr == NULL) return false;

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
		
	double fspan = elpot_high_val - elpot_low_val;
	
	if( fspan < 1E-10) 
	{
			ErrorInMod("ElectrostMod::ColorMolSurfElPot()",
				       "The difference between high and low boundary potential values is too small");
			return FALSE;
	}

	int i;
	for(i = 1; i <= numverts; i++ )
	{
		double x = sptr->verts(1,i);
		double y = sptr->verts(2,i);
		double z = sptr->verts(3,i);

		double fval = el_pot_map.GetInterpolValAtPoint(x,y,z);
	
		double didx_col = 1.0 + ((double) (ncol-1) ) * ( (fval - this->elpot_low_val)/fspan);

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
//>mikola 27 july 2006
bool ElectrostMod::ColorDotStruct(DotStruct* dotstruct)
{
  if( el_pot_map.GetNx() <= 0 || el_pot_map.GetNy() <= 0 ||
      el_pot_map.GetNz() <= 0)
  {
    return false;
  }
	
  MolSet* pmset = GetMolSet();
  if( pmset == NULL)
    return false;

  //HaDisplayedSurface* sptr = pmset->CalcMolSurface();
	
  if(dotstruct == NULL)
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

  int numverts = dotstruct->dots.size();
  /*if(sptr->colors.size() < numverts)
  {
    ErrorInMod("ElectrostMod::ColorMolSurfElPot()","Colors array size is smaller that number of verticies");
    sptr->ColourUniform(200,200,200);
    return FALSE;
  }*/
		
  double fspan = elpot_high_val - elpot_low_val;
	
  if( fspan < 1E-10) 
  {
    ErrorInMod("ElectrostMod::ColorMolSurfElPot()",
               "The difference between high and low boundary potential values is too small");
    return FALSE;
  }

  int i;
  for(i = 1; i <= numverts; i++ )
  {
    double x = dotstruct->dots[i].GetX();
    double y = dotstruct->dots[i].GetY();
    double z = dotstruct->dots[i].GetZ();

    double fval = el_pot_map.GetInterpolValAtPoint(x,y,z);
	
    double didx_col = 1.0 + ((double) (ncol-1) ) * ( (fval - this->elpot_low_val)/fspan);

    int idx_col = (int) didx_col;
    if( (didx_col - (double)idx_col ) > 0.5 )
      idx_col++;

    if( idx_col < 1) idx_col = 1;
    if( idx_col > ncol ) idx_col = ncol;

    col = col_num(idx_col);
    dotstruct->dots[i].col=col_num(ncol);
    //sptr->colors(i) = col;
  }

  return true;
}
//<mikola 27 july 2006
bool ElectrostMod::BuildPotVdwDots()
{
	if( el_pot_map.GetNx() <= 0 || el_pot_map.GetNy() <= 0 ||
		el_pot_map.GetNz() <= 0)
	{
		return false;
	}
	
	MolSet* pmset = GetMolSet();
	if( pmset == NULL)
		return false ;

	HaMolView* pView = pmset->GetActiveMolView();

	if(pView)pView->CalculateDotSurface(dots_number);
//>mikola 27 july 2006
 Object3D* obj=pmset->ViewObjects.back();
 if(obj->GetObjType()!=OBJ3D_DOT_SURFACE)return false;
 DotStruct *dotstructpmset=(DotStruct*)obj;
 if(!ColorDotStruct(dotstructpmset))return false;
 //<mikola 27 july 2006
	return true;
}


double ElectrostMod::CalcETReorgEne()
{
	double reorg_ene = -1.0;
	if(phost_mset == NULL)
		return reorg_ene;

	AtomGroup* pdon    = phost_mset->GetAtomGroupByID("DONOR");
	AtomGroup* pacc    = phost_mset->GetAtomGroupByID("ACCEPTOR");

	if( !pdon )
	{
		PrintLog(" Error in ElectrostMod::CalcETReorgEne() \n");
		PrintLog(" Donor is not set for the molecular set  \n");
		PrintLog(" ET reorganization energy won't be computed \n");
		return reorg_ene;
	}

	int no_acc_flag = FALSE;
    
	if(!pacc || pacc->size() == 0)
	{
		PrintLog(" Warning: \n");
		PrintLog(" Acceptor is not set or contains no atoms \n");
        PrintLog(" ET reorganization energy caclulations will be done only for the donor \n");
		no_acc_flag = TRUE;
	}

	HaAtom* aptr;
	AtomDoubleMap save_charges;
	
	AtomIteratorMolSet aitr(phost_mset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{	
		save_charges[aptr] = aptr->GetCharge();
		aptr->SetCharge(0.0);
	}
	
	AtomIteratorAtomGroup aitr_d(pdon);
	AtomIteratorAtomGroup aitr_a(pacc);

	double ch;
	ch = 1.0/pdon->size();
	for(aptr= aitr_d.GetFirstAtom(); aptr; aptr= aitr_d.GetNextAtom())
	{
		aptr->SetCharge(ch);
	}

    if( !no_acc_flag )
	{
		ch= -1.0/pacc->size();
		for(aptr= aitr_a.GetFirstAtom(); aptr; aptr= aitr_a.GetNextAtom())
		{
			aptr->SetCharge(ch);
		}
	}

	double xmin, xmax, ymin, ymax, zmin, zmax;
	phost_mset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

	SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);
	
	double ene1, ene2, old_ch;
	bool bres;

	bres = run(RUN_FOREGROUND);
	ene1 = this->tot_ene;

	if(!bres)
	{
	// restore charges
	    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{	
			old_ch = save_charges[aptr];
			aptr->SetCharge(old_ch);
		}
		return reorg_ene;
	}

	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom() )
	{
		aptr->UnSelect();
	}

	for(aptr= aitr_d.GetFirstAtom(); aptr; aptr= aitr_d.GetNextAtom())
	{
		aptr->Select();
	}
	
	if( !no_acc_flag)
	{
		for(aptr= aitr_a.GetFirstAtom(); aptr; aptr= aitr_a.GetNextAtom())
		{
			aptr->Select();
		}
	}

	phost_mset->p_save_opt_default->save_selected = TRUE;

	double epsout_save;
	epsout_save = epsout;
	epsout = 2.0;
		
	bres = run(RUN_FOREGROUND);

	ene2 = this->tot_ene;

	phost_mset->p_save_opt_default->save_selected = FALSE;
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom() )
	{
		aptr->Select();
	}
	epsout = epsout_save;

	if(!bres)
	{
	// restore charges
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{	
			old_ch = save_charges[aptr];
			aptr->SetCharge(old_ch);
		}
		return reorg_ene;
	}

	reorg_ene = (ene2 - ene1)/EV_TO_KT;
	PrintLog(" ene1= %12.3f ene2= %12.3f kT \n ",ene1,ene2);
	PrintLog(" ET Reorganization energy is %12.3f eV \n",reorg_ene);
	// restore charges
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{	
		old_ch = save_charges[aptr];
		aptr->SetCharge(old_ch);
	}

	return reorg_ene;
}

bool ElectrostMod::CalcAltStatePK(AltChemState* alt_res_st, AtomGroup* active_atoms )
{
	if(phost_mset == NULL)
		return false;
	
	if(alt_res_st == NULL) return false;	
	if(!alt_res_st->active_flag) return true;

	HaResidue* pres = (HaResidue*) alt_res_st->GetHostAtomGroup();

	if(pres == NULL)
	{
		ErrorInMod("ElectrostMod::CalcResPK()","Residue pointer is NULL");
		return false;
	}

	HaAtom* aptr; 
	map<void*,double,less<void*> > map_aptr_ch;
	AtomIteratorMolSet  aitr_mset(phost_mset);
	AtomIteratorResidue aitr_res(pres);
	
	for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
	{
		map_aptr_ch[(void*)aptr] = aptr->GetCharge();
	}

	double xmin, xmax, ymin, ymax, zmin, zmax;
	phost_mset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

	int iaxx = 0;
	
	for( aptr = aitr_mset.GetFirstAtom(); aptr; aptr = aitr_mset.GetNextAtom() )
	{
		if(active_atoms != NULL)
			aptr->UnSelect();
		else
			aptr->Select();
	}
	
	if( active_atoms != NULL )
	{
		AtomIteratorAtomGroup as_itr(active_atoms);
		for(aptr = as_itr.GetFirstAtom(); aptr; aptr = as_itr.GetNextAtom())
			aptr->Select();
	}
	
	std::string unprot_res_name;
	std::string prot_res_name;
	
	double ene1, ene2, ene3, ene4;
	bool bres;
	
	phost_mset->p_save_opt_default->save_selected = TRUE; // Only selected atoms are saved
	SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);
	
	alt_res_st->SetAltCharges(0.0);
	
	bres = run(RUN_FOREGROUND);
	ene1 = this->tot_ene;
	if(!bres) { return false; }
	
	alt_res_st->SetAltCharges(1.0);
	
	bres = run(RUN_FOREGROUND);
	ene2 = this->tot_ene;
	if(!bres) { return false; }
		
	for(aptr= aitr_mset.GetFirstAtom(); aptr; aptr= aitr_mset.GetNextAtom() )
	{
		aptr->UnSelect();
	}
	for(aptr= aitr_res.GetFirstAtom(); aptr; aptr= aitr_res.GetNextAtom() )
	{
		aptr->Select();
	}
	
	alt_res_st->SetAltCharges(0.0);
				
	bres = run(RUN_FOREGROUND);
	
	ene3 = this->tot_ene;
	if(!bres) { return false; }
				
	alt_res_st->SetAltCharges(1.0);
				
	bres = run(RUN_FOREGROUND);
	ene4 = this->tot_ene;
	if(!bres) { return false; }
												
	phost_mset->p_save_opt_default->save_selected = FALSE; // restore default mode of atoms saving
								
//
//  ene1 - full system (protein)      unmodified chemical state  (in kT)
//  ene2 - full system (protein)      modified   chemical state  (in kT)
//  ene3 - reference system (solvent) unmodified chemical state  (in kT)
//  ene4 - reference system (solvent) modified   chemical state  (in kT)
//
//  ddG - relative change of the energy of modified state vs unmodified compared to reference system
//
//  if ddG  < 0 for AltChemState::PROTONATED    - protonated state is stabilized delta_pK -> negative
//  if ddG  < 0 for AltChemState::UNPROTONATED  - unprotonated state is stabilized delta_pK -> positive
//  if ddG  < 0 for AltChemState::OXIDIZED      - oxidized state is stabilized delta_rp -> negative
//  if ddG  < 0 for AltChemState::REDUCED       - reduced  state is stabilized delta_rp -> negative
//
//	double ddG = (ene2 - ene1) - (ene4 - ene3); jose
//	double delta_pK, delta_rp; jose
	double delta_rp; //jose added
	ddG = (ene2 - ene1) - (ene4 - ene3); //added jose

	if(alt_res_st->alt_state_type == alt_res_st->alt_state_type.PROTONATED)
	{
	   delta_pK = -ddG/log(10.0);
	}
	else if(alt_res_st->alt_state_type == alt_res_st->alt_state_type.UNPROTONATED)
	{
	   delta_pK = ddG/log(10.0);
	}

	if(alt_res_st->alt_state_type == alt_res_st->alt_state_type.OXIDIZED)
	{
	   delta_rp = ddG/EV_TO_KT;
	}
	else if(alt_res_st->alt_state_type == alt_res_st->alt_state_type.REDUCED)
	{
       delta_rp = -ddG/EV_TO_KT;
	}
	
	if( alt_res_st->alt_state_type == alt_res_st->alt_state_type.UNPROTONATED || 
		alt_res_st->alt_state_type == alt_res_st->alt_state_type.PROTONATED)
	{
	   alt_res_st->pk = alt_res_st->std_pk + delta_pK;
				
	   PrintLog("Computing pK shift for %s  Alt Chem State %s \n", (pres->GetRef()).c_str(), alt_res_st->id.c_str());
	   PrintLog("ene1(in protein/unmodified)= %12.6f \n ene2(in protein/modified)= %12.6f \n ene3(in solvent/unmodified)= %12.6f \n ene4(in solvent/modified)= %12.6f \n",
	 				  ene1,ene2,ene3,ene4);
	   PrintLog("DeltaDelta G protonation/deprotonation: solvent -> protein ddG= %12.6f (kT) \n",ddG);
	   PrintLog("Calculated pKa shift = %12.6f  and predicted pKa = %12.6f \n", delta_pK, alt_res_st->pk);
	}
	else if( alt_res_st->alt_state_type == alt_res_st->alt_state_type.OXIDIZED || 
             alt_res_st->alt_state_type == alt_res_st->alt_state_type.REDUCED)
	{
       alt_res_st->pk = alt_res_st->std_pk + delta_rp;
	   PrintLog("ene1(in protein/unmodified)= %12.6f \n ene2(in protein/modified)= %12.6f \n ene3(in solvent/unmodified)= %12.6f \n ene4(in solvent/modified)= %12.6f \n",
	 				  ene1,ene2,ene3,ene4);
	   PrintLog("DeltaDelta G oxidation/reduction: solvent -> protein ddG= %12.6f (kT) \n",ddG);
	   PrintLog("Calculated redox pot shift = %12.6f  and predicted redox pot = %12.6f \n", delta_rp, alt_res_st->pk);
	}
							
// Restore charge state of the residue

	for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
	{
		double ch_at = map_aptr_ch[(void*)aptr];
		aptr->SetCharge(ch_at);
	}
	
	return true;
}


bool
ElectrostMod::CalcRedoxPotShft()
{
	if(phost_mset == NULL)
		return false;

    AtomGroup* pdon = phost_mset->GetAtomGroupByID("DONOR");
	if( !pdon )
	{
		PrintLog(" Error in ElectrostMod::CalcRedoxPotShft() \n");
		PrintLog(" Donor is not set for the molecular set \n");
		return false;
	}
	
	HaAtom* aptr;
	double ene1, ene2, ene3, ene4;
	bool bres;
    double epsout_save = epsout;
	double rionst_save = rionst;
    char buf[256]; int j;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	phost_mset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

	SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);

	int m;

	for( m = 0 ; m < 2; m++)
	{
		phost_mset->p_save_opt_default->save_selected = FALSE;
		
		bres = run(RUN_FOREGROUND);
		
		if( m == 0)
			ene1 = this->tot_ene;
		else
			ene3 = this->tot_ene;

		if(!bres)
			return false;
		
        AtomIteratorMolSet aitr(phost_mset);	
    	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom() )
		{
			aptr->UnSelect();
		}
		
		AtomIteratorAtomGroup aitr_d(pdon);
		for(aptr= aitr_d.GetFirstAtom(); aptr; aptr= aitr_d.GetNextAtom())
		{
			aptr->Select();
		}
		
		if(rdx_shft_mode & RDX_SHFT_VS_VAC)
		{
			epsout = 1.0;
			rionst = 0.0;
		}
		
		phost_mset->p_save_opt_default->save_selected = TRUE;

		bres = run(RUN_FOREGROUND);

		epsout = epsout_save;
		rionst = rionst_save;


		phost_mset->p_save_opt_default->save_selected = FALSE;
			
    	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom() )
		{
			aptr->Select();
		}
			
		if(m == 0)
			ene2 = this->tot_ene;
		else
			ene4 = this->tot_ene;

		if(!bres)
			return false;
		
		double delt_ch; 
		delt_ch= 1.0/pdon->size();
		
		for(aptr= aitr_d.GetFirstAtom(); aptr; aptr= aitr_d.GetNextAtom())
		{
			if(m == 0)
			{
				aptr->SetCharge(aptr->GetCharge() + delt_ch); // generated oxidized state charge distribution of the donor
			}
			else
			{
				aptr->SetCharge(aptr->GetCharge() - delt_ch); // restore reduced state charge distribution of the donor 
			}
		}
	}

	PrintLog(" ene1= %10.4f  ene2= %10.4f ene3= %10.4f ene4= %10.4f kT\n", ene1,ene2,ene3,ene4);
	PrintLog(" Electrostatic contributions to the redox-potential of the donor is : \n"); 

	rdx_shift = ((ene3 - ene1) - (ene4 -ene2))/EV_TO_KT ;
	

	if(rdx_shft_mode & RDX_SHFT_VS_VAC) 
		j= sprintf(buf,"%s", " Delta E0 (vs vaccuum) =  "); 
	else if( rdx_shft_mode & RDX_SHFT_VS_SOLV)
		j= sprintf(buf,"%s", " Delta E0 (vs solvent) =  ");

	sprintf(buf+j,"%12.6f V", rdx_shift); 
	PrintLog(buf);

	return true;
}

double
ElectrostMod::CalcAvgPotOn(PointContainer* ptlist)
{
	if(el_pot_map.GetNx() < 1 || el_pot_map.GetNx() < 1 || el_pot_map.GetNz() < 1)
	{
		PrintLog("Error in ElectrostMod::CalcAvgPotOn() \n");
        PrintLog("Electrostatic Poential Map is not set \n");
        return 0.0;
	}

	int n = ptlist->GetNumPt();
	if( n == 0)
	{
		PrintLog("Error in ElectrostMod::CalcAvgPotOn() \n");
        PrintLog("Points Collection is empty \n");
        return 0.0;		
	}

	PointIterator* pt_itr = ptlist->GetPointIteratorPtr();
	if(pt_itr == NULL) return 0.0;
	Vec3D* pt;
	
	double av = 0.0;
	for(pt = pt_itr->GetFirstPt(); pt; pt = pt_itr->GetNextPt())
	{
		double x = pt->GetX();
		double y = pt->GetY();
		double z = pt->GetZ();
		double val = el_pot_map.GetInterpolValAtPoint(x,y,z);
		av += val;
	}
	delete pt_itr;

	av = av/n;
	return av;
}


bool ElectrostMod::ReadTotEne(double &tot_ene, const char* fname)
{
//	PrintLog(" In ReadTotEne() pt 1 \n");
#if defined(INT_DELPHI)
//	PrintLog(" In ReadTotEne() pt 2\n");
	return true;
#endif
//	PrintLog(" In ReadTotEne() pt 3\n");

	std::string title;
	tot_ene = 0.0;
	std::ifstream is(fname);
	if(is.fail())
	{
		cerr << " Error in ElectrostMod::ReadTotEne() " << endl;
		cerr << " Cannot open file " << fname << endl;
		return false;
	}
	is >> tot_ene;
	if(is.fail())
	{
		cerr << " Error in ElectrostMod::ReadTotEne() " << endl;
		cerr << " Fail to read energy from file " << fname << endl;
		return false;		
	}
	return true;
}


bool
ElectrostMod::SetBoundaryAtoms(const double xmin, const double ymin, const double zmin, 
						   const double xmax, const double ymax, const double zmax)
{
	min_coord.SetX(xmin - 4.0);
    min_coord.SetY(ymin - 4.0);
	min_coord.SetZ(zmin - 4.0);
	max_coord.SetX(xmax + 4.0);
    max_coord.SetY(ymax + 4.0);
	max_coord.SetZ(zmax + 4.0);

	return true;
}

bool
ElectrostMod::ClearBoundaryAtoms()
{
	min_coord.SetX(0.0);
    min_coord.SetY(0.0);
	min_coord.SetZ(0.0);
	max_coord.SetX(0.0);
    max_coord.SetY(0.0);
	max_coord.SetZ(0.0);

	return true;
}


bool
ElectrostMod::AddBoundaryAtoms()
{	
// If boundary atoms are not set do nothing
    if( Vec3D::CalcDistance(&min_coord,&max_coord) < 1.0E-8) 
		return false;

	FILE*  fcoord = fopen ( coord_file_name.c_str(), "a");
	if(fcoord == NULL)
	{
		cerr << " Error in HaDelphi::AddBoundaryAtoms() " << endl;
		cerr << " Failed to open to append file " << coord_file_name.c_str() << endl; 
		return false;
	}
	
	int i_at=40000;
	int ires= 8000;

	double xmin1 = min_coord.GetX_Ang();
	double ymin1 = min_coord.GetY_Ang();
	double zmin1 = min_coord.GetZ_Ang();

	double xmax1 = max_coord.GetX_Ang();
	double ymax1 = max_coord.GetY_Ang();
	double zmax1 = max_coord.GetZ_Ang();
	
	fprintf(fcoord, "ATOM  %5d %.4s %.3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00 \n",i_at++," xx ","AXX",ires, 
		     xmin1, ymin1, zmin1);
	
	fprintf(fcoord, "ATOM  %5d %.4s %.3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00 \n",i_at++," xx ","AXX",ires, 
		     xmax1, ymin1, zmin1);

	fprintf(fcoord, "ATOM  %5d %.4s %.3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00 \n",i_at++," xx ","AXX",ires, 
		     xmin1, ymax1, zmin1);
	
	fprintf(fcoord, "ATOM  %5d %.4s %.3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00 \n",i_at++," xx ","AXX",ires, 
		     xmin1, ymin1, zmax1);


	fclose(fcoord);
	return true;
}

bool 
ElectrostMod::run(RunMode rmode)
{
#if !defined(INT_DELPHI)
	SaveParamFile();
	SaveChargeFile();
	SaveRadiusFile();
	SaveCoordFile();
#endif

	RunDelphi(rmode);
	bool bres = ReadTotEne(tot_ene);
	return bres;
}

bool ElectrostMod::RunDelphi(RunMode rmode)
{
	std::string cmd;
	std::string copy_cmd;

#if !defined(INT_DELPHI)
	ha_copy_file(param_file_name.c_str(), "fort.10");
	ha_copy_file(charge_file_name.c_str(),"fort.12");
	ha_copy_file(radius_file_name.c_str(),"fort.11");
	ha_copy_file(coord_file_name.c_str(), "fort.13");
#endif
	
#ifdef USE_CREATE_PROCESS
   PROCESS_INFORMATION pi;    
   STARTUPINFO si; 
   SECURITY_ATTRIBUTES sa;

   DWORD creation_flags=0;

   memset( &si, 0,sizeof(si) );
  memset( &pi, 0,sizeof(pi) );

	HANDLE log_file_handle;

   log_file_handle = CreateFile(log_file_name.c_str(),
	   GENERIC_WRITE | GENERIC_READ, 
	   FILE_SHARE_READ, &sa,
	   CREATE_ALWAYS, 
	   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING,
	   NULL);

//  int fi_log= _open(log_file_name.c_str(), _O_CREAT | _O_TRUNC,
//	                _S_IREAD | _S_IWRITE );
//  log_file_handle=  (HANDLE) _get_osfhandle(fi_log);

   si.cb = sizeof(si);    // Start the child process.
   si.dwFlags = STARTF_USESTDHANDLES;
   si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
//   si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
   si.hStdOutput = log_file_handle;
   si.hStdError  = GetStdHandle(STD_ERROR_HANDLE);
	

//   si.wShowWindow=SW_SHOW; // probably not nessessary
   

   char buf[256];
   int j=0;
   j=sprintf(buf+j,"%s","qdiffxs");
   j=sprintf(buf+j,"%s%s"," > ",log_file_name.c_str());
	
//  creation_flags = CREATE_NEW_CONSOLE;
//	creation_flags = DETACHED_PROCESS;
	
	if(!CreateProcess( NULL, buf ,
		NULL, NULL, TRUE, creation_flags, 
		NULL, NULL, &si, &pi) )
	{
		DWORD error_code = GetLastError();
		::AfxMessageBox("Error starting DELPHI program " );
	}
 
#endif
	
	StrVec delphi_args;

#if defined(_MSC_VER)
	delphi_args.push_back(" > ");
#endif
    StrVec prog_output;

#if defined(INT_DELPHI)
	MolSet* pmset = GetMolSet();
	int_4 mgrid = nx;
	int_4 na = pmset->GetNAtoms();

	size_t tot_size = sizeof(freal)*(mgrid*mgrid*mgrid + 1)*2 + 4*7*(mgrid*mgrid*mgrid) + 
		              4*sizeof(freal)*na + 7*sizeof(freal)*(mgrid*mgrid*2 + 1) +
					  4*(mgrid*mgrid*20 + 1)*15;

	tot_size += sizeof(freal)*(mgrid*mgrid*mgrid + 1); // sf1,sf2
	if(nnit > 0 ) 
	{
		tot_size += sizeof(freal)*(mgrid*mgrid*mgrid + 1)*2; // qmap1,qmap2,debmap1,debmap2
	}

	void* xptr = malloc(tot_size);
	bool alloc_phi = el_pot_map.SetDimensions(mgrid,mgrid,mgrid);

    if(xptr == NULL || !alloc_phi )
	{
		if(!xptr) free(xptr);
		el_pot_map.clear();
		PrintLog(" Error Allocating Memory in ElectrostMod::RunDelphi()");
		return false;
	}

	freal* phi   =  el_pot_map.GetFieldPtr();
    freal* phi1  = (freal*) xptr;
	freal* phi2  = phi1 + (mgrid*mgrid*mgrid +1)/2;
	freal* phi3  = phi2 + (mgrid*mgrid*mgrid +1)/2;
    freal* sf1   = phi3 + (mgrid*mgrid*mgrid +1);
	freal* sf2   = sf1  + (mgrid*mgrid*mgrid +1)/2;
	freal* qmap1 = sf2  + (mgrid*mgrid*mgrid +1)/2;
	freal* qmap2; freal* debmap1; freal* debmap2; freal* achrg;
	if(nnit >  0 )
	{
		qmap2   = qmap1   + (mgrid*mgrid*mgrid +1)/2;
		debmap1 = qmap2   + (mgrid*mgrid*mgrid +1)/2;
		debmap2 = debmap1 + (mgrid*mgrid*mgrid +1)/2;
		achrg   = debmap2 + (mgrid*mgrid*mgrid +1)/2; 	
	}
	else
	{
		qmap2   = qmap1;
		debmap1 = qmap2;
		debmap2 = debmap1;
		achrg   = debmap2;
	}
	
	int_4* ieps  = (int_4*)(achrg + 4*na);
	int_4* ieps2 = ieps  + 3*mgrid*mgrid*mgrid;
	int_4* ideb  = ieps2 + 3*mgrid*mgrid*mgrid;

	freal* bndx1 = (freal*)(ideb + mgrid*mgrid*mgrid);
	freal* bndx2 = bndx1 + mgrid*mgrid*2 + 1;
	freal* bndx3 = bndx2 + mgrid*mgrid*2 + 1;
	freal* bndx4 = bndx3 + mgrid*mgrid*2 + 1;
	freal* bndx  = bndx4 + mgrid*mgrid*2 + 1;
	freal* bndy  = bndx  + mgrid*mgrid*2 + 1;
	freal* bndz  = bndy  + mgrid*mgrid*2 + 1;
	freal* db    = bndz  + mgrid*mgrid*2 + 1;
	int_4* ibgrd = (int_4*)( db + (mgrid*mgrid*20 + 1)*6 ); // was *10
	int_4* idpos = ibgrd + (mgrid*mgrid*20 + 1)*4 ;
	int_4* ioff  = idpos + (mgrid*mgrid*20 + 1);
	int_4* iepsv = ioff + (mgrid*mgrid*20 + 1)*3;

//	if(phi3 == NULL || phi2 == NULL || phi1 == NULL || phi == NULL)
//	{
//		if(phi != NULL) free(phi);
//		if(phi != NULL) free(phi1);
//		if(phi != NULL) free(phi2);
//		if(phi != NULL) free(phi3);
//	}


	delphi_((float*)xptr,&mgrid,phi,phi1,phi2,phi3,achrg,
		    ieps, ieps2, ideb, 
			sf1, sf2, qmap1, qmap2, debmap1, debmap2,
			bndx1, bndx2, bndx3, bndx4, bndx, bndy, bndz,
			db, ibgrd, idpos, ioff, iepsv);

 //
 /*int i;
 int gridsize[3]={mgrid,mgrid,mgrid};
 VectorField3D VF(gridsize,1.0,3);
 float** V=VF.V;
 for(i=0;i<mgrid*mgrid*mgrid/2;i++)
 {
   V[0][i]=qmap1[i];
   V[0][i+(mgrid*mgrid*mgrid +1)/2]=qmap2[i];
   V[1][i]=debmap1[i];
   V[1][i+(mgrid*mgrid*mgrid +1)/2]=debmap2[i];
 }
 for(i=0;i<mgrid*mgrid*mgrid;i++)
 {
   V[2][i]=ieps[i];
 }
 VF.WriteToFile("delphi.gz");*/
 //
 
	free(xptr);
	if(phiwrt) el_pot_map.SaveToFile(elfield_fname.c_str());
#else
	HarlemApp::RunExternalProgram(rmode, "qdiffxs",delphi_args, prog_output,FALSE);
#endif
	return true;
	
}

#if defined(INT_DELPHI)


void rdprm2_(logical* iautocon,
		         freal* epsin, freal* epsout,
				 freal* rionst, freal* exrad, freal* radprb,
				 int_4* ibctyp, logical* iper, int_4* nlit, int_4* nnit,
				 logical* iconc, logical* imem, 
				 int_4* icon1, int_4* icon2)
{
	MolSet* pmset = GetCurMolSet();
	ElectrostMod* elmod = pmset->GetElectrostMod();
	if(elmod == NULL) return;
    
	*iautocon = FALSE;
	if(elmod->nlit < 0 ) 
	{
		*nlit = 0;
		*iautocon = TRUE;
	}
	else
	{
		*nlit = elmod->nlit;
	}
    
	scaleq_.igrid = elmod->nx;

	*epsin      = elmod->epsi;
	*epsout     = elmod->epsout;
	*rionst     = elmod->rionst;
	*exrad      = elmod->exrad;
	*radprb     = elmod->radprb;
	*ibctyp     = elmod->boundary+1;
	*iper       = elmod->iper[0];
	*(iper+1)   = elmod->iper[1];
	*(iper+2)   = elmod->iper[2];
	*nnit       = elmod->nnit;
	*iconc      = elmod->iconc;
	*imem       = elmod->imem;
	*icon1      = elmod->icon1;
	*icon2      = elmod->icon2;
}

void setatq_(freal* xn2, freal* rad3, freal* chrgv4, int_4* natom, freal* gr_cent, freal* pscale)
{

	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	ElectrostMod* elmod = pmset->GetElectrostMod();
	if(elmod == NULL) return;

	int na = 0;
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);

	bool calc_range = true;

	if( (elmod->max_coord.GetX() - elmod->min_coord.GetX()) > 0.001 &&
        (elmod->max_coord.GetY() - elmod->min_coord.GetY()) > 0.001 &&
		(elmod->max_coord.GetZ() - elmod->min_coord.GetZ()) > 0.001 ) 
	{
		calc_range = false;
	}

	double xmin1, ymin1,zmin1,xmax1,ymax1,zmax1;

	if(calc_range)
	{
		xmin1 = 10000.0;
		ymin1 = 10000.0;
		zmin1 = 10000.0;
		xmax1 = -10000.0;
		ymax1 = -10000.0;
		zmax1 = -10000.0;
	}
	else
	{
		xmin1 = elmod->min_coord.GetX_Ang();
		ymin1 = elmod->min_coord.GetY_Ang();
		zmin1 = elmod->min_coord.GetZ_Ang();

		xmax1 = elmod->max_coord.GetX_Ang();
		ymax1 = elmod->max_coord.GetY_Ang();
		zmax1 = elmod->max_coord.GetZ_Ang();
	}

	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if( pmset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
		xn2[3*na]   = aptr->GetX_Ang();
		xn2[3*na+1] = aptr->GetY_Ang();
		xn2[3*na+2] = aptr->GetZ_Ang();
		rad3[na]    = aptr->radius;
        chrgv4[na]  = aptr->GetCharge();
		
		if( calc_range)
		{
			double tmp;
			tmp = xn2[3*na]   - rad3[na]; if(tmp < xmin1) xmin1 = tmp;
			tmp = xn2[3*na+1] - rad3[na]; if(tmp < ymin1) ymin1 = tmp;
			tmp = xn2[3*na+2] - rad3[na]; if(tmp < zmin1) zmin1 = tmp;
			tmp = xn2[3*na]   + rad3[na]; if(tmp > xmax1) xmax1 = tmp;
			tmp = xn2[3*na+1] + rad3[na]; if(tmp > ymax1) ymax1 = tmp;
			tmp = xn2[3*na+2] + rad3[na]; if(tmp > zmax1) zmax1 = tmp;
		}
		na++;
	}


// Add boundary atoms


	gr_cent[0] = (xmin1 + xmax1)/2.0;
	gr_cent[1] = (ymin1 + ymax1)/2.0;
	gr_cent[2] = (zmin1 + zmax1)/2.0;

	double range = (xmax1 - xmin1);
	if( (ymax1 - ymin1) > range ) range = (ymax1 - ymin1);
    if( (zmax1 - zmin1) > range ) range = (zmax1 - zmin1);

	int ngrid = elmod->nx;

	range = range*100.0/elmod->perfil;

	float scale;

	if( elmod->perfil > 0.000001)  
	{
		scale = (ngrid - 1)/range;
	}
	else
	{
		scale = 0.5;
	}

	*pscale = scale;

	double gxmin = (gr_cent[0] - range*0.5);
    double gymin = (gr_cent[1] - range*0.5);
	double gzmin = (gr_cent[2] - range*0.5);
    double gxmax = (gr_cent[0] + range*0.5);
    double gymax = (gr_cent[1] + range*0.5);
	double gzmax = (gr_cent[2] + range*0.5);
	  
	elmod->el_pot_map.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax); 


/*    if( fabs(xmax1 - xmin1) > 0.01 || fabs(ymax1 - ymin1) > 0.01 || fabs(zmax1 - zmin1) > 0.01)
	{
		xn2[3*na]   = xmin1;
c		xn2[3*na+1] = ymin1;
c		xn2[3*na+2] = zmin1;
c		rad3[na]    = 0.0;
c       chrgv4[na]  = 0.0;
c		na++;
c		xn2[3*na]   = xmax1;
c		xn2[3*na+1] = ymin1;
c		xn2[3*na+2] = zmin1;
c		rad3[na]    = 0.0;
c       chrgv4[na]  = 0.0;
c		na++;
c		xn2[3*na]   = xmin1;
c		xn2[3*na+1] = ymax1;
c		xn2[3*na+2] = zmin1;
c		rad3[na]    = 0.0;
c       chrgv4[na]  = 0.0;
c		na++;
c		xn2[3*na]   = xmin1;
c		xn2[3*na+1] = ymin1;
c		xn2[3*na+2] = zmax1;
c		rad3[na]    = 0.0;
c       chrgv4[na]  = 0.0;
c		na++;
c	}
*/
	*natom = na;

}

void getpotm_(int_4* p_ngrid, freal* potmap, 
			 freal* pxmin,  freal* pymin, freal* pzmin, 
             freal* pxmax,  freal* pymax, freal* pzmax)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	ElectrostMod* elmod = pmset->GetElectrostMod();
	if(elmod == NULL) return;

	int ngrid = *p_ngrid;
 
	elmod->el_pot_map.SetGridCornersCoord( *pxmin, *pymin, *pzmin, *pxmax, *pymax, *pzmax); 

	elmod->el_pot_map.SetDimensions(ngrid, ngrid, ngrid);

	float* fpout = elmod->el_pot_map.GetFieldPtr();
	freal* fpin  = potmap;
	
	int np = ngrid*ngrid*ngrid;
	int i;
	for(i= 0; i < np; i++)
	{
		*fpout = *fpin;
		fpin++;
		fpout++;
	}
}

void geteneq_(freal* ptot_ene)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	ElectrostMod* elmod = pmset->GetElectrostMod();
	if(elmod == NULL) return;
	PrintLog(" Setting Total energy %12.6f \n",*ptot_ene);

	elmod->tot_ene = *ptot_ene;
}

#endif











