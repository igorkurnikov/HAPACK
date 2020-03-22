/*! \file hasurface.cpp

    Classes to define 3D-fields and surfaces 

    \author Igor Kurnikov 
	\date 1999-2002
*/

#define HASURFACE_CPP

#include <float.h>
#include <math.h>
#include "hasurface.h"
#include "hamolview.h"
//>mikola 29,July, 2006
#include "mapio.h"
#include "haobject.h"
#include <tinyxml.h>
#include "haxml.h"
 
//<mikola 29,July, 2006
/*HaField3D::HaField3D()
{
	m_Nx = m_Ny = m_Nz = 0;
}*/

HaField3D::HaField3D(float *fvec,int new_Nx, int new_Ny, int new_Nz, bool deligate_control)
{
  m_Nx = m_Ny = m_Nz = 0;
  if(new_Nx>0 && new_Ny> 0 && new_Nz > 0)
  {
    m_Nx = new_Nx;
    m_Ny = new_Ny;
    m_Nz = new_Nz;
  }
  if (fvec != NULL)
  {
	m_field_data.set_ext_alloc(fvec, m_Nx*m_Ny*m_Nz, deligate_control);
  }
}

HaField3D::~HaField3D()
{

}

void 
HaField3D::clear()
{
	SetDimensions(0,0,0);

	xmin = xmax = ymin = zmax = ymin = ymax = 0.0;
	xstep = ystep = zstep = 0.0;
}

void
HaField3D::copy_from(const HaField3D& ref_field)
{
	clear();
	SetDimensions(ref_field.GetNx(),ref_field.GetNy(),ref_field.GetNz());
	SetGridCornersCoord(ref_field.GetXmin(),ref_field.GetYmin(),ref_field.GetZmin(),
		                ref_field.GetXmax(),ref_field.GetYmax(),ref_field.GetZmax());
	
	m_field_data = ref_field.m_field_data;
}

void
HaField3D::FillZeros()
{
	m_field_data = 0.0;
}
void HaField3D::FillValues(float val)
{
	int ix=0,iy=0,iz=0;
	for(ix=0;ix<m_Nx;ix++)
		for(iy=0;iy<m_Ny;iy++)
			for(iz=0;iz<m_Nz;iz++)
	{
		m_field_data[ m_Nx*m_Ny*iz + m_Nx*iy + ix ]=val;
	}
}
void HaField3D::MultiplyByValues(float val)
{
	int ix=0,iy=0,iz=0;
	for(ix=0;ix<m_Nx;ix++)
		for(iy=0;iy<m_Ny;iy++)
			for(iz=0;iz<m_Nz;iz++)
	{
		m_field_data[ m_Nx*m_Ny*iz + m_Nx*iy + ix ]*=val;
	}
}
bool 
HaField3D::SetDimensions(int new_Nx, int new_Ny, int new_Nz)
{
	if( new_Nx < 0 || new_Ny < 0 || new_Nz < 0)
	{
		return false;
	}
	if( new_Nx == 0 || new_Ny == 0 || new_Nz == 0)
	{
		m_Nx = m_Ny = m_Nz = 0;
		m_field_data.newsize(0);
		return true;
	}

//  if no change in dimensions do nothing
	if( m_Nx == new_Nx && m_Ny == new_Ny && m_Nz == new_Nz)
	{
		return true;
	}

	m_Nx = new_Nx;
	m_Ny = new_Ny;
	m_Nz = new_Nz;

	m_field_data.newsize(m_Nx*m_Ny*m_Nz);

	if(m_Nx > 1 )
	{
		xstep = (xmax - xmin)/(m_Nx - 1);
	}

	if(m_Ny > 1 )
	{
		ystep = (ymax - ymin)/(m_Ny - 1);
	}

	if(m_Nz > 1 )
	{
		zstep = (zmax - zmin)/(m_Nz - 1);
	}

	return true;
}

bool HaField3D::SetGridCornersCoord(double xmin_new, double ymin_new, double zmin_new,  
		                            double xmax_new, double ymax_new, double zmax_new)
{
	if( xmin_new + DBL_EPSILON > xmax_new )
	{
		ErrorInMod("HaField3D::SetGridCornersCoord()",
			       " xmin >= xmax_new");
		xmin = 0.0;
		xmax = 1.0;
	}
	else
	{
		xmin = xmin_new;
		xmax = xmax_new;
	}

	if( ymin_new + DBL_EPSILON > ymax_new )
	{
		ErrorInMod("HaField3D::SetGridCornersCoord()",
			       " ymin >= ymax_new");
		ymin = 0.0;
		ymax = 1.0;
	}
	else
	{
		ymin = ymin_new;
		ymax = ymax_new;
	}


	if( zmin_new + DBL_EPSILON > zmax_new )
	{
		ErrorInMod("HaField3D::SetGridCornersCoord()",
			       " zmin >= zmax_new");
		zmin = 0.0;
		zmax = 1.0;
	}
	else
	{
		zmin = zmin_new;
		zmax = zmax_new;
	}

	if(m_Nx > 1 )
	{
		xstep = (xmax - xmin)/(m_Nx - 1);
	}

	if(m_Ny > 1 )
	{
		ystep = (ymax - ymin)/(m_Ny - 1);
	}

	if(m_Nz > 1 )
	{
		zstep = (zmax - zmin)/(m_Nz - 1);
	}

	return true;
}
//>mikola 11-19-06
bool HaField3D::ShiftGridCorners(double xsht, double ysht, double zsht)
{
  double xmin_new=xmin+xsht;
  double ymin_new=ymin+ysht;
  double zmin_new=zmin+zsht;
  double xmax_new=xmax+xsht;
  double ymax_new=ymax+ysht;
  double zmax_new=zmax+zsht;
  return SetGridCornersCoord(xmin_new, ymin_new, zmin_new, xmax_new, ymax_new, zmax_new);
}
bool HaField3D::SetCenterAsZero(float scale)
{
  int nxh=m_Nx/2;
  int nyh=m_Ny/2;
  int nzh=m_Nz/2;
  double xmin_new=-(float)nxh/scale;
  double ymin_new=-(float)nyh/scale;
  double zmin_new=-(float)nzh/scale;
  double xmax_new=(float)nxh/scale;
  double ymax_new=(float)nyh/scale;
  double zmax_new=(float)nzh/scale;
  return SetGridCornersCoord(xmin_new, ymin_new, zmin_new, xmax_new, ymax_new, zmax_new);
}
//<mikola 11-19-06
//>mikola 29,July, 2006
int HaField3D::SaveGZ(const char* fname,int Columns)
{
  TiXmlElement header("HaField3D");
  int GS[3]={m_Nx,m_Ny,m_Nz};
  float Rmin[3]={xmin,ymin,zmin},Rmax[3]={xmax,ymax,zmax};
  HaXML::SetAtributeV(&header,"Rmin",Rmin,3);
  HaXML::SetAtributeV(&header,"Rmax",Rmax,3);
  
  return HaWriteMapGZ(fname,&header,&m_field_data[0],GS,(float)1.0,"");
}

int HaField3D::LoadGZ(const char* fname)
{
  TiXmlElement* header;
  int GS[3]={-1,-1,-1};
  float Rmin[3],Rmax[3];
  float *vfield=NULL;
  HaReadMapGZ(fname,&header,&vfield,GS,1.0);
  
  
  bool bres = SetDimensions(GS[0], GS[1], GS[2]);
  if(!bres) { delete [] vfield; return EXIT_FAILURE; }
  
  if((header->GetArrOfFloatAttribute("Rmin",Rmin,3)==EXIT_SUCCESS)&&(header->GetArrOfFloatAttribute("Rmax",Rmax,3)==EXIT_SUCCESS))
  {
    bres = SetGridCornersCoord(Rmin[0],Rmin[1],Rmin[2],Rmax[0],Rmax[1],Rmax[2]);
    if(!bres) { delete [] vfield; return EXIT_FAILURE; }
  }
  else
  {
    bres = SetCenterAsZero(1.0);
    if(!bres) { delete [] vfield; return EXIT_FAILURE; }
  }
  delete header;

  size_t len = GS[0]*GS[1]*GS[2];
  int i;
  for(i = 0 ; i < len; i++)
  {
    m_field_data[i]=vfield[i];
  }
  delete [] vfield;
  return EXIT_SUCCESS;
}
float HaField3D::CompareHaField3D(HaField3D *CompareWith,float Prec)
{
  int i;
  int GS_XYZ=m_Nx*m_Ny*m_Nz;
  double rmsd=0.0;
  float dif;
  float *orig, *comp;
  if(m_Nx!=CompareWith->GetNx()||m_Nx!=CompareWith->GetNx()||m_Nx!=CompareWith->GetNx())
  {
    fprintf(stdout,"size of HaField3D is not coinside\n");
    return (float) -123.456;
  }
  orig=GetFieldPtr();
  comp=CompareWith->GetFieldPtr();
  if(Prec>0.0f)
    for(i=0;i<GS_XYZ;i++)
  {
    dif=fabs(comp[i]-orig[i]);
    if(dif>Prec)
    {
      fprintf(stdout,"dif>Prec orig[%d]=%g comp[%d]=%g\n",i,orig[i],i,comp[i]);
    }
    rmsd+=dif*dif;
  }
  else
    for(i=0;i<GS_XYZ;i++)
  {
    dif=comp[i]-orig[i];
    rmsd+=dif*dif;
  }
  rmsd=sqrt(rmsd/GS_XYZ);
  fprintf(stdout,"rmsd=%g\n",rmsd);
  return rmsd;
}
//<mikola 29,July, 2006
int 
HaField3D::SaveToFile(const char* fname, int binary)
{
	FILE* fp;
	if( binary )
	{
		fp = fopen(fname,"wb");
	}
	else
	{
		fp = fopen(fname,"w");
	}

	if(fp == NULL)
	{
		PrintLog("Error in HaField3D::SaveToFile() \n");
		PrintLog("Cannot open file for writing \n");
		return FALSE;
	}

	int_4 nx_l = GetNx();
	int_4 ny_l = GetNy();
	int_4 nz_l = GetNz();

	if(binary)
	{
		fwrite(&nx_l,4,1,fp);
		fwrite(&ny_l,4,1,fp);
		fwrite(&nz_l,4,1,fp);
		fwrite(&xmin,4,1,fp);
		fwrite(&ymin,4,1,fp);
		fwrite(&zmin,4,1,fp);
		fwrite(&xmax,4,1,fp);
		fwrite(&ymax,4,1,fp);
		fwrite(&zmax,4,1,fp);
	}
	else
	{
		fprintf(fp," %5d %5d %5d \n",GetNx(),GetNy(),GetNz());
		fprintf(fp," %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f \n",
			      GetXmin(),GetYmin(),GetZmin(),GetXmax(),GetYmax(),GetZmax() );
	}

	int len = GetNx()*GetNy()*GetNz();

	int i;

	if(binary)
	{
		fwrite(&m_field_data[0],4,len,fp);
	}
	else
	{
		for(i = 0 ; i < len; i++)
		{
			fprintf(fp,"%12.6e ",m_field_data[i]);
			if(i % 5 == 0) fprintf(fp," \n"); 
		}
	}
	fclose(fp);
	return TRUE;
}


int 
HaField3D::RestoreFromFile(const char* fname, int binary)
{  
	FILE* fp;
	if( binary )
	{
		fp = fopen(fname,"rb");
	}
	else
	{
		fp = fopen(fname,"r");
	}

	if(fp == NULL)
	{
		PrintLog("Error in HaField3D::RestoreFromFile() \n");
		PrintLog("Unable to open file %s for reading \n",fname);
		return FALSE;
	}

	float xmin_l,ymin_l,zmin_l,xmax_l,ymax_l,zmax_l;
	int_4 nx_l,ny_l,nz_l;

	size_t nread; 
	int ires;

	if(binary)
	{
		nread = fread(&nx_l,4,1,fp);
		nread = fread(&ny_l,4,1,fp);
		nread = fread(&nz_l,4,1,fp);
		nread = fread(&xmin_l,4,1,fp);
		nread = fread(&ymin_l,4,1,fp);
		nread = fread(&zmin_l,4,1,fp);
		nread = fread(&xmax_l,4,1,fp);
		nread = fread(&ymax_l,4,1,fp);
		nread = fread(&zmax_l,4,1,fp);

		if(nread != 1) { fclose(fp); return FALSE; }
	}
	else
	{
		ires = fscanf(fp,"%d %d %d \n",&nx_l,&ny_l,&nz_l);
		ires = fscanf(fp,"%16.9f %16.9f %16.9f %16.9f %16.9f %16.9f \n",
			       &xmin_l,&ymin_l,&zmin_l,&xmax_l,&ymax_l,&zmax_l );

		if(ires == EOF) { fclose(fp); return FALSE; }
	}

	bool bres = SetDimensions(nx_l, ny_l, nz_l);
	if(!bres) { fclose(fp); return FALSE; }

	bres = SetGridCornersCoord(xmin_l,ymin_l,zmin_l,xmax_l,ymax_l,zmax_l);
	if(!bres) { fclose(fp); return FALSE; }

	size_t len = nx_l*ny_l*nz_l;
	int i;

	if(binary)
	{
		nread = fread(&m_field_data[0],4,len,fp);
		if(nread != len) { fclose(fp); return FALSE; }
	}
	else
	{
		for(i = 0 ; i < len; i++)
		{
			ires = fscanf(fp,"%12.6e",&m_field_data[i]);
			if(ires == EOF) { fclose(fp); return FALSE; }
		}
	}
	fclose(fp);
	return TRUE;
}

bool HaField3D::GetXYZ(float& x, float& y , float& z,
		   const int ix, const int iy, const int iz)
{
	if(m_Nx < 2 ||  m_Ny < 2 || m_Nz < 2)
	{
		cerr << " Error in HaField3D::GetXYZ() " << endl;
		cerr << " Some of the 3D Field dimensions are less than 2" << endl;
		x=y=z=0.0;
		return false;
	}

	if( ix >= m_Nx || iy >= m_Ny || iz >= m_Nz )
	{
		cerr << "Error in HaField3D::GetXYZ() " << endl;
		cerr << " One of the indexes out of dimensions of the 3D field " << endl;
		cerr << " ix= " << ix << " iy= " << iy << " iz= " << iz << endl;
		cerr << " nx= " << m_Nx << " Nc= " << m_Ny << " Nl= " << m_Nz << endl;
		x=y=z=0.0;
		return false;
	}

	x = xmin + xstep*ix;
	y = ymin + ystep*iy;
	z = zmin + zstep*iz;

	return true;
}

Vec3D HaField3D::GetGridPtCrd(int ix, int iy, int iz)
{
	Vec3D vec;
	
	if(m_Nx < 2 ||  m_Ny < 2 || m_Nz < 2)
	{
		cerr << " Error in HaField3D::GetGridPtCrd() " << endl;
		cerr << " Some of the 3D Field dimensions are less than 2" << endl;
		vec[0] = 0.0; vec[1] = 0.0;vec[2] = 0.0;
		return vec;
	}

	if( ix >= m_Nx || iy >= m_Ny || iz >= m_Nz )
	{
		cerr << "Error in HaField3D::GetGridPtCrd() " << endl;
		cerr << " One of the indexes out of dimensions of the 3D field " << endl;
		cerr << " ix= " << ix << " iy= " << iy << " iz= " << iz << endl;
		cerr << " nx= " << m_Nx << " Nc= " << m_Ny << " Nl= " << m_Nz << endl;
		vec[0] = 0.0; vec[1] = 0.0;vec[2] = 0.0;
		return vec;
	}

	vec[0] = xmin + xstep*ix;
	vec[1] = ymin + ystep*iy;
	vec[2] = zmin + zstep*iz;

	return vec;
}

void 
HaField3D::GetMinMaxValue(float *ValueMin,float *ValueMax) const
{
	int ix=0,iy=0,iz=0;
	float val=m_field_data[ m_Nx*m_Ny*iz + m_Nx*iy + ix ];
	*ValueMin=val;
	*ValueMax=val;
	
	for(ix=0;ix<m_Nx;ix++)
		for(iy=0;iy<m_Ny;iy++)
			for(iz=0;iz<m_Nz;iz++)
	{
		val=m_field_data[ m_Nx*m_Ny*iz + m_Nx*iy + ix ];
		if(val<*ValueMin)*ValueMin=val;
		if(val>*ValueMax)*ValueMax=val;
	}
}
int 
HaField3D::GetLinIdx(int ix, int iy, int iz)
{
  return ( m_Nx*m_Ny*iz + m_Nx*iy + ix);
}

float* HaField3D::GetValPtr(int ix, int iy, int iz)
{
	if(m_Ny < 2 ||  m_Nx < 2 || m_Nz < 2)
	{
		cerr << " Error in HaField3D::GetValPtr() " << endl;
		cerr << " Some of the 3D Field dimensions less than 2 " << endl;
		return NULL;
	}

	if( ix >= m_Nx || iy >= m_Ny || iz >= m_Nz)
	{
		cerr << "Error in HaField3D::GetValPtr() " << endl;
		cerr << " One of the indexes out of dimensions of the 3D field " << endl;
		return NULL;
	}
	
  return &m_field_data[ m_Nx*m_Ny*iz + m_Nx*iy + ix ];
}

float HaField3D::GetValue( int ix, int iy, int iz)
{
	float* val_ptr = GetValPtr( ix,iy,iz);
	if(val_ptr == NULL)
		return 0.0;
	
	return(*val_ptr);
}

void HaField3D::SetValue( int ix, int iy, int iz, float value)
{
	float* val_ptr = GetValPtr( ix,iy,iz);
	if(val_ptr == NULL)
		return;
    *val_ptr = value;
}



int 
HaField3D::GetClosestGridPoint( double x, double y, double z, int& ix, int& iy, int& iz)
{
	ix = iy = iz = 0;
	if(m_Nx <= 0 || m_Ny <= 0 || m_Nz <= 0 )
	{
		return FALSE;
	}
	
	if( xstep == 0.0 || ystep == 0.0 || zstep == 0.0)
		return FALSE;
	
	ix = (int)((x - xmin)/xstep);
	iy = (int)((y - ymin)/ystep);
	iz = (int)((z - zmin)/zstep);
	
	if(ix < 0 ) ix = 0;
	if(iy < 0 ) iy = 0;
	if(iz < 0 ) iz = 0;

	if( ix >= m_Nx) ix = m_Nx - 1;
	if( iy >= m_Ny) ix = m_Ny - 1;
	if( iz >= m_Nz) ix = m_Nz - 1;

	return TRUE;
}
double HaField3D::CalcLinInter(double V0,double V1,double x)
{
  return V0+x*(V1-V0);
}
double 
HaField3D::GetInterpolValAtPoint( double x, double y, double z)
{
  int ir[3],i;
  double x1[3],fr[3];
  
  ir[0]=(int)((x-xmin)/xstep);
  ir[1]=(int)((y-ymin)/ystep);
  ir[2]=(int)((z-zmin)/zstep);
  fr[0]=(x-xmin)/xstep;
  fr[1]=(y-ymin)/ystep;
  fr[2]=(z-zmin)/zstep;
  //m_Nx, m_Ny, m_Nz;      
  //!@todo don't have good or any extrapolation
  for(i=0;i<3;i++)if(ir[i]<0)ir[i]=0;
  if(ir[0]>m_Nx-2)ir[0]=m_Nx-2;
  if(ir[1]>m_Ny-2)ir[1]=m_Ny-2;
  if(ir[2]>m_Nz-2)ir[2]=m_Nz-2;
  
  for(i=0;i<3;i++)x1[i]=fr[i]-ir[i];
  //printf("r [%f,%f,%f]\n",x,y,z);
  //printf("ri [%d,%d,%d]\n",ir[0],ir[1],ir[2]);
  //printf("rf [%e,%e,%e]\n",fr[0],fr[1],fr[2]);
  //printf("Phi %e\n",GetValue(ir[0],ir[1],ir[2]));
  double vx1,vx2,vx3,vx4,vy1,vy2;
  vx1=CalcLinInter(GetValue(ir[0],ir[1],ir[2]),GetValue(ir[0]+1,ir[1],ir[2]),x1[0]);
  vx2=CalcLinInter(GetValue(ir[0],ir[1]+1,ir[2]),GetValue(ir[0]+1,ir[1]+1,ir[2]),x1[0]);
  vx3=CalcLinInter(GetValue(ir[0],ir[1],ir[2]+1),GetValue(ir[0]+1,ir[1],ir[2]+1),x1[0]);
  vx4=CalcLinInter(GetValue(ir[0],ir[1]+1,ir[2]+1),GetValue(ir[0]+1,ir[1]+1,ir[2]+1),x1[0]);
  vy1=CalcLinInter(vx1,vx2,x1[1]);
  vy2=CalcLinInter(vx3,vx4,x1[1]);
  return CalcLinInter(vy1,vy2,x1[2]);
	/*if( m_Nx == 0 || m_Ny == 0 || m_Nz == 0)
		return 0.0;

	if( xstep == 0.0 || ystep == 0.0 || zstep == 0.0 )
		return 0.0;

	int ix1=1,iy1=1,iz1=1, ix2=1, iy2=1, iz2=1;

	if( x < xmin ) 
	{ 
		ix1 = ix2 = 1;
	}
	else if( x > xmax ) 
	{ 
		ix1 = ix2 = m_Nx;
	}
	else
	{
		ix1 = (int)((x - xmin)/xstep) + 1;
		ix2 = ix1 + 1;
		if(ix1 == m_Nx)
			ix2 = m_Nx;
	}

	if( y < ymin ) 
	{ 
		iy1 = iy2 = 1;
	}
	else if( y > ymax )
	{

		iy1 = iy2 = m_Ny;
	}
	else
	{
		iy1 = (int)((y - ymin)/ystep) + 1;
		iy2 = iy1 + 1;
		if(iy1 == m_Ny)
			iy2 = m_Ny;
	}


	if( z < zmin ) 
	{ 
		iz1 = iz2 = 1;
	}
	else if( z > zmax ) 
	{ 
		iz1 = iz2 = m_Nz;
	}
	else
	{
		iz1 = (int)((z - zmin)/zstep) + 1;
		iz2 = iz1 + 1;
		if(iz1 == m_Nz)
			iz2 = m_Nz;
	}

	double val = ( GetValue( ix1, iy1, iz1) + GetValue( ix2, iy1, iz1) + 
		           GetValue( ix1, iy2, iz1) + GetValue( ix2, iy2, iz1) +
				   GetValue( ix1, iy1, iz2) + GetValue( ix2, iy1, iz2) + 
		           GetValue( ix1, iy2, iz2) + GetValue( ix2, iy2, iz2) ) /8.0;

 return val;*/
}



bool 
HaField3D::grid_to_xyz_float(const int numverts, const float* vr, const float* vc, const float* vl,
							 float* xyz_coord)
// Input:
// numverts   - number of verticies
// vr, vc, vl - grid coordinates -
// Output:
// xyz_coord[numverts][3] - array of cartesian cooerdinates
// 
//
{
	if(m_Ny < 2 ||  m_Nx < 2 || m_Nz < 2)
	{
		cerr << " Error in HaField3D::grid_to_xyz_float " << endl;
		cerr << " m_Ny < 2 ||  m_Nx < 2 || m_Nz < 2 " << endl;
		return false;
	}

	for(int i= 0; i < numverts; i++)
	{
		xyz_coord[3*i] =   xmin + xstep*vc[i];
		xyz_coord[3*i+1] = ymin + ystep*vr[i];
		xyz_coord[3*i+2] = zmin + zstep*vl[i];  
	}
	return true;
}



HaSurface::HaSurface()
{

}

HaSurface::~HaSurface()
{
	clear();
}

void
HaSurface::clear()
{
	verts.newsize(0,0);
	norms.newsize(0,0);
	tr_indx.newsize(0,0);

	isolevel = 0.0;
	valid =0;
}

int
HaSurface::SetNumVerts( const int new_num_verts)
{
	if(new_num_verts < 0)
		return FALSE;
	verts.newsize(3,new_num_verts);
	return TRUE;
}




bool HaSurface::calc_isosurf( HaField3D* field, float iso_level)

//  Input: Nx, Ny, Nz - dimensions of the 3D grid  pointed by float* grid
//         LowLev - Lowest level of 3D grid
//         iso_level - level to construct contour at.
//         colorvar - which variable to color with or -1
//
//  Output: resulting poly-triangle strip is saved in surf
//  
{
   int LowLev = 0;

   /* marching cubes parameters */
   float arx=1.0, ary=1.0, arz=1.0;
   float *vc,  *vr,  *vl;
   float *nx,  *ny,  *nz;
   int *vpts;
   int numverts, numindexes, ipoly, itri;
 
      /* compute the isosurface */

      if (!field) return false;

      vc = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      vr = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      vl = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      nx = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      ny = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      nz = (float *) malloc(sizeof(float)*MAX_ISO_VERTS);
      vpts = (int *) malloc(sizeof(int)*2*MAX_ISO_VERTS);
      if ( !vc || !vr || !vl || !nx || !ny || !nz || !vpts ){
         printf(" You do not have enough memory to create isosurfaces.\n");
         if (vc)
		 {
            free(vc);
         }
         if (vr)
		 {
            free(vr);
         }
         if (vl)
		 {
            free(vl);
         }
         if (nx)
		 {
            free(nx);
         }
         if (ny)
		 {
            free(ny);
         }
         if (nz)
		 {
            free(nz);
         }
         if (vpts)
		 {
            free(vpts);
         }
         return false;
      }
    
      /* Pass number of levels of parameter. main_march is not changed */
      main_march( field->GetFieldPtr(), field->GetNy(), field->GetNx(), field->GetNz(), 
		          LowLev, iso_level, arx, ary, arz, MAX_ISO_VERTS,
                  vr,vc,vl, nx,ny,nz, 2*MAX_ISO_VERTS, vpts,
                  &numverts, &numindexes, &ipoly, &itri );

      if (numindexes>MAX_ISO_VERTS)
         numindexes = MAX_ISO_VERTS;

      /*************************** Compress data ***************************/

      if (numverts>0 && numindexes>0) 
	  {
		  
		  SetNumVerts(numverts);
		  field->grid_to_xyz_float( numverts, vr,vc,vl, verts.begin() );

//		  cerr << endl << " vertices coordinates: " << endl;
//		  for( int j= 1; j <= numverts; j++)
//		  {
//			  cout << " vert# " << j << " grd coord=" << "(" << vr[j-1] << "," << vc[j-1] << "," << vl[j-1] << ")" << endl;
//			  cout << "          cartesian coord = " << "(" << verts(1,j) << "," <<  verts(2,j) << "," << verts(3,j) << ")" << endl;
//		  }
		  
		  norms.newsize(3,numverts);
	      for(int ip = 0; ip < numverts; ip++)
		  {
			 norms(1,ip+1) = ny[ip];
			 norms(2,ip+1) = nx[ip];
			 norms(3,ip+1) = nz[ip];
		  }

//
//		  field->grid_to_xyz_float( numverts, nx,ny,nz, norms.begin());
		  //         project_normals( ctx, numverts, vr,vc,vl, nx,ny,nz, (void*) cnorms );
		  
//		  this->SetNumIndex( numindexes );
//		  memcpy( index.begin(), vpts, numindexes * sizeof(int) );
		  this->isolevel = iso_level;
		  this->valid = 1;

		  bool subst_vert_1= true;
		  int idx1, idx2, idx3;

		  list<int> tr_idx_list;
		  int i;
		  for(i = 1; i <= numindexes; i++ )
		  {
			  idx3 = vpts[i-1];
			  if( (i > 2) &&  idx1 != idx2 &&   idx2 != idx3 && idx1 != idx3 )
			  {
				  tr_idx_list.push_back(idx1);
				  tr_idx_list.push_back(idx2);
				  tr_idx_list.push_back(idx3);
			  }
			  if(subst_vert_1) // New point substitute the first point
			  {
				  idx1 = idx3;
				  subst_vert_1 = false;
			  }
			  else     // New point substitute the second point
			  {
				  idx2 = idx3;
				  subst_vert_1 = true;
			  }       
		  }
		  int ntr = tr_idx_list.size()/3; // the number of triangles
		  tr_indx.newsize(3,ntr);
		  list<int>::iterator itr = tr_idx_list.begin();
		  i = 0;
		  int idx;
		  for(; itr != tr_idx_list.end(); itr++)
		  {
			  i++;
			  idx = (*itr); 
			  tr_indx(1,i) = idx;
			  itr++; 
			  idx = (*itr); 
			  tr_indx(2,i) = idx;
			  itr++; 
			  idx = (*itr); 
			  tr_indx(3,i) = idx;
		  }
	  }
	  else
	  {
		  this->valid = 0;
	  }
      free(vc);
      free(vr);
      free(vl);
      free(nx);
      free(ny);
      free(nz);
      free(vpts);

	  return true;
   
}

bool
HaSurface::Print_info(ostream& sout, const int level) const
{
	int i, nv;
	sout << endl;
	sout << " Description of the surface " << endl;
	nv = GetNumVerts();
	sout << " Number of Verticies " << nv << endl;
	sout << " List of Verticies " << endl;
	for(i = 1; i <= nv; i++)
	{
		sout << i << "  " << verts(1,i) << "  " << verts(2,i) << "  " << verts(3,i) << endl;
	}
	sout << endl;
	return true;
}


HaDisplayedSurface::HaDisplayedSurface():
Object3D(OBJ3D_SURFACE, "Contour_Surface")
{

}

HaDisplayedSurface::~HaDisplayedSurface()
{
	clear();
}

void
HaDisplayedSurface::clear()
{
	HaSurface::clear();
	colors.newsize(0);
}


int
HaDisplayedSurface::SetNumVerts( const int new_num_verts)
{
	int result;
	result = HaSurface::SetNumVerts( new_num_verts );
	if(!result)
		return FALSE;
	if(colors.size() != new_num_verts)
	colors.newsize(new_num_verts);
	colors= 0;
	ColourUniform(255,0,0); // Colour Red by default
	return true;
}


bool HaDisplayedSurface::ColourUniform(int r, int g, int b)
{
	int i;
	HaColor color(r,g,b);

	for(i = 1; i <= colors.size(); i++)
	{
		colors(i) = color.cidx;
	}
	return true;
}


int HaDisplayedSurface::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt)
{
	Object3D::RotateObj( rot_mat, cnt );
	
	int nv = GetNumVerts();

	for(int i = 1; i <= nv; i++)
	{	
		double x = verts(1,i) - cnt.GetX(); 
		double y = verts(2,i) - cnt.GetY();
		double z = verts(3,i) - cnt.GetZ();
						
		verts(1,i) = cnt.GetX() + rot_mat(1,1) * x + rot_mat(1,2) * y + rot_mat(1,3) * z;
		verts(2,i) = cnt.GetY() + rot_mat(2,1) * x + rot_mat(2,2) * y + rot_mat(2,3) * z;
		verts(3,i) = cnt.GetZ() + rot_mat(3,1) * x + rot_mat(3,2) * y + rot_mat(3,3) * z;
	}
	return True;
}


int 
HaDisplayedSurface::Translate( const Vec3D& tr_vec )
{
	bool add_x = ( fabs(tr_vec[0]) > DBL_EPSILON );
	bool add_y = ( fabs(tr_vec[1]) > DBL_EPSILON );
	bool add_z = ( fabs(tr_vec[2]) > DBL_EPSILON );

	int nv = GetNumVerts();

	for(int i = 1; i <= nv; i++ )
	{
		if( add_x ) verts(1,i) += tr_vec[0]; 
		if( add_y ) verts(2,i) += tr_vec[1];
		if( add_z ) verts(3,i) += tr_vec[2];
	}
	
	return True;
}

int 
HaDisplayedSurface::SetTransparency(double transp_new)
{
	transparency = transp_new;
    return TRUE;	
}


HaDot::HaDot(const double new_x, const double new_y,
		  const double new_z, const int new_col)
{
	SetX(new_x);
	SetY(new_y);
	SetZ(new_z);
	col = new_col;
}

HaDot::~HaDot()
{

}

DotStruct::DotStruct():
Object3D(OBJ3D_DOT_SURFACE)
{

}

DotStruct::~DotStruct()
{

}

void
DotStruct::AddDot( double x, double y, double z, int col )
{
	dots.push_back(HaDot(x,y,z,col));
}


int 
DotStruct::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt)
{
	Object3D::RotateObj( rot_mat, cnt );
	
	vector<HaDot>::iterator ditr;;

	for(ditr = dots.begin(); ditr != dots.end(); ditr++)
	{
		double x = (*ditr).GetX() - cnt.GetX(); 
		double y = (*ditr).GetY() - cnt.GetY();
		double z = (*ditr).GetZ() - cnt.GetZ();
						
		(*ditr).SetX( cnt.GetX() + rot_mat(1,1) * x + rot_mat(1,2) * y + rot_mat(1,3) * z);
		(*ditr).SetY( cnt.GetY() + rot_mat(2,1) * x + rot_mat(2,2) * y + rot_mat(2,3) * z);
		(*ditr).SetZ( cnt.GetZ() + rot_mat(3,1) * x + rot_mat(3,2) * y + rot_mat(3,3) * z);
	}
	return True;
}


int 
DotStruct::Translate(const Vec3D& tr_vec )
{	
	bool add_x = ( fabs(tr_vec[0]) > DBL_EPSILON );
	bool add_y = ( fabs(tr_vec[1]) > DBL_EPSILON );
	bool add_z = ( fabs(tr_vec[2]) > DBL_EPSILON );

	vector<HaDot>::iterator ditr;;

	for(ditr = dots.begin(); ditr != dots.end(); ditr++)
	{
		if( add_x ) (*ditr).SetX( (*ditr).GetX() + tr_vec[0]); 
		if( add_y ) (*ditr).SetY( (*ditr).GetY() + tr_vec[1]);
		if( add_z ) (*ditr).SetZ( (*ditr).GetZ() + tr_vec[2]);
	}
	
	return True;
}


HaNonLocField3D_2::HaNonLocField3D_2()
{
	depth = 0;
	nl_cube_size = 0;
}

HaNonLocField3D_2::~HaNonLocField3D_2()
{

}

void 
HaNonLocField3D_2::clear()
{
	HaField3D::clear();
	m_nonl_field_data.newsize(0,0);
}

bool
HaNonLocField3D_2::SetDimensions(int new_Nx, int new_Ny, int new_Nz)
{
	bool same_dim = (new_Nx == m_Nx) && (new_Ny == m_Ny ) && (new_Nz == m_Nz);
	bool bres = HaField3D::SetDimensions(new_Nx, new_Ny, new_Nz);
	if(!bres)
		return false;

	if(depth == 0 || same_dim)
		return true;

	if( depth < 0 )
		depth = -depth;

	if( depth > 10) 
	{
		ErrorInMod("HaNonLocField::SetDimensions()", 
			       " Depth is too large " );
		return false;
	}
		

	m_nonl_field_data.newsize( nl_cube_size, m_Nx*m_Ny*m_Nz );
	 
	return true;
}


bool
HaNonLocField3D_2::SetDepth(int new_depth)
{
	if(depth < 0 ) depth = -depth;
	depth = new_depth;
	nl_cube_size = 2*depth + 1; 
	bool bres = SetDimensions(m_Nx,m_Ny,m_Nz);
	return bres;
}

int
HaNonLocField3D_2::GetDepth() const
{
	return depth;
}

float 
HaNonLocField3D_2::GetValue_nloc(int ix, int iy, int iz,
            		   int ir_shift, int ic_shift, int il_shift)
{
	if(depth == 0)
		return GetValue(ix,iy,iz);

		int indx_pt = m_Nx*m_Ny*iz + m_Ny*iy + ix;
		int idx_shft = nl_cube_size*nl_cube_size* (il_shift + depth)+ 
			           nl_cube_size*(ic_shift + depth) +
					   (depth + ir_shift); 
	return( m_nonl_field_data(idx_shft,indx_pt) );
}


ValAtPoint::ValAtPoint(short ir_new,short ic_new, short il_new, double new_val)
{
	ix = ir_new;
	iy = ic_new;
	iz = il_new;
	val = new_val;
}

HaNonLocField3D::HaNonLocField3D()
{

}

HaNonLocField3D::~HaNonLocField3D()
{

}

int HaNonLocField3D::SaveField(const std::string& fname)
{
	FILE* fp = fopen(fname.c_str(),"w");

	if(fp == NULL)
		return FALSE;

	int nx = GetNx(); 
	int nc = GetNy();
	int nl = GetNz();

	int ngrid_size = nx*nc*nl;

	fprintf(fp," %5d %5d %5d \n", nx,nc,nl);
	fprintf(fp," %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
	 	       xmin,xmax,ymin,ymax,zmin,zmax);
		 

	if(ngrid_size != fvals.size() )
	{
		ErrorInMod("HaNonLocField3D::SaveField()",
		" The number of Grid points doesn't equal to the size of fvals array");
		return FALSE;
	}
	
	int ix, iy, iz;
	int lidx = 0;

	for(iz = 1; iz <= nl ; iz++)
		for(iy = 1; iy <= nc; iy++)
			for( ix = 1; ix <= nx; ix++)
			{
				lidx++;
				list<ValAtPoint>& val_list = fvals[lidx-1];

				fprintf(fp," %5d%5d%5d %8d \n", ix, iy, iz, val_list.size());
				
				list<ValAtPoint>::iterator vitr;
				for(vitr = val_list.begin(); vitr != val_list.end(); vitr++)
				{
					fprintf(fp," %5d%5d%5d %16.9e \n", 
						    (*vitr).ix, (*vitr).iy, (*vitr).iz,(*vitr).val);
				}

			}

	fclose(fp);

	return TRUE;
}

const int jvt1_c[3*60] = { 1,   6,   2,  1,   2,   3,  1,   3,   4,
						   1,   4,   5,  1,   5,   6,  7,   2,   6,
                           8,   3,   2,  9,   4,   3, 10,   5,   4,
                          11,   6,   5,  8,   2,  12,  9,   3,  13,
                          10,   4,  14, 11,   5,  15,  7,   6,  16,
                           7,  12,   2,  8,  13,   3,  9,  14,   4,
                          10,  15,   5, 11,  16,   6,  8,  12,  18,
                           9,  13,  19, 10,  14,  20, 11,  15,  21,
                           7,  16,  17,  7,  17,  12,  8,  18,  13,
                           9,  19,  14, 10,  20,  15, 11,  21,  16,
                          22,  12,  17, 23,  13,  18, 24,  14,  19,
                          25,  15,  20, 26,  16,  21, 22,  18,  12,
                          23,  19,  13, 24,  20,  14, 25,  21,  15,
                          26,  17,  16, 22,  17,  27, 23,  18,  28,
                          24,  19,  29, 25,  20,  30, 26,  21,  31,
                          22,  28,  18, 23,  29,  19, 24,  30,  20,
                          25,  31,  21, 26,  27,  17, 22,  27,  28,
                          23,  28,  29, 24,  29,  30, 25,  30,  31,
                          26,  31,  27, 32,  28,  27, 32,  29,  28,
                          32,  30,  29, 32,  31,  30, 32,  27,  31 };

const int jvt2_c[3*4] = { 1,   5,   4,
                          5,   2,   6,
                          4,   6,   3,
                          6,   4,   5 };

class HaSurface::GEPOLParams HaSurface::gepol_prm; 

HaSurface::GEPOLParams::GEPOLParams()
{
	rmin_axx_sph = 0.5; 
//  rmin_axx_sph = 0.2;
	overlap_axx_sph = (float) 0.8;
//  overlap_axx_sph = 0.95;
//  ndiv_sph = 3;
	ndiv_sph = 2;
}

HaSurface::GEPOLParams::~GEPOLParams()
{

}

int HaSurface::CalcMolSurf( HaSurface* sptr, int surf_type, float solv_rad, AtomContainer& at_coll )
{
	HaMat_float cnt_crd; 
	HaVec_float cnt_rad;
	AtomGroup at_arr;

    AtomIteratorGen aitr(&at_coll);

	int natom = at_coll.GetNAtoms();
	cnt_crd.resize(3,natom);
	cnt_rad.resize(natom);
	at_arr.resize(natom);

	HaAtom* aptr;
	int i = 0;
	logical ghost = 0;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->solv_access_area = 0.0;
		at_arr[i] = aptr;
		
		csfe_.xe[i] = aptr->GetX();
		csfe_.ye[i] = aptr->GetY();
		csfe_.ze[i] = aptr->GetZ();
		double rad  = aptr->radius;

		if(rad > -0.00001 && rad < 0.00001)
		{
			csfe_.iuse[i] = 1;
			rad = 0.0;
		}
		else if( rad > -0.00001 )
		{
			csfe_.iuse[i] = 6;
		}
		else if( rad < -0.00001 )
		{
			csfe_.iuse[i] = 3;
			rad = fabs(rad);
			ghost = 1;
		}
		csfe_.re[i] = rad;
		i++;
	}
	
	int np = 0;
	double svolume = 0.0; 

	if(cnt_rad.size() != natom)
	{
		PrintLog(" Error in HaSurface::CalcMolSurf() \n");
		PrintLog(" Dimensions of coordinate and radii array are not compatible \n"); 
	}

//	penta_.Init();

	for(i=0; i < 180; i++)
	{
		penta_.jvt1[i] = jvt1_c[i];
	}
	for(i=0; i < 12; i++)
	{
		penta_.jvt2[i] = jvt2_c[i];
	}
	logical assign_to_atoms = 1;

	gsurf_(&surf_type,&gepol_prm.rmin_axx_sph,&gepol_prm.overlap_axx_sph,&solv_rad,
		   &gepol_prm.ndiv_sph,&assign_to_atoms,&natom,&ghost,&np,&svolume);

	sptr->surf_volume = svolume;
	PrintLog(" The number of triangles in the surface - %6d \n",np); 

	sptr->SetNumVerts( np*3 );
	sptr->norms.newsize(3,np*3);
	sptr->tr_indx.newsize(3,np);

	for(i = 0; i < np; i++) // cycle on surface triangle tesserae
	{
		int jtr = pun_.ito[i] - 1; // triangle index on a standard sphere

		int idx_sph = pun_.iso[i] - 1;
		float r_sph = csfe_.re[idx_sph];
		float x_sph = csfe_.xe[idx_sph];
		float y_sph = csfe_.ye[idx_sph];
		float z_sph = csfe_.ze[idx_sph];

		int nv[3];

		nv[0] = penta_.jvt1[3*jtr];  // triangle verticies on a standard sphere
		nv[1] = penta_.jvt1[3*jtr + 1];
		nv[2] = penta_.jvt1[3*jtr + 2];

		sptr->tr_indx(1,i+1) = 3*i;
		sptr->tr_indx(2,i+1) = 3*i+1;
		sptr->tr_indx(3,i+1) = 3*i+2;

		int k;
		for(k=0; k <3; k++)
		{
			sptr->verts(1,3*i+k+1) = (poli_.cv[ nv[k] - 1]       * r_sph + x_sph);
			sptr->verts(2,3*i+k+1) = (poli_.cv[ nv[k] - 1 + 32]  * r_sph + y_sph); 
			sptr->verts(3,3*i+k+1) = (poli_.cv[ nv[k] - 1 + 32*2]* r_sph + z_sph);

			sptr->norms(1,3*i+k+1) = poli_.cv[ nv[k] - 1];
			sptr->norms(2,3*i+k+1) = poli_.cv[ nv[k] - 1 + 32]; 
			sptr->norms(3,3*i+k+1) = poli_.cv[ nv[k] - 1 + 32*2]; 
		}
		int idx_at = pun_.isa[i] - 1;
		at_arr[idx_at]->solv_access_area += pun_.ap[i];
	}

	return TRUE;
}

int HaSurface::CalcMolSurfAlpha(int calc_d,double solv_rad, HaMat_double& cnt_crd_alpha, HaVec_double& cnt_rad_alpha)
{
	int natom =  cnt_crd_alpha.num_cols();

	if(cnt_rad_alpha.size() != natom)
	{ 
		ErrorInMod(" Error in HaSurface::CalcMolSurfAlpha() ",
			       " Dimensions of coordinate and radii array are not compatible "); 
		return FALSE;
	} 
	
	surface_alpha.newsize(natom);
	volume_alpha.newsize(natom);
	if (calc_d)
	{
		d_surface_alpha.newsize(3,natom);
		d_volume_alpha.newsize(3,natom);
	}
	int i;

//	GeoBallSurface(&calc_d, &natom, &solv_rad, cnt_crd_alpha.v(),cnt_rad_alpha.v(), surface_alpha.v(), volume_alpha.v(), 
//		           d_surface_alpha.v(), d_volume_alpha.v() );

	surface_alpha_total =0;
	for(i=0; i < natom; i++)
	{ 
		surface_alpha_total = surface_alpha_total + surface_alpha[i] ;
	}  
	return TRUE;
} 

