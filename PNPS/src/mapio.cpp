#ifdef MPI_PARALLEL
#  include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "pnpdebug.h"
#include "mapio.h"
#include "math.h"
#include <iostream>
#include <sstream>

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

#define CHUNK 262144

#ifndef MAP_IO_STRING_LENGTH
#define MAP_IO_STRING_LENGTH 256
#endif
#ifndef MAP_IO_FORCE
#define MAP_IO_FORCE 1
#endif
#ifndef MAP_IO_ADD
#define MAP_IO_ADD 2
#endif

#include "pnps.h"
#include <tinyxml.h>
#include <string>

#ifdef HARLEM_MOD
	#include "hasurface.h"
	#include "halinalg.h"
#endif

using namespace std;

int TypeOfMapFile(const char *filename)
{
	/*!
	Define what tipe of file we have gz, gbn, gb2
	*/
	int i=0,j;
	while(filename[i]!='\0'){
		i++;
		if(i>MAP_IO_STRING_LENGTH)break;
	}
	if(filename[i-1]=='z'&&filename[i-2]=='g'&&filename[i-3]=='.')return FILE_MAP_GZ;
	if(filename[i-1]=='n')return FILE_MAP_MBN;
	if(filename[i-1]=='x'&&filename[i-2]=='d'&&filename[i-3]=='.')return FILE_MAP_DX;
	if(filename[i-1]=='2')return FILE_MAP_MB2;
	if(filename[i-1]=='b')return FILE_MAP_IGB;
	return 0;
}
int ConvertToPBC(float *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ)
{
	if(pbcX==false&&pbcY==false&&pbcZ==false)
		return EXIT_SUCCESS;
	int GS[3]={GS_X, GS_Y, GS_Z};
	int GSwPBC[3]={GS_X+pbcX*2, GS_Y+pbcY*2, GS_Z+pbcZ*2};
	int i,j,k,gridPoint1,gridPoint2;
	int itmp1=GS[0];
	int itmp2=itmp1*GS[1];
	int itmp3=GSwPBC[0];
	int itmp4=itmp3*GSwPBC[1];
	int itmp6,itmp5;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	int GSwPBC_XYZ=GSwPBC[0]*GSwPBC[1]*GSwPBC[2];
	//If periodic Boundary condition then rearrange map
	if(GSwPBC_XYZ!=GS_XYZ)
	{
		//Move data
		for(k=GS[2]-1;k>-1;k--)
			for(j=GS[1]-1;j>-1;j--)
				for(i=GS[0]-1;i>-1;i--)
		{
					gridPoint1 = i + itmp1 * j + itmp2 *k;
					gridPoint2 = i + itmp3 * j + itmp4 *k;
					if(pbcX)gridPoint2+=1;
					if(pbcY)gridPoint2+=itmp3;
					if(pbcZ)gridPoint2+=itmp4;
					V[gridPoint2]=V[gridPoint1];
		}
		//copy sides;
		if(pbcX)
			for(k=0;k<GSwPBC[2];k++)
				for(j=0;j<GSwPBC[1];j++)
				 {
					 itmp1=itmp3 * j + itmp4 * k;
					 itmp2=itmp1 + itmp3 - 2;
					 V[itmp1]=V[itmp2];
					 V[itmp2+1]=V[itmp1+1];
				 }
		if(pbcY)
			for(k=0;k<GSwPBC[2];k++)
				for(i=0;i<GSwPBC[0];i++)
				 {
					 itmp1=i + itmp4 *k;
					 itmp2=itmp1+itmp4 - itmp3;
					 V[itmp1]=V[itmp2 - itmp3];
					 V[itmp2]=V[itmp1 + itmp3];
				 }
		if(pbcZ)
		{
			itmp1=itmp4*(GSwPBC[2]-2);
			itmp2=itmp1+itmp4;
			for(i=0;i<itmp4;i++)
			{
				V[i]=V[i+itmp1];
				V[i+itmp2]=V[i+itmp4];
			}
		}
	}
	return EXIT_SUCCESS;
}
int RemovePBC(float *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ)
{
	int i,j,k,gridPoint1,gridPoint2;
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	int GS_X_PBC=GS_X+pbcX*2;
	int GS_Y_PBC=GS_Y+pbcY*2;
	int GS_Z_PBC=GS_Z+pbcZ*2;
	int GS_XY_PBC=GS_X_PBC*GS_Y_PBC;
	int GS_XYZ_PBC=GS_X_PBC*GS_Y_PBC*GS_Z_PBC;
	//If periodic Boundary condition then rearrange map
	if(GS_XYZ!=GS_XYZ_PBC)
	{
		//Move data
		for(k=0;k<GS_Z;k++)
			for(j=0;j<GS_Y;j++)
				for(i=0;i<GS_X;i++)
		{
			gridPoint1 = i + GS_X * j +  GS_XY*k;
			gridPoint2 = i + pbcX + GS_X_PBC * (j+pbcY) +  GS_XY_PBC*(k+pbcZ);
			V[gridPoint1]=V[gridPoint2];
		}
	}
	return EXIT_SUCCESS;
}
int WriteMapGZinTwoColumns(gzFile file, float * fmap, int GSXYZ, float coef) 
{
	/*!
	@author Warren Dukes 
	@date 2002-10-27
	and I
	*/
	int count;
	float current;
	int gridPoint;
	
	//Write
	count = 0;
	current = fmap[0];
	for(gridPoint=0;gridPoint<GSXYZ;gridPoint++) {
		if(current==fmap[gridPoint])
			count++;
		else {
			if(!gzprintf(file,"%.7e %d\n",current*coef,count)) {
				fprintf(stderr,"ERROR 106: Problems writing to file\n");
				return EXIT_FAILURE;
			}
			count = 1;
			current = fmap[gridPoint];
		}
	}
	if(!gzprintf(file,"%.7e %d\n",current*coef,count)) {
		fprintf(stderr,"ERROR 106: Problems writing to file\n");
		exit(106);
	}
	return EXIT_SUCCESS;
}
int WriteMapGZ(const char *filename, TiXmlElement *header,float * fmap, int* gridsize, float coef,char* Comments) 
{
	//Build XML Header
	header->SetArrOfIntAttribute("GridSize",gridsize,3);
	header->SetAttribute("Comments",Comments);
	header->SetIntAttribute("Columns",2);
	
	gzFile file;
	file = gzopen(filename,"wb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	//header.Print(file);
	int GSXYZ=gridsize[0]*gridsize[1]*gridsize[2];
	WriteMapGZinTwoColumns(file,fmap,GSXYZ,coef);
	gzclose(file);
	return EXIT_SUCCESS;
}
int WriteMapGZ(const char *filename, float * fmap, int* gridsize, float coef,int Columns,char* Comments) 
{
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("Field");
	header->SetArrOfIntAttribute("GridSize", gridsize, 3);
	header->SetAttribute("Comments", Comments);
	header->SetIntAttribute("Columns", 2);
	
	gzFile file;
	file = gzopen(filename,"wb1");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	//header.Print(file);
	int GSXYZ=gridsize[0]*gridsize[1]*gridsize[2];
	if(Columns==1)WriteMapGZOneColumns(file,fmap,GSXYZ,coef);
	else WriteMapGZinTwoColumns(file,fmap,GSXYZ,coef);
	gzclose(file);
	delete header;
	return EXIT_SUCCESS;
}

int ReadMapGZ(gzFile file, float *fmap,int GSXYZ,float coef,int Columns)
{
	int count;
	float current;
	int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int coll=MAP_IO_FORCE;

	if(Columns==1)
	{
		ReadMapGZOneColumns(file,fmap,GSXYZ,coef);
	}
	else
	{
		gridPoint=0;
		while(gridPoint<GSXYZ)
		{
			if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH))
			{
				fprintf(stderr,"readMap: problem reading from file\n");
				return EXIT_FAILURE;
			}
			if(2!=sscanf(str,"%f %i\n",&current,&count)) 
			{
				fprintf(stderr,"readMap: unkown format\n");
				return EXIT_FAILURE;
			}
			//gzscanf(file,"%f %i\n",&current,&count);
			if(gridPoint+count>GSXYZ) {
				fprintf(stderr,"readMap: too many points file\n");
				return EXIT_FAILURE;
			}
			for(i=gridPoint;i<gridPoint+count;i++) {
				if(coll==MAP_IO_ADD)
					fmap[i]+=current*coef;
				else if(coll==MAP_IO_FORCE)
						fmap[i] = current*coef;
				else {
								if(fmap[i]==-1)
							fmap[i] = current*coef;
				}
			}
						gridPoint+=count;
		}
	}
	return EXIT_SUCCESS;
}

int ReadMapGZ(const char *filename, float *fmap,int* gridsize,float coef)
{
	//! Read Map from file and return values *fmap
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GS[3],GSXYZ;
	int Columns;
	
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement* header=new TiXmlElement("Field");
	ins>>*header;
	header->GetArrOfIntAttribute("GridSize",GS,3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=2;
	if(gridsize[0]!=GS[0]||gridsize[1]!=GS[1]||gridsize[2]!=GS[2])
	{
		fprintf(stderr,"readMap: Grid size of map in file do not coinside with request\n");
		return EXIT_FAILURE;
	}
	GSXYZ=GS[0]*GS[1]*GS[2];
	ReadMapGZ(file,fmap,GSXYZ,coef,Columns);
	gzclose(file);
	return EXIT_SUCCESS;
}
int ReadMapGZ(const char *filename,TiXmlElement** pheader, float **pfmap,int* gridsize,float coef)
{
	//! Read Map from file and return values *fmap
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GS[3],GSXYZ;
	int Columns;
	float *vfield=*pfmap;
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("Field");
	*pheader=header;
	ins>>*header;
	header->GetArrOfIntAttribute("GridSize", GS, 3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=2;
	if(gridsize[0]>=0)
	{
		if(gridsize[0]!=GS[0]||gridsize[1]!=GS[1]||gridsize[2]!=GS[2])
		{
			fprintf(stderr,"readMap: Grid size of map in file do not coinside with request\n");
			return EXIT_FAILURE;
		}
	}
	else
	{
		gridsize[0]=GS[0];
		gridsize[1]=GS[1];
		gridsize[2]=GS[2];
	}
	GSXYZ=GS[0]*GS[1]*GS[2];
	
	if(vfield==NULL)vfield=new float[GSXYZ];
	ReadMapGZ(file,vfield,GSXYZ,coef,Columns);
	*pfmap=vfield;
	gzclose(file);
	return EXIT_SUCCESS;
}
int WriteIndexGZinTwoColumns(gzFile file, unsigned	int	* nmap,unsigned	int GSXYZ)
{
	/*!
	@author Warren Dukes 
	@date 2002-10-27
						and me
	 */
				 unsigned int count;
				unsigned	int	current;
				unsigned int gridPoint;
	
	//Write
				count = 0;
				current = nmap[0];
				for(gridPoint=0;gridPoint<GSXYZ;gridPoint++) {
					if(current==nmap[gridPoint])
						count++;
					else {
						if(!gzprintf(file,"%X %d\n",current,count)) {
							fprintf(stderr,"ERROR 106: Problems writing to file\n");
							return EXIT_FAILURE;
						}
						count = 1;
						current = nmap[gridPoint];
					}
				}
				if(!gzprintf(file,"%X %d\n",current,count)) {
					fprintf(stderr,"ERROR 106: Problems writing to file\n");
					exit(106);
				}
				return EXIT_SUCCESS;
}
int WriteIndexGZTwoColumns(gzFile file, int	* nmap,unsigned	int GSXYZ)
{
	/*!
	@author Warren Dukes 
	@date 2002-10-27
						and me
	 */
				unsigned int count;
				int	current;
				unsigned int gridPoint;
	
	//Write
				count = 0;
				current = nmap[0];
				for(gridPoint=0;gridPoint<GSXYZ;gridPoint++) {
					if(current==nmap[gridPoint])
						count++;
					else {
						if(!gzprintf(file,"%d %d\n",current,count)) {
							fprintf(stderr,"ERROR 106: Problems writing to file\n");
							return EXIT_FAILURE;
						}
						count = 1;
						current = nmap[gridPoint];
					}
				}
				if(!gzprintf(file,"%d %d\n",current,count)) {
					fprintf(stderr,"ERROR 106: Problems writing to file\n");
					exit(106);
				}
				return EXIT_SUCCESS;
}
int ReadIndexGZfromTwoColumns(gzFile file,unsigned	int	* nmap,unsigned	int GSXYZ)
{
	int count;
	unsigned	int	current;
	unsigned int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	unsigned int i;
	/*!
	@author Warren Dukes 
	@date 2002-10-27
						and me
	 */
						gridPoint=0;
				while(gridPoint<GSXYZ) {
					if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH))
					{
						fprintf(stderr,"ERR:readMap: problem reading from file\n");
						return EXIT_FAILURE;
					}
					if(2!=sscanf(str,"%X %i\n",&current,&count))
					{
						fprintf(stderr,"ERR:readMap: unkown format\n");
						return EXIT_FAILURE;
					}
					if(gridPoint+count>GSXYZ) {
						fprintf(stderr,"ERR:readMap: too many points file\n");
						return EXIT_FAILURE;
					}
			
					for(i=gridPoint;i<gridPoint+count;i++)nmap[i] = current;
					gridPoint+=count;
				}

				return EXIT_SUCCESS;
}
int ReadIndexGZTwoColumns(gzFile file, int	* nmap,unsigned	int GSXYZ)
{
	int count;
	int	current;
	unsigned int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	unsigned int i;
	/*!
	@author Warren Dukes 
	@date 2002-10-27
						and me
	 */
						gridPoint=0;
				while(gridPoint<GSXYZ) {
					if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH))
					{
						fprintf(stderr,"ERR:readMap: problem reading from file\n");
						return EXIT_FAILURE;
					}
					if(2!=sscanf(str,"%d %d\n",&current,&count))
					{
						fprintf(stderr,"ERR:readMap: unkown format\n");
						return EXIT_FAILURE;
					}
					if(gridPoint+count>GSXYZ) {
						fprintf(stderr,"ERR:readMap: too many points file\n");
						return EXIT_FAILURE;
					}
			
					for(i=gridPoint;i<gridPoint+count;i++)nmap[i] = current;
					gridPoint+=count;
				}

				return EXIT_SUCCESS;
}
int WriteMapGZOneColumns(gzFile file, float	* nmap,unsigned int N,float coef)
{
	unsigned int i;
	float ftmp;
	for(i=0;i<N;i++)
	{
		ftmp=coef*nmap[i];
		//if(ftmp>0.01)
		//	fprintf(stdout,"%.7e\n",ftmp);
		if(!gzprintf(file,"%.7e\n",ftmp)) 
		{
			fprintf(stderr,"ERROR 106: Cannot write to file\n");
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}
int ReadMapGZOneColumns(gzFile file, float	* nmap,unsigned	int N,float coef)
{
	float current;
	unsigned int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	unsigned int i;
	
	gridPoint=0;
	while(gzgets(file,str,MAP_IO_STRING_LENGTH)!=Z_NULL&&gridPoint<N)
	{
		nmap[gridPoint] = coef*atof(str);
		gridPoint++;
	}

	return EXIT_SUCCESS;
}
int WriteMapGZTwoColumns(gzFile file, float	* nmap,unsigned int N,float coef)
{
	return WriteMapGZinTwoColumns(file,nmap,N,coef);
}
int ReadMapGZTwoColumns(gzFile file, float	* nmap,unsigned	int N,float coef)
{
	return ReadMapGZ(file,nmap,N,coef,2);
}
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////VectorField3D
VectorField3D::VectorField3D()
{
	amode=INTERNAL_ALLOC;
	InitZero();
}
VectorField3D::VectorField3D(int *gridsize, float gridscale, int nelem)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	SetVectorField3D(gridsize,gridscale,nelem);
}
VectorField3D::VectorField3D(const char* filename,float coef)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	ReadFromFile(filename,coef);
}
VectorField3D::VectorField3D(int Nx,int Ny,int Nz, float gridscale, int nelem)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	int gridsize[3]={Nx, Ny, Nz};
	SetVectorField3D(gridsize,gridscale,nelem);
}
VectorField3D::VectorField3D(int *gridsize, float gridscale, int nelem,float **v)
{
	amode=EXTERNAL_ALLOC;
	InitZero();
	Clear();
	
	GridSize[0]=gridsize[0];
	GridSize[1]=gridsize[1];
	GridSize[2]=gridsize[2];
	GridScale=gridscale;
	Nelem=nelem;
	V=v;
}

VectorField3D::~VectorField3D()
{
	if(amode==INTERNAL_ALLOC)Clear();
}
int VectorField3D::InitZero()
{
	GridSize[0]=0;
	GridSize[1]=0;
	GridSize[2]=0;
	GridScale=0.0f;
	Nelem=0;
	V=NULL;
	return EXIT_SUCCESS;
}
int VectorField3D::Clear()
{
	int i;
	
	if(V!=NULL)
	{
		for(i=0;i<Nelem;i++)
		{
			if(V[i]!=NULL)delete [] V[i];
		}
		delete [] V;
	}
	InitZero();
	return EXIT_SUCCESS;
}
int VectorField3D::SetVectorField3D(int *gridsize, float gridscale, int nelem)
{
	Clear();
	
	GridSize[0]=gridsize[0];
	GridSize[1]=gridsize[1];
	GridSize[2]=gridsize[2];
	GridScale=gridscale;
	Nelem=nelem;
	int i,j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	V=new float*[Nelem];
	for(i=0;i<Nelem;i++)
	{
		V[i]=new float[GS_XYZ];
		for(j=0;j<GS_XYZ;j++)
		{
			 V[i][j]=0.0f;
		}
	}
	return EXIT_SUCCESS;
}

float VectorField3D::GetInterpolatedValueGrid(int in, float fx, float fy, float fz)
{
	int ir[3],i;
	float x1[3],fr[3];
	ir[0]=fx;
	ir[1]=fy;
	ir[2]=fz;
	fr[0]=fx;
	fr[1]=fy;
	fr[2]=fz;
	//!@todo don't have good or any extrapolation
	for(i=0;i<3;i++)if(ir[i]<0)ir[i]=0;
	if(ir[0]>GridSize[0]-2)ir[0]=GridSize[0]-2;
	if(ir[1]>GridSize[1]-2)ir[1]=GridSize[1]-2;
	if(ir[2]>GridSize[2]-2)ir[2]=GridSize[2]-2;
	
	for(i=0;i<3;i++)x1[i]=fr[i]-ir[i];
	float vx1,vx2,vx3,vx4,vy1,vy2;
	vx1=CalcLinInterFloat(GetValue(in,ir[0],ir[1],ir[2]),GetValue(in,ir[0]+1,ir[1],ir[2]),x1[0]);
	vx2=CalcLinInterFloat(GetValue(in,ir[0],ir[1]+1,ir[2]),GetValue(in,ir[0]+1,ir[1]+1,ir[2]),x1[0]);
	vx3=CalcLinInterFloat(GetValue(in,ir[0],ir[1],ir[2]+1),GetValue(in,ir[0]+1,ir[1],ir[2]+1),x1[0]);
	vx4=CalcLinInterFloat(GetValue(in,ir[0],ir[1]+1,ir[2]+1),GetValue(in,ir[0]+1,ir[1]+1,ir[2]+1),x1[0]);
	vy1=CalcLinInterFloat(vx1,vx2,x1[1]);
	vy2=CalcLinInterFloat(vx3,vx4,x1[1]);
	return CalcLinInterFloat(vy1,vy2,x1[2]);
}
int VectorField3D::InterpolateFromExt(VectorField3D* VExt)
{
	int i;
	int CentralNode[3];
	int CentralNodeExt[3];
	float fxExt,fyExt,fzExt;
	float fx,fy,fz;
	float GSr=VExt->GridScale/GridScale;
	int Start[3];
	int End[3];
	
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=VExt->GridSize[i]/2;
		Start[i]=0;
		End[i]=GridSize[i]-1;
		fxExt=(Start[i]-CentralNode[i])*GSr+CentralNodeExt[i];
		if(fxExt<0.0)
		{
			fx=(0-CentralNodeExt[i])/GSr+CentralNode[i];
			Start[i]=int(fx)+2;
		}
		fxExt=(End[i]-CentralNode[i])*GSr+CentralNodeExt[i];
		if(fxExt>VExt->GridSize[i]-1)
		{
			fx=(VExt->GridSize[i]-1-CentralNodeExt[i])/GSr+CentralNode[i];
			End[i]=int(fx)-2;
		}
	}
	if(Start[0]!=0||Start[1]!=0||Start[2]!=0||End[0]!=GridSize[0]-1||End[1]!=GridSize[1]-1||End[2]!=GridSize[2]-1)
		pnpPrint("Will interpolate internal part: Start %d %d %d End %d %d %d\n",Start[0],Start[1],Start[2],End[0],End[1],End[2]);
	
	int in, ix,iy,iz,GrdPnt;
	for(in=0;in<VExt->Nelem;in++)
		for(ix=Start[0];ix<=End[0];ix++)
			for(iy=Start[1];iy<=End[1];iy++)
				for(iz=Start[2];iz<=End[2];iz++)
	{
		GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
		fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
		fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
		fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
		V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
	}
	return EXIT_SUCCESS;
}
int VectorField3D::InterpolateInternalBoarderFromExt(VectorField3D* VExt,int *InternalBoarderMinGridLev0,int *InternalBoarderMaxGridLev0)
{
	int i;
	int CentralNode[3];
	int CentralNodeExt[3];
	float fxExt,fyExt,fzExt;
	float GSr=VExt->GridScale/GridScale;
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=VExt->GridSize[i]/2;
	}
	
	int in, ix,iy,iz,GrdPnt;
	for(in=0;in<VExt->Nelem;in++)
	{
		//
		iz=InternalBoarderMinGridLev0[2];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iz=InternalBoarderMaxGridLev0[2];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=InternalBoarderMinGridLev0[1];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=InternalBoarderMaxGridLev0[1];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=InternalBoarderMinGridLev0[0];
		for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=InternalBoarderMaxGridLev0[0];
		for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::InterpolateInternalBoarderFromExtWithMultipl(VectorField3D* VExt,int *InternalBoarderMinGridLev0,int *InternalBoarderMaxGridLev0,float Multipl)
{
	int i;
	int CentralNode[3];
	int CentralNodeExt[3];
	float fxExt,fyExt,fzExt;
	float GSr=VExt->GridScale/GridScale;
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=VExt->GridSize[i]/2;
	}
	
	int in, ix,iy,iz,GrdPnt;
	for(in=0;in<VExt->Nelem;in++)
	{
		//
		iz=InternalBoarderMinGridLev0[2];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iz=InternalBoarderMaxGridLev0[2];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=InternalBoarderMinGridLev0[1];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=InternalBoarderMaxGridLev0[1];
		for(ix=InternalBoarderMinGridLev0[0];ix<=InternalBoarderMaxGridLev0[0];ix++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=InternalBoarderMinGridLev0[0];
		for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=InternalBoarderMaxGridLev0[0];
		for(iy=InternalBoarderMinGridLev0[1];iy<=InternalBoarderMaxGridLev0[0];iy++)
			for(iz=InternalBoarderMinGridLev0[2];iz<=InternalBoarderMaxGridLev0[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::InterpolateBoarderFromExt(VectorField3D* VExt)
{
	int i;
	int CentralNode[3];
	int CentralNodeExt[3];
	float fxExt,fyExt,fzExt;
	float GSr=VExt->GridScale/GridScale;
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=VExt->GridSize[i]/2;
	}
	/*pnpPrint("CentralNode = [%d, %d, %d] [g, g, g]\n",
					 CentralNode[0],CentralNode[1],CentralNode[2]);
	pnpPrint("CentralNodeExt = [%d, %d, %d] [g, g, g]\n",
					 CentralNodeExt[0],CentralNodeExt[1],CentralNodeExt[2]);
	fxExt=(0-CentralNode[0])*GSr+CentralNodeExt[0];
	fyExt=(0-CentralNode[1])*GSr+CentralNodeExt[1];
	fzExt=(0-CentralNode[2])*GSr+CentralNodeExt[2];
	pnpPrint("Ext0 = [%g, %g, %g] [g, g, g]\n",fxExt,fyExt,fzExt);
	fxExt=(GridSize[0]-1-CentralNode[0])*GSr+CentralNodeExt[0];
	fyExt=(GridSize[1]-1-CentralNode[1])*GSr+CentralNodeExt[1];
	fzExt=(GridSize[2]-1-CentralNode[2])*GSr+CentralNodeExt[2];
	pnpPrint("Ext1 = [%g, %g, %g] [g, g, g]\n",fxExt,fyExt,fzExt);*/
	
	int in, ix,iy,iz,GrdPnt;
	for(in=0;in<VExt->Nelem;in++)
	{
		//
		iz=0;
		for(ix=0;ix<GridSize[0];ix++)
			for(iy=0;iy<GridSize[1];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iz=GridSize[2]-1;
		for(ix=0;ix<GridSize[0];ix++)
			for(iy=0;iy<GridSize[1];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=0;
		for(ix=0;ix<GridSize[0];ix++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=GridSize[1]-1;
		for(ix=0;ix<GridSize[0];ix++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=0;
		for(iy=0;iy<GridSize[1];iy++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=GridSize[0]-1;
		for(iy=0;iy<GridSize[1];iy++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::InterpolateBoarderFromExtWithMultipl(VectorField3D* VExt,float Multipl)
{
	int i;
	int CentralNode[3];
	int CentralNodeExt[3];
	float fxExt,fyExt,fzExt;
	float GSr=VExt->GridScale/GridScale;
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=VExt->GridSize[i]/2;
	}
	
	int in, ix,iy,iz,GrdPnt;
	for(in=0;in<VExt->Nelem;in++)
	{
		//
		iz=0;
		for(ix=0;ix<GridSize[0];ix++)
			for(iy=0;iy<GridSize[1];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iz=GridSize[2]-1;
		for(ix=0;ix<GridSize[0];ix++)
			for(iy=0;iy<GridSize[1];iy++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=0;
		for(ix=0;ix<GridSize[0];ix++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		iy=GridSize[1]-1;
		for(ix=0;ix<GridSize[0];ix++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=0;
		for(iy=0;iy<GridSize[1];iy++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
		ix=GridSize[0]-1;
		for(iy=0;iy<GridSize[1];iy++)
			for(iz=0;iz<GridSize[2];iz++)
		{
			GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			fxExt=(ix-CentralNode[0])*GSr+CentralNodeExt[0];
			fyExt=(iy-CentralNode[1])*GSr+CentralNodeExt[1];
			fzExt=(iz-CentralNode[2])*GSr+CentralNodeExt[2];
			V[in][GrdPnt]=Multipl*VExt->GetInterpolatedValueGrid(in,fxExt,fyExt,fzExt);
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileAddPBCandSplitSDMPI4FDaZ(const char *filename, float coef, bool *pbc,int LocaliZ0,int LocaliZ1)
{
#ifdef MPI_PARALLEL
	VectorField3D *VExt=NULL;
	if(pnpsapp->GetMyRankInGroup()==pnpsapp->GetMyGroupLeader())
	{
		VExt = new VectorField3D();
		VExt->ReadFromFileAddPBC(filename,coef,pbc[0],pbc[1],pbc[2]);
		
		DbgPrint0("ReadFromFileAddPBCandSplitSDMPI4FDaZ Read from file GS %d %d %d GScale %f Nelem %d\n",VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2],VExt->GridScale,VExt->Nelem);
	}
	SplitExtVF3DSDMPI4FDaZ(VExt,LocaliZ0,LocaliZ1);
	DeleteObjByPnt(VExt);
#endif
	return EXIT_SUCCESS;
}
int VectorField3D::SplitExtVF3DSDMPI4FDaZ(VectorField3D *VExt,int LocaliZ0,int LocaliZ1)
{
#ifdef MPI_PARALLEL
	int i,elm,GS_XY,GS_XYZ,proc;
	int GS[3],nelem;
	float gridscale;
	
	if(pnpsapp->GetMyRankInGroup()==pnpsapp->GetMyGroupLeader())
	{
		gridscale=VExt->GridScale;
		nelem=VExt->Nelem;
		
		//itself
		GS[0]=VExt->GridSize[0];
		GS[1]=VExt->GridSize[1];
		
		GS[2]=LocaliZ1-LocaliZ0+1;
		
		SetVectorField3D(GS,gridscale,nelem);
		GS_XY=this->GridSize[0]*this->GridSize[1];
		GS_XYZ=this->GridSize[0]*this->GridSize[1]*this->GridSize[2];
		for(elm=0;elm<this->Nelem;elm++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				this->V[elm][i]=VExt->V[elm][i];
			}
		}
		DbgPrint0("VExt->GridSize=%d %d %d\n", VExt->GridSize[0], VExt->GridSize[1], VExt->GridSize[2]);
		DbgPrint0("ReadFromFileAddPBCandSplitSDMPI4FDaZ GS %d %d %d GScale %f Nelem %d LocaliZ0 %d LocaliZ1 %d\n",GridSize[0],GridSize[1],GridSize[2],GridScale,Nelem,LocaliZ0,LocaliZ1);
		//others proseses
		for(proc=1;proc<pnpsapp->GetNumProcsInMyGroup();proc++)
		{
			pnpsapp->MyComGroup.Send(&VExt->GridScale, 1, MPI::FLOAT, proc, 1);
			pnpsapp->MyComGroup.Send(&VExt->Nelem, 1, MPI::INT, proc, 2);
			pnpsapp->MyComGroup.Send(&VExt->GridSize, 3, MPI::INT, proc, 3);
			
			pnpsapp->MyComGroup.Recv(GS, 3, MPI::INT, proc, 4);
			DbgPrint0("From proc %d LocaliZ0=%d LocaliZ1=%d GS[2]=%d\n",proc,GS[0],GS[1],GS[2]);
			for(elm=0;elm<this->Nelem;elm++)
			{
				DbgPrint0("R %d\n",elm);
				pnpsapp->MyComGroup.Send(VExt->V[elm]+GS[0]*GS_XY, GS[2]*GS_XY, MPI::FLOAT, proc, elm);
			}
		}
		
	}
	else
	{
		pnpsapp->MyComGroup.Recv(&gridscale, 1, MPI::FLOAT, 0, 1);
		pnpsapp->MyComGroup.Recv(&nelem, 1, MPI::INT, 0, 2);
		pnpsapp->MyComGroup.Recv(GS, 3, MPI::INT, 0, 3);
		
		GS[2]=LocaliZ1-LocaliZ0+1;
		
		SetVectorField3D(GS,gridscale,nelem);
		GS_XY=this->GridSize[0]*this->GridSize[1];
		GS_XYZ=this->GridSize[0]*this->GridSize[1]*this->GridSize[2];
		
		DbgPrint0("ReadFromFileAddPBCandSplitSDMPI4FDaZ GS %d %d %d GScale %f Nelem %d LocaliZ0 %d LocaliZ1 %d\n",GridSize[0],GridSize[1],GridSize[2],GridScale,Nelem,LocaliZ0,LocaliZ1);
		GS[0]=LocaliZ0;
		GS[1]=LocaliZ1;
		
		DbgPrint0("R\n");
		pnpsapp->MyComGroup.Send(GS, 3, MPI::INT, 0, 4);
		
		for(elm=0;elm<Nelem;elm++)
		{
			DbgPrint0("R %d\n",elm);
			pnpsapp->MyComGroup.Recv(V[elm], GS[2]*GS_XY, MPI::FLOAT, 0, elm);
		}
	}
#endif
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileAddPBC(const char *filename, float coef, bool pbcX, bool pbcY, bool pbcZ)
{
	if(pbcX==false&&pbcY==false&&pbcZ==false)
		return ReadFromFile(filename,coef);
	else
	{
		if(TypeOfMapFile(filename)==FILE_MAP_GZ)
			return ReadFromFileAddPBC_GZ(filename, coef, pbcX, pbcY, pbcZ);
		else
			return ReadFromFileAddPBC_BIN(filename, coef, pbcX, pbcY, pbcZ);
	}
}
int VectorField3D::ReadFromFile(const char *filename,float coef)
{
	if(TypeOfMapFile(filename)==FILE_MAP_GZ)
		return ReadFromFileGZ(filename, coef);
	else if(TypeOfMapFile(filename)==FILE_MAP_MBN)
		return ReadFromFileBIN(filename, coef);
	else if(TypeOfMapFile(filename)==FILE_MAP_DX)
		return ReadFromFileDX(filename, coef);
}
int VectorField3D::ReadFromFileDX(const char *filename,float coef)
{
	FILE *file;
	int i,j;
	int GS[3],nelem=1;
	float gridscale;
	
	file = fopen(filename,"rt");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	/*# PME potential (kT/e, T=300K)
object 1 class gridpositions counts 176 176 144
origin 0.473442 0.603447 0.463593
delta 0.968875 0 0
delta -0 0.968875 0
delta -0 -0 0.948201
object 2 class gridconnections counts 176 176 144
object 3 class array type double rank 0 items 4460544 data follows
	*/

	char str[MAP_IO_STRING_LENGTH];
	//comment
	fgets(str,MAP_IO_STRING_LENGTH,file);
	//object 1 class gridpositions counts 176 176 144
	fgets(str,MAP_IO_STRING_LENGTH,file);
	//origin
	fgets(str,MAP_IO_STRING_LENGTH,file);
	//delta x
	fgets(str,MAP_IO_STRING_LENGTH,file);
	sscanf(str,"%*s %f",&gridscale);
	//delta y
	fgets(str,MAP_IO_STRING_LENGTH,file);
	//delta z
	fgets(str,MAP_IO_STRING_LENGTH,file);
	//object 2 class gridconnections counts 176 176 144
	fgets(str,MAP_IO_STRING_LENGTH,file);
	sscanf(str,"%*s %*d %*s %*s %*s %d %d %d",&(GS[0]),&(GS[1]),&(GS[2]));
	//object 3 class array type double rank 0 items 4460544 data follows
	fgets(str,MAP_IO_STRING_LENGTH,file);
	
	DbgPrint0("GS [%d %d %d] %f\n",GS[0],GS[1],GS[2],gridscale);
	
	SetVectorField3D(GS,gridscale,nelem);
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	for(i=0;i<GS_XYZ;i++)
	{
		fscanf(file,"%f",&(V[0][i]));
	}
	fclose(file);
	DbgPrint0("done\n");
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileBIN(const char *filename,float coef)
{
	FILE *file;
	int i,j;
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = fopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	
	
	int size;
	fread(&size,sizeof(int),1,file);
	DbgPrint0("Size %d %d\n",size,sizeof(int));
	char *str;
	str=new char[size+1];
	fread(str,sizeof(char),size,file);
	str[size]='\0';
	
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GS,3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	delete [] str;
	SetVectorField3D(GS,gridscale,nelem);
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	for(i=0;i<Nelem;i++)
	{
		fread(V[i],sizeof(float),GS_XYZ,file);
		for(j=0;j<GS_XYZ;j++)
		{
			V[i][j]*=coef;
		}
	}
	fclose(file);
	
	return EXIT_SUCCESS;
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileGZ(const char *filename,float coef)
{
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GS,3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	
	SetVectorField3D(GS,gridscale,nelem);
	
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *VPar;
	for(i=0;i<Nelem;i++)
	{
		VPar=V[i];
		if(Columns==1)
			ReadMapGZOneColumns(file,V[i],GS_XYZ,coef);
		else
			ReadMapGZTwoColumns(file,V[i],GS_XYZ,coef);
	}
	gzclose(file);
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileAddPBC_BIN(const char *filename,float coef, bool pbcX, bool pbcY, bool pbcZ)
{
	FILE *file;
	int i,j;
	int GSorig[3];
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = fopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	
	
	int size;
	fread(&size,sizeof(int),1,file);
	DbgPrint0("Size %d %d\n",size,sizeof(int));
	char *str;
	str=new char[size+1];
	fread(str,sizeof(char),size,file);
	str[size]='\0';
	
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GSorig,3);
	GS[0]=GSorig[0]+2*pbcX;
	GS[1]=GSorig[1]+2*pbcY;
	GS[2]=GSorig[2]+2*pbcZ;
	
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	delete [] str;
	SetVectorField3D(GS,gridscale,nelem);
	int GS_XYZ=GS[0]*GS[1]*GS[2];
	int GSorig_XYZ=GSorig[0]*GSorig[1]*GSorig[2];
	for(i=0;i<Nelem;i++)
	{
		fread(V[i],sizeof(float),GSorig[0]*GSorig[1]*GSorig[2],file);
		for(j=0;j<GSorig_XYZ;j++)
		{
			V[i][j]*=coef;
		}
		ConvertToPBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
	}
	fclose(file);
	
	return EXIT_SUCCESS;
}
int VectorField3D::ReadFromFileAddPBC_GZ(const char *filename,float coef, bool pbcX, bool pbcY, bool pbcZ)
{
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GSorig[3];
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GSorig,3);
	GS[0]=GSorig[0]+2*pbcX;
	GS[1]=GSorig[1]+2*pbcY;
	GS[2]=GSorig[2]+2*pbcZ;
	
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	
	SetVectorField3D(GS,gridscale,nelem);
	
	int GSorig_XYZ=GSorig[0]*GSorig[1]*GSorig[2];
	for(i=0;i<Nelem;i++)
	{
		if(Columns==1)
			ReadMapGZOneColumns(file,V[i],GSorig_XYZ,coef);
		else
			ReadMapGZTwoColumns(file,V[i],GSorig_XYZ,coef);
		ConvertToPBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
	}
	
	gzclose(file);
	return EXIT_SUCCESS;
}
int VectorField3D::WriteToFileCombineSDMPIDistrRemovePBC(const char *filename,float coef,int Columns, int *GridSizeGlobal,bool *pbc,int LocaliZ0,int LocaliZ1)
{
#ifdef MPI_PARALLEL
	int i,elm,GS_XY,GS_XYZ,proc;
	int GS[3];
	
	if(pnpsapp->GetMyRankInGroup()==pnpsapp->GetMyGroupLeader())
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSizeGlobal,GridScale,Nelem);
		
		GS_XY=GridSize[0]*GridSize[1];
		
		GS[0]=LocaliZ0;
		GS[1]=LocaliZ1;
		GS[2]=LocaliZ1-LocaliZ0+1;
		GS_XYZ=GridSize[0]*GridSize[1]*GS[2];
		
		//itself
		for(elm=0;elm<this->Nelem;elm++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				VField3DContr->V[elm][i]=this->V[elm][i];
			}
		}
		
		for(proc=1;proc<pnpsapp->GetNumProcsInMyGroup();proc++)
		{
			pnpsapp->MyComGroup.Recv(GS, 3, MPI::INT, proc, 4);
			
			for(elm=0;elm<this->Nelem;elm++)
			{
				pnpsapp->MyComGroup.Recv(VField3DContr->V[elm]+GS[0]*GS_XY, GS[2]*GS_XY, MPI::FLOAT, proc, elm);
			}
		}
		
		VField3DContr->WriteToFileRemovePBC(filename,coef,Columns,pbc[0],pbc[1],pbc[2]);
		delete VField3DContr;VField3DContr=NULL;
	}
	else
	{
		GS[0]=LocaliZ0;
		GS[1]=LocaliZ1;
		GS[2]=LocaliZ1-LocaliZ0+1;
		GS_XY=GridSize[0]*GridSize[1];
		
		pnpsapp->MyComGroup.Send(GS, 3, MPI::INT, 0, 4);
		
		for(elm=0;elm<Nelem;elm++)
		{
			pnpsapp->MyComGroup.Send(V[elm], GS[2]*GS_XY, MPI::FLOAT, 0, elm);
		}
	}
#endif
	return EXIT_SUCCESS;
}

int VectorField3D::WriteToFileRemovePBC(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ)
{
	if(pbcX==false&&pbcY==false&&pbcZ==false)
		return WriteToFile(filename,coef,Columns);
	else
	{
		if(TypeOfMapFile(filename)==FILE_MAP_GZ)
			return WriteToFileRemovePBC_GZ(filename, coef, Columns, pbcX, pbcY, pbcZ);
		else
			return WriteToFileRemovePBC_BIN(filename, coef, Columns, pbcX, pbcY, pbcZ);
	}
}
int VectorField3D::WriteToFileRemovePBC_GZ(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	int GSorig[3]={GridSize[0],GridSize[1],GridSize[2]};
	int GS[3]={GridSize[0],GridSize[1],GridSize[2]};
	
	if(pbcX)
		GSorig[0]=GSorig[0]-2;
	if(pbcY)
		GSorig[1]=GSorig[1]-2;
	if(pbcY)
		GSorig[2]=GSorig[2]-2;
	
	int GS_XYZ=GS[0]*GS[1]*GS[2];
	int GSorig_XYZ=GSorig[0]*GSorig[1]*GSorig[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GSorig,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	header->SetIntAttribute("Columns",Columns);
	
	gzFile file;
	file = gzopen(filename,"wb1");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	
	for(i=0;i<Nelem;i++)
	{
		RemovePBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
		if(Columns==1)
			WriteMapGZOneColumns(file,V[i],GSorig_XYZ,coef);
		else
			WriteMapGZTwoColumns(file,V[i],GSorig_XYZ,coef);
		ConvertToPBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
	}
	gzclose(file);
	delete header;
	return EXIT_SUCCESS;
}
int VectorField3D::WriteToFileRemovePBC_BIN(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	int GSorig[3] = { GridSize[0] - 2 * pbcX,GridSize[1] - 2 * pbcY,GridSize[2] - 2 * pbcZ };
	int GS[3]={GridSize[0],GridSize[1],GridSize[2]};
	
	int GS_XYZ=GS[0]*GS[1]*GS[2];
	int GSorig_XYZ=GSorig[0]*GSorig[1]*GSorig[2];
	
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GSorig,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	
	FILE *file;
	file = fopen(filename,"wb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	int size=StrHeader.size();
	fwrite(&size,sizeof(int),1,file);
	fwrite(StrHeader.c_str(),sizeof(char),StrHeader.size(),file);
	
	MultiplyBy(coef);
	for(i=0;i<Nelem;i++)
	{
		RemovePBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
		fwrite(V[i],sizeof(float),GSorig_XYZ,file);
		ConvertToPBC(V[i], GSorig[0], GSorig[1], GSorig[2], pbcX, pbcY, pbcZ);
	}
	MultiplyBy(1.0/coef);
	fclose(file);
	delete header;
	return EXIT_SUCCESS;
}
int VectorField3D::WriteToFile(const char *filename,float coef,int Columns)
{
	if(TypeOfMapFile(filename)==FILE_MAP_GZ)return WriteToFileGZ(filename, coef, Columns);
	else return WriteToFileBIN(filename, coef, Columns);
}
int VectorField3D::WriteToFileBIN(const char *filename,float coef,int Columns)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GridSize,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	
	FILE *file;
	file = fopen(filename,"wb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	int size=StrHeader.size();
	fwrite(&size,sizeof(int),1,file);
	fwrite(StrHeader.c_str(),sizeof(char),StrHeader.size(),file);
	
	MultiplyBy(coef);
	for(i=0;i<Nelem;i++)
	{
		fwrite(V[i],sizeof(float),GS_XYZ,file);
	}
	MultiplyBy(1.0/coef);
	fclose(file);
	delete header;
	return EXIT_SUCCESS;
}
int VectorField3D::WriteToFileGZ(const char *filename,float coef,int Columns)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GridSize,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	header->SetIntAttribute("Columns",Columns);
	
	gzFile file;
	file = gzopen(filename,"wb1");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	
	for(i=0;i<Nelem;i++)
	{
		if(Columns==1)
			WriteMapGZOneColumns(file,V[i],GS_XYZ,coef);
		else
			WriteMapGZTwoColumns(file,V[i],GS_XYZ,coef);
	}
	gzclose(file);
	delete header;
	return EXIT_SUCCESS;
}
int VectorField3D::MultiplyBy(float coef)
{
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)V[i][j]*=coef;
	}
	return EXIT_SUCCESS;
}
int VectorField3D::MultiplyOneElementBy(float coef,int Elm2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot multiply since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(j=0;j<GS_XYZ;j++)V[Elm2set][j]*=coef;
	return EXIT_SUCCESS;
}
int VectorField3D::MaskWithVectorField3D(VectorField3D* VExt,float Value, float dvalue,float value2set)
{
	if(VExt->Nelem!=Nelem)
	{
		pnpError("Cannot mask since number of elements is not coinside (this->Nelem=%d)!=(VExt->Nelem=%d)\n",Nelem,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot mask since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)
		{
			if(fabs(VExt->V[i][j]-Value)<=dvalue)
				V[i][j]=value2set;
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::MaskWithVectorField3DUseOneElement(VectorField3D* VExt,int VExtElm,int Elm2set,float Value, float dvalue,float value2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot multiply since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	if(VExtElm>=VExt->Nelem)
	{
		pnpError("Cannot mask since (VExtElm=%d)>=(VExt->Nelem=%d)\n",VExtElm,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot mask since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	for(j=0;j<GS_XYZ;j++)
	{
		if(fabs(VExt->V[VExtElm][j]-Value)<=dvalue)
			V[Elm2set][j]=value2set;
	}
	return EXIT_SUCCESS;
}
int VectorField3D::Copy(VectorField3D* VExt)
{
	if(VExt->Nelem!=Nelem)
	{
		pnpError("Cannot copy since number of elements is not coinside (this->Nelem=%d)!=(VExt->Nelem=%d)\n",Nelem,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot copy since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)V[i][j]=VExt->V[i][j];
	}
	return EXIT_SUCCESS;
}
int VectorField3D::CopyOneElement(VectorField3D* VExt,int VExtElm,int Elm2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot copy since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	if(VExtElm>=VExt->Nelem)
	{
		pnpError("Cannot copy since (VExtElm=%d)>=(VExt->Nelem=%d)\n",VExtElm,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot copy since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(j=0;j<GS_XYZ;j++)V[Elm2set][j]=VExt->V[VExtElm][j];
	
	return EXIT_SUCCESS;
}
int VectorField3D::AddVectorField3D(VectorField3D* VExt)
{
	if(VExt->Nelem!=Nelem)
	{
		pnpError("Cannot add since number of elements is not coinside (this->Nelem=%d)!=(VExt->Nelem=%d)\n",Nelem,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot add since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)V[i][j]+=VExt->V[i][j];
	}
	return EXIT_SUCCESS;
}
int VectorField3D::AddOneElementOfVectorField3D(VectorField3D* VExt,int VExtElm,int Elm2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot add since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	if(VExtElm>=VExt->Nelem)
	{
		pnpError("Cannot add since (VExtElm=%d)>=(VExt->Nelem=%d)\n",VExtElm,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot add since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(j=0;j<GS_XYZ;j++)V[Elm2set][j]+=VExt->V[VExtElm][j];
	
	return EXIT_SUCCESS;
}
int VectorField3D::SubtractVectorField3D(VectorField3D* VExt)
{
	if(VExt->Nelem!=Nelem)
	{
		pnpError("Cannot subtract since number of elements is not coinside (this->Nelem=%d)!=(VExt->Nelem=%d)\n",Nelem,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot subtract since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)V[i][j]-=VExt->V[i][j];
	}
	return EXIT_SUCCESS;
}
int VectorField3D::SubtractOneElementOfVectorField3D(VectorField3D* VExt,int VExtElm,int Elm2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot subtract since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	if(VExtElm>=VExt->Nelem)
	{
		pnpError("Cannot subtract since (VExtElm=%d)>=(VExt->Nelem=%d)\n",VExtElm,VExt->Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GridSize[0]!=GridSize[0])||(VExt->GridSize[1]!=GridSize[1])||(VExt->GridSize[2]!=GridSize[2]))
	{
		pnpError("Cannot subtract since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GridSize[0],VExt->GridSize[1],VExt->GridSize[2]);
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	for(j=0;j<GS_XYZ;j++)V[Elm2set][j]-=VExt->V[VExtElm][j];
	return EXIT_SUCCESS;
}
#ifdef HARLEM_MOD
int VectorField3D::MaskWithHaField3D(HaField3D* VExt,float Value, float dvalue,float value2set)
{
	if((VExt->GetNx()!=GridSize[0])||(VExt->GetNy()!=GridSize[1])||(VExt->GetNz()!=GridSize[2]))
	{
		pnpError("Cannot mask since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GetNx(),VExt->GetNy(),VExt->GetNz());
		return EXIT_FAILURE;
	}
	unsigned int i,j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *hpnt=VExt->GetFieldPtr();
	for(i=0;i<Nelem;i++)
	{
		for(j=0;j<GS_XYZ;j++)
		{
			if(fabs(hpnt[j]-Value)<=dvalue)
				V[i][j]=value2set;
		}
	}
	return EXIT_SUCCESS;
}
int VectorField3D::MaskWithHaField3DUseOneElement(HaField3D* VExt,int Elm2set,float Value, float dvalue,float value2set)
{
	if(Elm2set>=Nelem)
	{
		pnpError("Cannot mask since (Elm2set=%d)>=(Nelem=%d)\n",Elm2set,Nelem);
		return EXIT_FAILURE;
	}
	if((VExt->GetNx()!=GridSize[0])||(VExt->GetNy()!=GridSize[1])||(VExt->GetNz()!=GridSize[2]))
	{
		pnpError("Cannot mask since grids are not coinside (this->GridSize=[%d, %d, %d])!=(VExt->GridSize=[%d, %d, %d])\n",GridSize[0],GridSize[1],GridSize[2],VExt->GetNx(),VExt->GetNy(),VExt->GetNz());
		return EXIT_FAILURE;
	}
	unsigned int j;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *hpnt=VExt->GetFieldPtr();
	for(j=0;j<GS_XYZ;j++)
	{
		if(fabs(hpnt[j]-Value)<=dvalue)
			V[Elm2set][j]=value2set;
	}
	return EXIT_SUCCESS;
}
#endif
void VectorField3D::MinMax(int Elem,float *Min,float *Max,int *MinPnt,int *MaxPnt)
{
	if(Elem>=Nelem)return;
	int j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float _Min=V[Elem][0],_Max=V[Elem][0];
	int _MinPnt=0, _MaxPnt=0;
	
	for(j=1;j<GS_XYZ;j++)
	{
		if(V[Elem][j]>_Max)
		{
			_Max=V[Elem][j];
			_MaxPnt=j;
		}
		if(V[Elem][j]<_Min)
		{
			_Min=V[Elem][j];
			_MinPnt=j;
		}
	}
	*Min=_Min;
	*Max=_Max;
	*MinPnt=_MinPnt;
	*MaxPnt=_MaxPnt;
}
float VectorField3D::RMSD(VectorField3D* vcomp)
{
	if(GridSize[0]!=vcomp->GridSize[0]
		||GridSize[1]!=vcomp->GridSize[1]
		||GridSize[2]!=vcomp->GridSize[2]
		||GridScale!=vcomp->GridScale
		||Nelem!=vcomp->Nelem)
	{
		pnpError("Size of current VectorField3D is not coinside with size of VectorField3D to compare with\n");
		return -1.0;
	}
	double dC;
	double dCMin,dCMax;
	double RMSD;
	int i,j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	
	
	for(j=0;j<Nelem;j++)
	{
		RMSD=0.0;
		i=1+1*GridSize[0]+1*GridSize[0]*GridSize[1];
		dC=fabs(V[j][i]-vcomp->V[j][i]);
		dCMin=dC;
		dCMax=dC;
		int ix,iy,iz;
		for(ix=0;ix<GridSize[0];ix++)
			for(iy=0;iy<GridSize[1];iy++)
				for(iz=0;iz<GridSize[2];iz++)
		{
			i=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			dC=V[j][i]-vcomp->V[j][i];
			RMSD+=dC*dC;
			dC=fabs(dC);
			if(dC<dCMin)dCMin=dC;
			if(dC>dCMax)dCMax=dC;
		}
		RMSD=sqrt(RMSD/GS_XYZ);
		pnpPrint("RMSD[%d] = %14.6e dCMin = %14.6e dCMax = %14.6e\n",j,RMSD,dCMin,dCMax);
	}
	return 0.0;
}
float VectorField3D::RMSDInternal(VectorField3D* vcomp)
{
	if(GridSize[0]!=vcomp->GridSize[0]
			||GridSize[1]!=vcomp->GridSize[1]
			||GridSize[2]!=vcomp->GridSize[2]
			||GridScale!=vcomp->GridScale
			||Nelem!=vcomp->Nelem)
	{
		pnpError("Size of current VectorField3D is not coinside with size of VectorField3D to compare with\n");
		return -1.0;
	}
	double dC;
	double dCMin,dCMax;
	double RMSD;
	int i,j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	
	
	for(j=0;j<Nelem;j++)
	{
		RMSD=0.0;
		i=1+1*GridSize[0]+1*GridSize[0]*GridSize[1];
		dC=fabs(V[j][i]-vcomp->V[j][i]);
		dCMin=dC;
		dCMax=dC;
		int ix,iy,iz;
		for(ix=1;ix<GridSize[0]-1;ix++)
			for(iy=1;iy<GridSize[1]-1;iy++)
				for(iz=1;iz<GridSize[2]-1;iz++)
		{
			i=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			dC=V[j][i]-vcomp->V[j][i];
			RMSD+=dC*dC;
			dC=fabs(dC);
			if(dC<dCMin)dCMin=dC;
			if(dC>dCMax)dCMax=dC;
		}
		RMSD=sqrt(RMSD/((GridSize[0]-2)*(GridSize[1]-2)*(GridSize[2]-2)));
		pnpPrint("RMSD[%d] = %14.6e dCMin = %14.6e dCMax = %14.6e\n",j,RMSD,dCMin,dCMax);
	}
	return 0.0;
}
int VectorField3D::AverageThrLaplas(VectorIntField3D *_Mask,int iter,int MaskNotToDo)
{
  //!Not very right now. Solve FDM laplasian iter times for RF field
	DbgPrint2("RFAverage. iter:%d\n",iter);
	int i,j,k,i1,GridPoint,BnW,iElm;
	int GS_X = GridSize[0];
	int GS_Y = GridSize[1];
	int GS_Z = GridSize[2];
	int GS_XY = GridSize[0]*GridSize[1];
	int GS_XYZ = GridSize[0]*GridSize[1]*GridSize[2];
	for(iElm=0;iElm<Nelem;iElm++)
	{
		for(i1=0;i1<iter;i1++)
		{
			for(BnW=0;BnW<2;BnW++)
				for(GridPoint=BnW;GridPoint<GS_XYZ;GridPoint=GridPoint+2)
			{
				i=GridPoint%GS_X;
				k=GridPoint/GS_XY;
				j=(GridPoint%GS_XY)/GS_X;
				
				if(i>0&&j>0&&k>0&&i<GS_X-1&&j<GS_Y-1&&k<GS_Z-1)
					if(_Mask->V[iElm][GridPoint]!=MaskNotToDo)
				{
					V[iElm][GridPoint]=(V[iElm][GridPoint]*6.0+V[iElm][GridPoint+1]+V[iElm][GridPoint-1]+V[iElm][GridPoint-GS_X]+V[iElm][GridPoint+GS_X]+V[iElm][GridPoint+GS_XY]+V[iElm][GridPoint-GS_XY])/12.0;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
#ifdef HARLEM_MOD
float VectorField3D::RMSDinPore(VectorField3D* vcomp,HaVec_float Rlim,float x0,float y0)
{
	if(GridSize[0]!=vcomp->GridSize[0]
			||GridSize[1]!=vcomp->GridSize[1]
			||GridSize[2]!=vcomp->GridSize[2]
			||GridScale!=vcomp->GridScale
			||Nelem!=vcomp->Nelem)
	{
		pnpError("Size of current VectorField3D is not coinside with size of VectorField3D to compare with\n");
		return -1.0;
	}
	double dC;
	double dCMin,dCMax;
	double RMSD;
	int i,j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	float RSQ,R;
	float RlimSQ;
	
	
	
	for(j=0;j<Nelem;j++)
	{
		RMSD=0.0;
		int count=0;
		bool dCMinInit=false,dCMaxInit=false;
		int ix,iy,iz;
		for(iz=1;iz<GridSize[2]-1;iz++)
			for(ix=1;ix<GridSize[0]-1;ix++)
				for(iy=1;iy<GridSize[1]-1;iy++)
		{
			i=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
			
			RSQ=(ix-x0)*(ix-x0)+(iy-y0)*(iy-y0);
			RlimSQ=Rlim[iz]*Rlim[iz]*GridScale*GridScale;
			if(RSQ<=RlimSQ && (V[j][i]!=0.0 || vcomp->V[j][i]!=0.0) )
			{
				dC=V[j][i]-vcomp->V[j][i];
				RMSD+=dC*dC;
				dC=fabs(dC);
				if(!dCMinInit){dCMin=dC;dCMinInit=true;}
				if(!dCMaxInit){dCMax=dC;dCMaxInit=true;}
				if(dC<dCMin)dCMin=dC;
				if(dC>dCMax)dCMax=dC;
				count++;
			}
		}
		RMSD=sqrt(RMSD/count);
		pnpPrint("RMSD[%d] = %14.6e dCMin = %14.6e dCMax = %14.6e count=%d GS_XYZ=%d\n",j,RMSD,dCMin,dCMax,count,GS_XYZ);
	}
	return 0.0;
}
HaField3D* VectorField3D::GetHaField3D(int Nion,AllocMode AMode)
{
	HaField3D* Field3D=NULL;
	int i;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *HVec=NULL;
	float *VPar=V[Nion];
	if(AMode==INTERNAL_ALLOC)
	{
		Field3D=new HaField3D();
		Field3D->SetDimensions(GridSize[0],GridSize[1],GridSize[2]);
		Field3D->SetCenterAsZero(GridScale);
		HVec=Field3D->GetFieldPtr();
		for(i=0;i<GS_XYZ;i++)HVec[i]=VPar[i];
	}
	else
	{
		 Field3D=new HaField3D(VPar,GridSize[0],GridSize[1],GridSize[2]);
		 Field3D->SetCenterAsZero(GridScale);
	}
	return Field3D;
}
void VectorField3D::SetValuesFromTable(HaVec_int* GrdV,HaVec_double* ValV)
{
	fprintf(stdout,"SetValuesFromTable\n");
	int i,j,count;
	int Size=GrdV->size()/Nelem;
	count=0;
	for(j=0;j<Nelem;j++)
	{
		for(i=0;i<Size;i++)
		{
			V[j][GrdV->GetVal_idx0(count)]=ValV->GetVal_idx0(count);
			count++;
		}
	}
	
}
#endif
////////////////////VectorField3D
////////////////////VectorIntField3D
VectorIntField3D::VectorIntField3D()
{
	amode=INTERNAL_ALLOC;
	InitZero();
}
VectorIntField3D::VectorIntField3D(int *gridsize, float gridscale, int nelem)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	SetVectorField3D(gridsize,gridscale,nelem);
}
VectorIntField3D::VectorIntField3D(int Nx,int Ny,int Nz, float gridscale, int nelem)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	int gridsize[3]={Nx, Ny, Nz};
	SetVectorField3D(gridsize,gridscale,nelem);
}
VectorIntField3D::VectorIntField3D(const char* filename)
{
	amode=INTERNAL_ALLOC;
	InitZero();
	ReadFromFile(filename);
}
VectorIntField3D::VectorIntField3D(int *gridsize, float gridscale, int nelem,int **v)
{
	amode=EXTERNAL_ALLOC;
	InitZero();
	Clear();
	
	GridSize[0]=gridsize[0];
	GridSize[1]=gridsize[1];
	GridSize[2]=gridsize[2];
	GridScale=gridscale;
	Nelem=nelem;
	V=v;
}

VectorIntField3D::~VectorIntField3D()
{
	if(amode==INTERNAL_ALLOC)Clear();
}
int VectorIntField3D::InitZero()
{
	GridSize[0]=0;
	GridSize[1]=0;
	GridSize[2]=0;
	GridScale=0.0f;
	Nelem=0;
	V=NULL;
	return EXIT_SUCCESS;
}
int VectorIntField3D::Clear()
{
	int i;
	
	if(V!=NULL)
	{
		for(i=0;i<Nelem;i++)
		{
			if(V[i]!=NULL)delete [] V[i];
		}
		delete [] V;
	}
	InitZero();
	return EXIT_SUCCESS;
}
int VectorIntField3D::SetVectorField3D(int *gridsize, float gridscale, int nelem)
{
	Clear();
	
	GridSize[0]=gridsize[0];
	GridSize[1]=gridsize[1];
	GridSize[2]=gridsize[2];
	GridScale=gridscale;
	Nelem=nelem;
	int i,j;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	V=new int*[Nelem];
	for(i=0;i<Nelem;i++)
	{
		V[i]=new int[GS_XYZ];
		for(j=0;j<GS_XYZ;j++)
		{
			V[i][j]=0;
		}
	}
	return EXIT_SUCCESS;
}
int VectorIntField3D::ReadFromFile(const char *filename)
{
	if(TypeOfMapFile(filename)==FILE_MAP_GZ)return ReadFromFileGZ(filename);
	else return ReadFromFileBIN(filename);
}
int VectorIntField3D::ReadFromFileBIN(const char *filename)
{
	FILE *file;
	int i,j;
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = fopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	
	
	int size;
	fread(&size,sizeof(int),1,file);
	DbgPrint0("Size %d %d\n",size,sizeof(int));
	char *str;
	str=new char[size+1];
	fread(str,sizeof(char),size,file);
	str[size]='\0';
	
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GS,3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	delete [] str;
	SetVectorField3D(GS,gridscale,nelem);
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	for(i=0;i<Nelem;i++)
	{
		fread(V[i],sizeof(int),GS_XYZ,file);
	}
	fclose(file);
	
	return EXIT_SUCCESS;
	return EXIT_SUCCESS;
}
int VectorIntField3D::ReadFromFileGZ(const char *filename)
{
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GS[3],nelem,Columns;
	float gridscale;
	
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement *header=new TiXmlElement("VectorField3D");
	ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GS,3);
	if(header->GetIntAttribute("Columns",&Columns)==EXIT_FAILURE)Columns=1;
	if(header->GetIntAttribute("Nelem",&nelem)==EXIT_FAILURE)nelem=1;

	if(header->GetFloatAttribute("GridScale",&gridscale)==EXIT_FAILURE)gridscale=1.0f;
	
	delete header;
	
	SetVectorField3D(GS,gridscale,nelem);
	
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int *VPar;
	for(i=0;i<Nelem;i++)
	{
		VPar=V[i];
		if(Columns==1)
		{
		//	ReadMapGZOneColumns(file,V[i],GS_XYZ,coef);
		}
		else
		{
			ReadIndexGZTwoColumns(file,V[i],GS_XYZ);
		}
	}
	gzclose(file);
	return EXIT_SUCCESS;
}
int VectorIntField3D::WriteToFile(const char *filename,int Columns)
{
	if(TypeOfMapFile(filename)==FILE_MAP_GZ)return WriteToFileGZ(filename, Columns);
	else return WriteToFileBIN(filename, Columns);
}
int VectorIntField3D::WriteToFileBIN(const char *filename,int Columns)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GridSize,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	
	FILE *file;
	file = fopen(filename,"wb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	int size=StrHeader.size();
	fwrite(&size,sizeof(int),1,file);
	fwrite(StrHeader.c_str(),sizeof(char),StrHeader.size(),file);
	
	for(i=0;i<Nelem;i++)
	{
		fwrite(V[i],sizeof(int),GS_XYZ,file);
	}
	fclose(file);
	delete header;
	return EXIT_SUCCESS;
}
int VectorIntField3D::WriteToFileGZ(const char *filename,int Columns)
{
	DbgPrint0("VectorField3D::WriteIndexFile(%s)\n",filename);
	//Prepare Values
	unsigned int i;
	float fpoh= 4*M_PI*GridScale;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("VectorField3D");
	header->SetArrOfIntAttribute("GridSize",GridSize,3);
	header->SetFloatAttribute("GridScale",GridScale);
	header->SetIntAttribute("Nelem",Nelem);
	header->SetAttribute("Comment","OneColumns");
	header->SetIntAttribute("Columns",Columns);
	gzFile file;
	file = gzopen(filename,"wb1");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	
	for(i=0;i<Nelem;i++)
	{
		 if(Columns==1)
		 {
			 //WriteMapGZOneColumns(file,V[i],GS_XYZ,coef);
		 }
		 else
		 {
			 WriteIndexGZTwoColumns(file,V[i],GS_XYZ);
		 }
	}
	gzclose(file);
	delete header;
	return EXIT_SUCCESS;
}

#ifdef HARLEM_MOD
HaField3D* VectorIntField3D::GetHaField3D(int Nion)
{
	HaField3D* Field3D=NULL;
	int i;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *HVec=NULL;
	int *VPar=V[Nion];
	Field3D=new HaField3D();
	Field3D->SetDimensions(GridSize[0],GridSize[1],GridSize[2]);
	Field3D->SetCenterAsZero(GridScale);
	HVec=Field3D->GetFieldPtr();
	for(i=0;i<GS_XYZ;i++)HVec[i]=(float)VPar[i];
	return Field3D;
}

#endif
