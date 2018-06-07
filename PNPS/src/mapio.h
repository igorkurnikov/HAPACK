#ifndef MAPIO_H
#define MAPIO_H

#include <zlib.h>

#define FILE_MAP_GZ 4001///< map gz pnp_c style
#define FILE_MAP_MBN 4002///< map binary
#define FILE_MAP_MB2 4003///< map binary 2 coloumn style
#define FILE_MAP_IGB 4004///< map binary igb style for flounder
#define FILE_MAP_DX 4005///< dx

#ifdef HARLEM_MOD
class HaField3D;
//class HaVec_int;
//class HaVec_double;
#include "halinalg.h"
#endif

class TiXmlElement;

int TypeOfMapFile(const char *filename);
//! Write field with additional information in header, TiXmlElement should have only attributes
int WriteMapGZ(const char *filename, TiXmlElement *header,float * fmap, int* gridsize, float coef,char* Comments);
//! Write array from file *.gz of size Size to map, and each element is multiply be coef
int WriteMapGZ(const char *filename, float *fmap, int* gridsize, float coef,int Columns=2,char* Comments="Comments");
//float* ReadMapGZ(const char *filename, int* gridsize,float coef);
int    ReadMapGZ(const char *filename, float *fmap,int* gridsize,float coef);
int    ReadMapGZ(const char *filename,TiXmlElement** header, float **pfmap,int* gridsize,float coef);


int WriteIndexGZinTwoColumns(gzFile file, unsigned int  * nmap,unsigned int GSXYZ);
int  ReadIndexGZfromTwoColumns(gzFile file,unsigned int  * nmap,unsigned int GSXYZ);
int WriteIndexGZTwoColumns(gzFile file, int  * nmap,unsigned int GSXYZ);
int  ReadIndexGZTwoColumns(gzFile file, int  * nmap,unsigned int GSXYZ);

int WriteMapGZOneColumns(gzFile file, float  * nmap,unsigned int N,float coef=1.0);
int ReadMapGZOneColumns(gzFile file, float  * nmap,unsigned  int N,float coef=1.0);
int WriteMapGZTwoColumns(gzFile file, float  * nmap,unsigned int N,float coef=1.0);
int ReadMapGZTwoColumns(gzFile file, float  * nmap,unsigned  int N,float coef=1.0);

/** Convert Array *V to array with periodic boundary boarders. 
 * Do not really aaffect the size of V, size of V should be right from the very beggining.
 * ConvertToPBC will move date to center and copy bouders for PBC calculation
 * *V should be size of [(GS_X+pbcX)*(GS_Y+pbcY)*(GS_Z+pbcZ)]
 * and have [(GS_X)*(GS_Y)*(GS_Z)] written points.
*/
int ConvertToPBC(float *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ);
/** Remove periodic boundaries
 * Convert Arrays from calculative box (GridSize[]) to original box (GridSizeOriginal[]). The difference appears because of periodic boundary condition. mapIn and mapOut can be the same array, after rearragement just mapOut have sence.
*/
int RemovePBC(float *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ);

class VectorIntField3D
{
	public:
		VectorIntField3D();
		VectorIntField3D(int *gridsize, float gridscale, int nelem);
		VectorIntField3D(int Nx,int Ny,int Nz, float gridscale, int nelem);
		VectorIntField3D(int *gridsize, float gridscale, int nelem,int **v);
		VectorIntField3D(const char* filename);
		~VectorIntField3D();
		
		enum AllocMode {INTERNAL_ALLOC=0, EXTERNAL_ALLOC=1};
		AllocMode amode;
		
		int InitZero();
		int Clear();
			
		int GridSize[3];
		float GridScale;
		int Nelem;//!<Number of elements in Vector
		int **V;//!<Vector Field
	
		int SetVectorField3D(int *gridsize, float gridscale, int nelem);
		int ReadFromFile(const char *filename);
			//!will write index table in semi XML way
		int WriteToFile(const char *filename,int Columns=2);
#ifdef HARLEM_MOD
		/*!AMode allocation for HaField3D if INTERNAL_ALLOC will create a new array in HaField3D.
		*/
		HaField3D* GetHaField3D(int Nion);
#endif
		int GetNelem() const { return Nelem; }
		int GetNx() const { return GridSize[0]; } //!< return number of grid points along X axes 
		int GetNy() const { return GridSize[1]; } //!< return number of grid points along Y axes 
		int GetNz() const { return GridSize[2]; } //!< return number of grid points along Z axes 
		int GetValue(int in, int ix, int iy, int iz)
		{
			return V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]];
		}
		void SetValue(int in, int ix, int iy, int iz,int val)
		{
			V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]]=val;
		}
		int GetValueByGrdNmb(int in, int igrd)
		{
			return V[in][igrd];
		}
		void SetValueByGrdNmb(int in, int igrd,int val)
		{
			V[in][igrd]=val;
		}
		void FillValue(int val)
		{
			int i,j,gs=GridSize[0]*GridSize[1]*GridSize[2];
			for(i=0;i<Nelem;i++)
				for(j=0;j<gs;j++)
					V[i][j]=val;
		}
	private:
		int WriteToFileGZ(const char *filename,int Columns=2);
		int WriteToFileBIN(const char *filename,int Columns=2);
		int ReadFromFileGZ(const char *filename);
		int ReadFromFileBIN(const char *filename);
};
template <class TREAL> class HaVectorField3D
{
	public:
		HaVectorField3D()
		{
			amode=INTERNAL_ALLOC;
			InitZero();
		}
		HaVectorField3D(int *gridsize, TREAL gridscale, int nelem)
		{
			amode=INTERNAL_ALLOC;
			InitZero();
			SetVectorField3D(gridsize,gridscale,nelem);
		}
		HaVectorField3D(int Nx,int Ny,int Nz, TREAL gridscale, int nelem)
		{
			amode=INTERNAL_ALLOC;
			InitZero();
			int gridsize[3]={Nx, Ny, Nz};
			SetVectorField3D(gridsize,gridscale,nelem);
		}
		HaVectorField3D(int *gridsize, TREAL gridscale, int nelem,TREAL **v)
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
		~HaVectorField3D()
		{
			if(amode==INTERNAL_ALLOC)Clear();
		}
		
		enum AllocMode {INTERNAL_ALLOC=0, EXTERNAL_ALLOC=1};
		AllocMode amode;
		
		int InitZero()
		{
			GridSize[0]=0;
			GridSize[1]=0;
			GridSize[2]=0;
			GridScale=0.0f;
			Nelem=0;
			V=NULL;
			return EXIT_SUCCESS;
		}
		int Clear()
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
		
		int GridSize[3];
		TREAL GridScale;
		int Nelem;//!<Number of elements in Vector
		TREAL **V;//!<Vector Field
		
		int SetVectorField3D(int *gridsize, TREAL gridscale, int nelem)
		{
			Clear();
			
			GridSize[0]=gridsize[0];
			GridSize[1]=gridsize[1];
			GridSize[2]=gridsize[2];
			GridScale=gridscale;
			Nelem=nelem;
			int i,j;
			int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
			V=new TREAL*[Nelem];
			for(i=0;i<Nelem;i++)
			{
				V[i]=new TREAL[GS_XYZ];
				for(j=0;j<GS_XYZ;j++)
				{
					V[i][j]=0.0f;
				}
			}
			return EXIT_SUCCESS;
		}
		int GetNelem() const { return Nelem; }
		int GetNx() const { return GridSize[0]; } //!< return number of grid points along X axes 
		int GetNy() const { return GridSize[1]; } //!< return number of grid points along Y axes 
		int GetNz() const { return GridSize[2]; } //!< return number of grid points along Z axes 
		TREAL GetValue(int in, int ix, int iy, int iz)
		{
			return V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]];
		}
		void SetValue(int in, int ix, int iy, int iz,TREAL val)
		{
			V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]]=val;
		}
		TREAL GetValueByGrdNmb(int in, int igrd)
		{
			return V[in][igrd];
		}
		void SetValueByGrdNmb(int in, int igrd,TREAL val)
		{
			V[in][igrd]=val;
		}
		void FillValue(TREAL val)
		{
			int i,j,gs=GridSize[0]*GridSize[1]*GridSize[2];
			for(i=0;i<Nelem;i++)
				for(j=0;j<gs;j++)
					V[i][j]=val;
		}
		void FillElementWithValue(int i,TREAL val)
		{
			int j,gs=GridSize[0]*GridSize[1]*GridSize[2];
			for(j=0;j<gs;j++)
				V[i][j]=val;
		}
		int MultiplyBy(TREAL coef)
		{
			unsigned int i,j;
			unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
			for(i=0;i<Nelem;i++)
			{
				for(j=0;j<GS_XYZ;j++)V[i][j]*=coef;
			}
			return EXIT_SUCCESS;
		}
		void MinMax(int Elem,TREAL *Min,TREAL *Max,int *MinPnt,int *MaxPnt)
		{
			if(Elem>=Nelem)return;
			int j;
			int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
			TREAL _Min=V[Elem][0],_Max=V[Elem][0];
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
};
typedef HaVectorField3D<double> HaPyDoubleVectorField3D;

/*class HaDoubleVectorField3D:public HaDoubleVectorField3Dparent
{
	public:
		HaDoubleVectorField3D()
		{
			amode=INTERNAL_ALLOC;
			InitZero();
		}
		HaDoubleVectorField3D(int *gridsize, double gridscale, int nelem)
		{
			amode=INTERNAL_ALLOC;
			InitZero();
			SetVectorField3D(gridsize,gridscale,nelem);
		}
		HaDoubleVectorField3D(int Nx,int Ny,int Nz, double gridscale, int nelem)
		{
			amode=INTERNAL_ALLOC;
			InitZero();
			int gridsize[3]={Nx, Ny, Nz};
			SetVectorField3D(gridsize,gridscale,nelem);
		}
		HaDoubleVectorField3D(int *gridsize, double gridscale, int nelem,double **v)
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
		~HaDoubleVectorField3D()
		{
			if(amode==INTERNAL_ALLOC)Clear();
		}
};*/
class VectorField3D
{
	public:
		VectorField3D();
		VectorField3D(int *gridsize, float gridscale, int nelem);
		VectorField3D(int Nx,int Ny,int Nz, float gridscale, int nelem);
		VectorField3D(int *gridsize, float gridscale, int nelem,float **v);
		VectorField3D(const char* filename,float coef=1.0f);
		~VectorField3D();
		
		enum AllocMode {INTERNAL_ALLOC=0, EXTERNAL_ALLOC=1};
		AllocMode amode;
		
		int InitZero();
		int Clear();
		
		int GridSize[3];
		float GridScale;
		int Nelem;//!<Number of elements in Vector
		float **V;//!<Vector Field
		
		int SetVectorField3D(int *gridsize, float gridscale, int nelem);
		int ReadFromFile(const char *filename,float coef=1.0f);
		//!will write index table in semi XML way
		int WriteToFile(const char *filename,float coef=1.0f,int Columns=1);
		
		float RMSD(VectorField3D* vcomp);
		float RMSDInternal(VectorField3D* vcomp);
#		ifdef HARLEM_MOD
		float RMSDinPore(VectorField3D* vcomp,HaVec_float Rlim,float x0,float y0);
		/*!AMode allocation for HaField3D if INTERNAL_ALLOC will create a new array in HaField3D.
		*/
		HaField3D* GetHaField3D(int Nion,AllocMode AMode=INTERNAL_ALLOC);
		void SetValuesFromTable(HaVec_int* GrdV,HaVec_double* ValV);
#		endif
		int GetNelem() const { return Nelem; }
		int GetNx() const { return GridSize[0]; } //!< return number of grid points along X axes 
		int GetNy() const { return GridSize[1]; } //!< return number of grid points along Y axes 
		int GetNz() const { return GridSize[2]; } //!< return number of grid points along Z axes 
		float GetValue(int in, int ix, int iy, int iz)
		{
			return V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]];
		}
		//!trilinear interpolation or something like this, (0.0,0.0,0.0) is located at (0,0,0) and fx,fy,fz is in grids
		float GetInterpolatedValueGrid(int in, float fx, float fy, float fz);
		void SetValue(int in, int ix, int iy, int iz,float val)
		{
			V[in][ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1]]=val;
		}
		float GetValueByGrdNmb(int in, int igrd)
		{
			return V[in][igrd];
		}
		void SetValueByGrdNmb(int in, int igrd,float val)
		{
			V[in][igrd]=val;
		}
		void FillValue(float val)
		{
			int i,j,gs=GridSize[0]*GridSize[1]*GridSize[2];
			for(i=0;i<Nelem;i++)
				for(j=0;j<gs;j++)
					V[i][j]=val;
		}
		int MultiplyBy(float coef);
		int MultiplyOneElementBy(float coef,int Elm2set);
		int MaskWithVectorField3D(VectorField3D* VExt,float Value, float dvalue,float value2set);
		int MaskWithVectorField3DUseOneElement(VectorField3D* VExt,int VExtElm,int Elm2set,float Value, float dvalue,float value2set);
		int Copy(VectorField3D* VExt);
		int CopyOneElement(VectorField3D* VExt,int VExtElm,int Elm2set);
		int AddVectorField3D(VectorField3D* VExt);
		int AddOneElementOfVectorField3D(VectorField3D* VExt,int VExtElm,int Elm2set);
		int SubtractVectorField3D(VectorField3D* VExt);
		int SubtractOneElementOfVectorField3D(VectorField3D* VExt,int VExtElm,int Elm2set);
#		ifdef HARLEM_MOD
		int MaskWithHaField3D(HaField3D* VExt,float Value, float dvalue,float value2set);
		int MaskWithHaField3DUseOneElement(HaField3D* VExt,int Elm2set,float Value, float dvalue,float value2set);
#		endif
		void MinMax(int Elem,float *Min,float *Max,int *MinPnt,int *MaxPnt);

		int ReadFromFileAddPBC(const char *filename, float coef, bool pbcX, bool pbcY, bool pbcZ);
		//!Read map from file add PBC if needed and split among processes
		int ReadFromFileAddPBCandSplitSDMPI4FDaZ(const char *filename, float coef, bool *pbc,int LocaliZ0,int LocaliZ1);
		int SplitExtVF3DSDMPI4FDaZ(VectorField3D *VExt,int LocaliZ0,int LocaliZ1);
		int WriteToFileRemovePBC(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ);
		int WriteToFileCombineSDMPIDistrRemovePBC(const char *filename,float coef,int Columns, int *GridSizeGlobal,bool *pbc,int LocaliZ0,int LocaliZ1);
		
		inline float CalcLinInterFloat(float V0,float V1,float x)
		{return V0+x*(V1-V0);}
		//!interpolate arrays from external array, physical size of VExt should be bigger
		int InterpolateFromExt(VectorField3D* VExt);
		int InterpolateBoarderFromExt(VectorField3D* VExt);
		int InterpolateBoarderFromExtWithMultipl(VectorField3D* VExt,float Multipl);
		int InterpolateInternalBoarderFromExt(VectorField3D* VExt,int *InternalBoarderMinGridLev0,int *InternalBoarderMaxGridLev0);
		int InterpolateInternalBoarderFromExtWithMultipl(VectorField3D* VExt,int *InternalBoarderMinGridLev0,int *InternalBoarderMaxGridLev0,float Multipl);
	public:
		inline float ConvFloatToIntUnitsX(float Xext)
		{
			return Xext*GridScale+(float)int(GridSize[0]/2);
		}
		inline float ConvFloatToIntUnitsY(float Yext)
		{
			return Yext*GridScale+(float)int(GridSize[1]/2);
		}
		inline float ConvFloatToIntUnitsZ(float Zext)
		{
			return Zext*GridScale+(float)int(GridSize[2]/2);
		}
		inline int ConvIntXYZToGrdPnt(float Xint,float Yint,float Zint)
		{
			return int((Xint+0.5))+int((Yint+0.5))*GridSize[0]+int((Zint+0.5))*GridSize[0]*GridSize[1];
		}
		inline int ConvExtXYZToGrdPnt(float Xext,float Yext,float Zext)
		{
			return ConvIntXYZToGrdPnt(ConvFloatToIntUnitsX(Xext), ConvFloatToIntUnitsY(Yext), ConvFloatToIntUnitsZ(Zext));
		}
		inline float ConvIntToExtUnitsX(float Xint)
		{
			return (Xint-(float)int(GridSize[0]/2))/GridScale;
		}
		inline float ConvIntToExtUnitsY(float Yint)
		{
			return (Yint-(float)int(GridSize[1]/2))/GridScale;
		}
		inline float ConvIntToExtUnitsZ(float Zint)
		{
			return (Zint-(float)int(GridSize[2]/2))/GridScale;
		}
	private:
		int WriteToFileGZ(const char *filename,float coef=1.0f,int Columns=1);
		int WriteToFileBIN(const char *filename,float coef=1.0f,int Columns=1);
		int ReadFromFileGZ(const char *filename,float coef=1.0f);
		int ReadFromFileBIN(const char *filename,float coef=1.0f);
		int ReadFromFileDX(const char *filename,float coef=1.0f);
		
		int ReadFromFileAddPBC_GZ(const char *filename, float coef, bool pbcX, bool pbcY, bool pbcZ);
		int ReadFromFileAddPBC_BIN(const char *filename, float coef, bool pbcX, bool pbcY, bool pbcZ);
		int WriteToFileRemovePBC_GZ(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ);
		int WriteToFileRemovePBC_BIN(const char *filename,float coef,int Columns, bool pbcX, bool pbcY, bool pbcZ);
	public:
		int AverageThrLaplas(VectorIntField3D *_Mask ,int iter,int MaskNotToDo);
};
#endif
