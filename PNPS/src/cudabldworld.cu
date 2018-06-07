#ifndef _CUDA_BLD_WORLD_CU_
#define _CUDA_BLD_WORLD_CU_

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

#include "pnpstructs.h"
#include "pnpdebug.h"


#define FMUL __dmul_rn
#define FADD __dadd_rn
#define FMAF __dmaf_rn

#define BIGDISTANSE 10000

int checkCUDAError(const char* msg);

extern "C" int DoBldDielMapsOnCuda()
{
	return 1;
}

//BuildWorldNI *Builder,ContWorld* world,
extern "C" int BuildAtomsDielPreMapsOnCUDA0(GOAtomsStruct* atms,float *Displ)
{
	int iValue=atms->iDiel;
	int iBulkValue=atms->iDielBulk;
	float Rsmoth=atms->Rsmoth;
	int NAtoms=atms->Natoms;
	float *r[3]={atms->r[0],atms->r[1],atms->r[2]};
	float *R=atms->R;
	float *Surf=atms->Surf;
	int *Field=atms->iVtmp;
	DbgPrint2("GOAtoms::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
// 	DbgPrint2("\t Rion=%f[grids] Rsmoth=%f[grids] \nDispl=[%f,%f,%f][grids,grids,grids]\n", Rion, Rsmoth, Displ[0], Displ[1], Displ[2]);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][0],r[1][0],r[2][0],R[0],q[0],NAtoms);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][NAtoms-1],r[1][NAtoms-1],r[2][NAtoms-1],R[NAtoms-1],q[NAtoms-1],NAtoms);
	int i,j,k,gridpoint,rint[3];
	float RSQ,RsmSQ,RtmpSQ,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	int GridSize[3]={atms->GS_X,atms->GS_Y,atms->GS_Z};
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	
	int start[3];
	int end[3];
	
	if(Rsmoth==0.0f)
	{
		DbgPrint0("Rsmoth==0.0f");
		for(i=0;i<NAtoms;i++)
		{
			iR=(int)(R[i]+0.5);
			RSQ=R[i]*R[i];
			for(k=0;k<3;k++){
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iR;
				end[k]=rint[k]+iR;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
				}
			}
		}
	}
	else
	{
		for(i=0;i<NAtoms;i++)
		{
			Rsm=R[i]+Rsmoth;
			iR=(int)(R[i]+0.5);
			iRsm=(int)(Rsm+1.0);
			RleftSQ=Rsm-0.707106781f;
			RrightSQ=Rsm+0.707106781f;
			RleftSQ*=RleftSQ;
			RrightSQ*=RrightSQ;
			RSQ=R[i]*R[i];
			RsmSQ=Rsm*Rsm;
			
			for(k=0;k<3;k++)
			{
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iRsm;
				end[k]=rint[k]+iRsm;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
			}
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
				else//RSQ<RtmpSQ
				{
					Rtmp=sqrt(RtmpSQ);
					if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
					{
						if(Field[gridpoint]==iBulkValue)
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-R[i];
						}
						else if(Field[gridpoint]<0&&Rtmp-R[i]<vtmp[3])
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-R[i];
						}
					}
					if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
					{
						if(vtmp[0]>-100.0f)//e.i. is intersection
						{
							if(Rtmp-R[i]<vtmp[3])
							{
								ftmp=Rsm/Rtmp;
								vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
								vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
								vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
								vtmp[3]=Rtmp-R[i];
							}
						}
						else if(Field[gridpoint]==iBulkValue)
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-R[i];
						}
					}
				}
			}
		}
		/*if(MakePreRoll)
		{
			
			DbgPrint0("MakePreRoll\n");
			int lim=(int)(Rsmoth+0.5);
			int itmp=GSXYZ*4;
			int iSolAcFlag=1000000000;
			RsmSQ=Rsmoth*Rsmoth;
			for(i=0;i<itmp;i=i+4)
			{
				if(Surf[i]>-100.0)
				{
					vtmp=Surf+i;
					for(k=0;k<3;k++)
					{
						rint[k]=(int)(vtmp[k]+0.5);
						start[k]=rint[k]-lim;
						end[k]=rint[k]+lim;
						if(start[k]<0)start[k]=0;
						if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
					}
					for(ix=start[0];ix<=end[0];ix++)
						for(iy=start[1];iy<=end[1];iy++)
							for(iz=start[2];iz<=end[2];iz++)
					{
						gridpoint=ix+iy*GSX+iz*GSXY;
						RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
						if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
						{
							Field[gridpoint]=iSolAcFlag;
						}
					}
				}
			}
			for(i=0;i<GSXYZ;i++)
			{
				if(Field[i]<0)
				{
					Field[i]=-Field[i];
				}
				else if(Field[i]==iSolAcFlag)
				{
					Field[i]=-iValue;
				}
			}
		}*/
	}
	return EXIT_SUCCESS;
}
__global__ void KBuildAtomsDielPreMaps(float *x,float *y,float *z,float *R, float *Surf,int *Field,int GS_X,int GS_Y,int GS_Z,float Rsmoth,float DisplX,float DisplY,float DisplZ,int iValue, int iBulkValue,int Natoms)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	
	
	
	if(i<Natoms)
	{
		float *r[3]={x,y,z};
		float Displ[3]={DisplX, DisplY, DisplZ};
		int GridSize[3]={GS_X, GS_Y, GS_Z};
		int GS_XY=GS_X*GS_Y;
		float Rsm=R[i]+Rsmoth;
		int iRsm=(int)(Rsm+1.0);
		float RleftSQ=Rsm-0.707106781f;
		float RrightSQ=Rsm+0.707106781f;
		RleftSQ*=RleftSQ;
		RrightSQ*=RrightSQ;
		float RSQ=R[i]*R[i];
		float RsmSQ=Rsm*Rsm;
		
		int k;
		int start[3];
		int end[3];
		float rf[3];
		int rint[3];
		int ix,iy,iz;
		int gridpoint;
		
		float Rtmp;
		float ftmp;
		float RtmpSQ;
		float *vtmp;
		
		
		for(k=0;k<3;k++)
		{
			rf[k]=r[k][i]+Displ[k];
			rint[k]=(int)(rf[k]+0.5);
			start[k]=rint[k]-iRsm;
			end[k]=rint[k]+iRsm;
			if(start[k]<0)start[k]=0;
			if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
		}
		
		for(ix=start[0];ix<=end[0];ix++)
			for(iy=start[1];iy<=end[1];iy++)
				for(iz=start[2];iz<=end[2];iz++)
		{
			gridpoint=ix+iy*GS_X+iz*GS_XY;
			RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
			vtmp=Surf+4*gridpoint;
			if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
			{
				Field[gridpoint]=iValue;
				vtmp[0]=-BIGDISTANSE;
				vtmp[3]=0.0;
			}
			else//RSQ<RtmpSQ
			{
				Rtmp=sqrt(RtmpSQ);
				if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
				{
					if(Field[gridpoint]==iBulkValue)
					{
						vtmp[0]=-BIGDISTANSE;
						Field[gridpoint]=-iValue;
						vtmp[3]=Rtmp-R[i];
					}
					else if(Field[gridpoint]<0&&Rtmp-R[i]<vtmp[3])
					{
						vtmp[0]=-BIGDISTANSE;
						Field[gridpoint]=-iValue;
						vtmp[3]=Rtmp-R[i];
					}
				}
				if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
				{
					if(vtmp[0]>-100.0f)//e.i. is intersection
					{
						if(Rtmp-R[i]<vtmp[3])
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-R[i];
						}
					}
					else if(Field[gridpoint]==iBulkValue)
					{
						ftmp=Rsm/Rtmp;
						vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
						vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
						vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
						vtmp[3]=Rtmp-R[i];
					}
				}
			}
		}
	}
}
extern "C" int BuildAtomsDielPreMapsOnCUDA(GOAtomsStruct* atms,float *Displ)
{
	int iValue=atms->iDiel;
	int iBulkValue=atms->iDielBulk;
	float Rsmoth=atms->Rsmoth;
	int Natoms=atms->Natoms;
	float *r[3]={atms->r[0],atms->r[1],atms->r[2]};
	float *R=atms->R;
	float *Surf=atms->Surf;
	int *Field=atms->iVtmp;
	DbgPrint2("GOAtoms::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
// 	DbgPrint2("\t Rion=%f[grids] Rsmoth=%f[grids] \nDispl=[%f,%f,%f][grids,grids,grids]\n", Rion, Rsmoth, Displ[0], Displ[1], Displ[2]);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][0],r[1][0],r[2][0],R[0],q[0],NAtoms);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][NAtoms-1],r[1][NAtoms-1],r[2][NAtoms-1],R[NAtoms-1],q[NAtoms-1],NAtoms);
	int i,j,k,gridpoint,rint[3];
	float RSQ,RsmSQ,RtmpSQ,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	int GridSize[3]={atms->GS_X,atms->GS_Y,atms->GS_Z};
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	int start[3];
	int end[3];
	
	int DoGPU=1;
	if(DoGPU)
	{
		int BlockSize=256;
		int GridOfBlocksSize=Natoms/BlockSize;
		if(Natoms%BlockSize!=0)
			GridOfBlocksSize++;
		
		float *cuda_r[3];
		float *cuda_R;
		float *cuda_Surf;
		int *cuda_Field;
		
		cudaMalloc((void**)&cuda_r[0], Natoms*sizeof(float));
		cudaMalloc((void**)&cuda_r[1], Natoms*sizeof(float));
		cudaMalloc((void**)&cuda_r[2], Natoms*sizeof(float));
		
		cudaMalloc((void**)&cuda_R, GSXYZ*sizeof(float));
		
		cudaMalloc((void**)&cuda_Field, GSXYZ*sizeof(int));
		cudaMalloc((void**)&cuda_Surf, GSXYZ*4*sizeof(float));
			
		cudaMemcpy(cuda_r[0], r[0], Natoms*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(cuda_r[1], r[1], Natoms*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(cuda_r[2], r[2], Natoms*sizeof(float), cudaMemcpyHostToDevice);
		
		cudaMemcpy(cuda_R, R, Natoms*sizeof(float), cudaMemcpyHostToDevice);
		
		cudaMemcpy(cuda_Field, Field, GSXYZ*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(cuda_Surf, Surf, GSXYZ*4*sizeof(float), cudaMemcpyHostToDevice);
		
		KBuildAtomsDielPreMaps<<<GridOfBlocksSize, BlockSize>>>(cuda_r[0],cuda_r[1],cuda_r[2],cuda_R,cuda_Surf,cuda_Field, atms->GS_X,atms->GS_Y,atms->GS_Z,Rsmoth, Displ[0], Displ[1], Displ[2],iValue, iBulkValue, Natoms);
		cudaThreadSynchronize();
		
		cudaMemcpy(Field, cuda_Field, GSXYZ*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(Surf, cuda_Surf, GSXYZ*4*sizeof(float), cudaMemcpyDeviceToHost);
		
		cudaFree(cuda_Surf);
		cudaFree(cuda_Field);
		cudaFree(cuda_R);
		cudaFree(cuda_r[2]);
		cudaFree(cuda_r[1]);
		cudaFree(cuda_r[0]);
		
		printf("GridSize [%d %d %d]\n",GridSize[0],GridSize[1],GridSize[2]);
		printf("\tGridSize [%d %d %d]\n",Field[0],Field[1],Field[2]);
		printf("r [%f %f %f]\n",r[0][0],r[1][0],r[2][0]);
		printf("\tr [%f %f %f]\n",Surf[0],Surf[1],Surf[2]);
	}
	else
	{
		for(i=0;i<Natoms;i++)
		{
			Rsm=R[i]+Rsmoth;
			iRsm=(int)(Rsm+1.0);
			RleftSQ=Rsm-0.707106781f;
			RrightSQ=Rsm+0.707106781f;
			RleftSQ*=RleftSQ;
			RrightSQ*=RrightSQ;
			RSQ=R[i]*R[i];
			RsmSQ=Rsm*Rsm;
			
			for(k=0;k<3;k++)
			{
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iRsm;
				end[k]=rint[k]+iRsm;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
			}
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
				else//RSQ<RtmpSQ
				{
					Rtmp=sqrt(RtmpSQ);
					if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
					{
						if(Field[gridpoint]==iBulkValue)
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-R[i];
						}
						else if(Field[gridpoint]<0&&Rtmp-R[i]<vtmp[3])
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-R[i];
						}
					}
					if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
					{
						if(vtmp[0]>-100.0f)//e.i. is intersection
						{
							if(Rtmp-R[i]<vtmp[3])
							{
								ftmp=Rsm/Rtmp;
								vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
								vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
								vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
								vtmp[3]=Rtmp-R[i];
							}
						}
						else if(Field[gridpoint]==iBulkValue)
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-R[i];
						}
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
__global__ void KFinalazeSEVOnCUDA(float3 *surf_points,int *Field, float Rsmooth, int iBulkValue,int4 GS,int Nsurf_points)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	int ix,iy,iz,gridpoint;
	float RsmSQ,RtmpSQ;
	
	int3 rint;
	int3 start;
	int3 end;
	
	if(i<Nsurf_points)
	{
		float RsmSQ=Rsmooth*Rsmooth;
		int lim=(int)(Rsmooth+0.5);
		
		rint.x=(int)(surf_points[i].x+0.5);
		rint.y=(int)(surf_points[i].y+0.5);
		rint.z=(int)(surf_points[i].z+0.5);
		
		start.x=rint.x-lim;
		start.y=rint.y-lim;
		start.z=rint.z-lim;
		
		end.x=rint.x+lim;
		end.y=rint.y+lim;
		end.z=rint.z+lim;
		
		if(start.x<0)start.x=0;
		if(start.y<0)start.y=0;
		if(start.z<0)start.z=0;
		
		if(end.x>GS.x-1)end.x=GS.x-1;
		if(end.y>GS.y-1)end.x=GS.y-1;
		if(end.z>GS.z-1)end.x=GS.z-1;
		
		for(ix=start.x;ix<=end.x;ix++)
			for(iy=start.y;iy<=end.y;iy++)
				for(iz=start.z;iz<=end.z;iz++)
		{
			gridpoint=ix+iy*GS.x+iz*GS.w;
			RtmpSQ=(surf_points[i].x-ix)*(surf_points[i].x-ix)+(surf_points[i].y-iy)*(surf_points[i].y-iy)+(surf_points[i].z-iz)*(surf_points[i].z-iz);
			if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
			{
				Field[gridpoint]=iBulkValue;
			}
		}
	}
}
extern "C" int FinalazeSEVOnCUDA2(int *GridSize,int *Field,int iBulkValue,float Rsmooth,float3 *surf_points, int Nsurf_points)
{
	int i;
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int itmp=GS_XYZ*4;
	
	int BlockSize=512;
	int GridOfBlocksSize=Nsurf_points/BlockSize;
	if(GridOfBlocksSize>65535)
		printf("Error: GridOfBlocksSize>65535\n");
	if(GS_XYZ%BlockSize!=0)
		GridOfBlocksSize++;
		
	printf("BlockSize=%d GridOfBlocksSize=%d\n",BlockSize,GridOfBlocksSize);
		
	float3 *cuda_surf_points;
	int *cuda_Field;
			
	cudaMalloc((void**)&cuda_Field, GS_XYZ*sizeof(int));
	cudaMalloc((void**)&cuda_surf_points, Nsurf_points*sizeof(float3));
				
	cudaMemcpy(cuda_Field, Field, GS_XYZ*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_surf_points, surf_points, Nsurf_points*sizeof(float3), cudaMemcpyHostToDevice);
		
	int4 GS={GridSize[0],GridSize[1],GridSize[2],GridSize[0]*GridSize[1]};
		
	KFinalazeSEVOnCUDA<<<GridOfBlocksSize, BlockSize>>>(cuda_surf_points,cuda_Field,Rsmooth, iBulkValue,GS,Nsurf_points);
	cudaThreadSynchronize();
			
	cudaMemcpy(Field, cuda_Field, GS_XYZ*sizeof(int), cudaMemcpyDeviceToHost);
			
	cudaFree(cuda_surf_points);
	cudaFree(cuda_Field);
	checkCUDAError("FinalazeSEVOnCUDA");
}
extern "C" int FinalazeSEVOnCUDA(int *GridSize,int *Field,int iBulkValue,float Rsmooth,float *Surf)
{
	DbgPrint2("BuildSES iBulkValue=%d\n",iBulkValue);
	int DoGPU=1;
	if(DoGPU)
	{
		int i;
		int Nsurf_points=0;
		int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
		int itmp=GS_XYZ*4;
		
		for(i=0;i<itmp;i=i+4)
			if(Surf[i]>-100.0)
				Nsurf_points++;
		
		float3* surf_points=new float3[Nsurf_points];
		printf("Nsurf_points=%d\n",Nsurf_points);
		
		int count=0;
		for(i=0;i<itmp;i=i+4)
			if(Surf[i]>-100.0)
		{
			surf_points[count].x=Surf[i];
			surf_points[count].y=Surf[i+1];
			surf_points[count].z=Surf[i+2];
			count++;
		}
		
		int BlockSize=512;
		int GridOfBlocksSize=Nsurf_points/BlockSize;
		if(GridOfBlocksSize>65535)
			printf("Error: GridOfBlocksSize>65535\n");
		if(GS_XYZ%BlockSize!=0)
			GridOfBlocksSize++;
		
		printf("BlockSize=%d GridOfBlocksSize=%d\n",BlockSize,GridOfBlocksSize);
		
		float3 *cuda_surf_points;
		int *cuda_Field;
			
		cudaMalloc((void**)&cuda_Field, GS_XYZ*sizeof(int));
		cudaMalloc((void**)&cuda_surf_points, Nsurf_points*sizeof(float3));
				
		cudaMemcpy(cuda_Field, Field, GS_XYZ*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(cuda_surf_points, surf_points, Nsurf_points*sizeof(float3), cudaMemcpyHostToDevice);
		
		int4 GS={GridSize[0],GridSize[1],GridSize[2],GridSize[0]*GridSize[1]};
		
		KFinalazeSEVOnCUDA<<<GridOfBlocksSize, BlockSize>>>(cuda_surf_points,cuda_Field,Rsmooth, iBulkValue,GS,Nsurf_points);
		cudaThreadSynchronize();
			
		cudaMemcpy(Field, cuda_Field, GS_XYZ*sizeof(int), cudaMemcpyDeviceToHost);
			
		cudaFree(cuda_surf_points);
		cudaFree(cuda_Field);
		delete [] surf_points;
		checkCUDAError("FinalazeSEVOnCUDA");
		
		for(i=0;i<GS_XYZ;i++)
		{
			if(Field[i]<0)
			{
				Field[i]=-Field[i];
			}
		}
	}
	else
	{
		int i,j,k,ix,iy,iz,gridpoint,rint[3];
		float RsmSQ,RtmpSQ;
		int GS_X=GridSize[0];
		int GS_XY=GridSize[0]*GridSize[1];
		int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
		int start[3];
		int end[3];
		int lim=(int)(Rsmooth+0.5);
		int itmp=GS_XYZ*4;
		float *vtmp;
		
		RsmSQ=Rsmooth*Rsmooth;
		for(i=0;i<itmp;i=i+4)
		{
			if(Surf[i]>-100.0)
			{
				//DbgPrint2("BuildSES %d\n",i/4);
				vtmp=Surf+i;
				for(k=0;k<3;k++)
				{
					rint[k]=(int)(vtmp[k]+0.5);
					start[k]=rint[k]-lim;
					end[k]=rint[k]+lim;
					if(start[k]<0)start[k]=0;
					if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
				}
				for(ix=start[0];ix<=end[0];ix++)
					for(iy=start[1];iy<=end[1];iy++)
						for(iz=start[2];iz<=end[2];iz++)
				{
					gridpoint=ix+iy*GS_X+iz*GS_XY;
					RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
					if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
					{
						Field[gridpoint]=iBulkValue;
					}
				}
			}
		}
		for(i=0;i<GS_XYZ;i++)
		{
			if(Field[i]<0)
			{
				Field[i]=-Field[i];
			}
		}
	}
	//RemovingCavitiesOnDielectricMap
// 	if(RemovingCavitiesOnDielectricMap)
// 	{
// 		DbgPrint2("RemovingCavitiesOnDielectricMap\n");
// 		RemovingCavitiesAtValues(GS_X,GS_Y,GS_Z,Field, RemCavOnDielWhere2Look ,RemCavOnDielFillWith);
// 	}
	return EXIT_SUCCESS;
}
__global__ void KGOAtomsCoulBCXYOnCUDA(float4 *rq,int Nq, int4 GS,float eps,float *Potential)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	int j=blockIdx.y;
	int iq;
	int k,GrPnt;
	float r1;
	
	
	if(i<GS.x && j<GS.y)
	{
		k=0;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
		k=GS.z-1;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
	}
}
__global__ void KGOAtomsCoulBCYZOnCUDA(float4 *rq,int Nq, int4 GS,float eps,float *Potential, char3 cBldBCatPlane)
{
	int j=blockIdx.x*blockDim.x+threadIdx.x;
	int k=blockIdx.y+cBldBCatPlane.z;
	int iq;
	int i,GrPnt;
	float r1;
	
	if(j<GS.y && k<GS.z)
	{
		i=0;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
		i=GS.x-1;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
	}
}
__global__ void KGOAtomsCoulBCXZOnCUDA(float4 *rq,int Nq, int4 GS,float eps,float *Potential, char3 cBldBCatPlane)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x+cBldBCatPlane.x;
	int j;
	int k=blockIdx.y+cBldBCatPlane.z;
	int iq;
	int GrPnt;
	float r1;
	if(i<GS.x && k<GS.z)
	{
		j=0;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
		j=GS.y-1;
		GrPnt=i+j*GS.x+k*GS.w;
					
		for(iq=0;iq < Nq;iq++)
		{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
	}
}
extern "C" int GOAtoms_SetCoulombicBC(int *GridSize,int *BldBCatPlane,float4* rq,int Nq,float eps,float *Potential)
{
	int DoGPU=1;
	if(DoGPU)
	{
		int4 GS;
		GS.x=GridSize[0];
		GS.y=GridSize[1];
		GS.z=GridSize[2];
		GS.w=GridSize[0]*GridSize[1];
		
		float4 *cuda_rq;
		
		cudaMalloc((void**)&cuda_rq, Nq*sizeof(float4));
		
		
		cudaMemcpy(cuda_rq, rq, Nq*sizeof(float4), cudaMemcpyHostToDevice);
		
		dim3 dimBlock(512,1,1);
		dim3 dimGrid(1,1,1);
		char3 cBldBCatPlane;
		
		cBldBCatPlane.x=BldBCatPlane[0];
		cBldBCatPlane.y=BldBCatPlane[1];
		cBldBCatPlane.z=BldBCatPlane[2];
		
		if(BldBCatPlane[2])
		{
			float *cuda_pot;
			cudaMalloc((void**)&cuda_pot, GS.w*GS.z*sizeof(float));
			cudaMemcpy(cuda_pot, Potential, GS.w*GS.z*sizeof(float), cudaMemcpyHostToDevice);
			
			if(GS.x<=512)
			{
				dimBlock.x=GS.x;
				dimGrid.x=1;
			}
			else
			{
				dimBlock.x=512;
				dimGrid.x=GS.x/512;
				if(GS.x%512!=0)dimGrid.x++;
			}
			dimGrid.y=GS.y;
			
			printf("dimBlock [%d,%d,%d]\n",dimBlock.x,dimBlock.y,dimBlock.z);
			printf("dimGrid [%d,%d,%d]\n",dimGrid.x,dimGrid.y,dimGrid.z);
			KGOAtomsCoulBCXYOnCUDA<<<dimGrid, dimBlock>>>(cuda_rq,Nq, GS,eps, cuda_pot);
			cudaThreadSynchronize();
			
			cudaMemcpy(Potential, cuda_pot, GS.w*GS.z*sizeof(float), cudaMemcpyDeviceToHost);
			
			cudaFree(cuda_pot);
		}
		if(BldBCatPlane[0])
		{
			float *cuda_pot;
			cudaMalloc((void**)&cuda_pot, GS.w*GS.z*sizeof(float));
			cudaMemcpy(cuda_pot, Potential, GS.w*GS.z*sizeof(float), cudaMemcpyHostToDevice);
			
			if(GS.y<=512)
			{
				dimBlock.x=GS.y;
				dimGrid.x=1;
			}
			else
			{
				dimBlock.x=512;
				dimGrid.x=GS.y/512;
				if(GS.y%512!=0)dimGrid.x++;
			}
			dimGrid.y=GS.z;
			if(BldBCatPlane[2])
			{
				dimGrid.y=GS.z-2;
			}
			printf("dimBlock [%d,%d,%d]\n",dimBlock.x,dimBlock.y,dimBlock.z);
			printf("dimGrid [%d,%d,%d]\n",dimGrid.x,dimGrid.y,dimGrid.z);
			KGOAtomsCoulBCYZOnCUDA<<<dimGrid, dimBlock>>>(cuda_rq,Nq, GS,eps, cuda_pot,cBldBCatPlane);
			cudaThreadSynchronize();
			
			cudaMemcpy(Potential, cuda_pot, GS.w*GS.z*sizeof(float), cudaMemcpyDeviceToHost);
			
			cudaFree(cuda_pot);
		}
		if(BldBCatPlane[1])
		{
			float *cuda_pot;
			cudaMalloc((void**)&cuda_pot, GS.w*GS.z*sizeof(float));
			cudaMemcpy(cuda_pot, Potential, GS.w*GS.z*sizeof(float), cudaMemcpyHostToDevice);
			
			if(GS.y<=512)
			{
				if(BldBCatPlane[0])
				{
					dimBlock.x=GS.x-2;
					dimGrid.x=1;
				}
				else
				{
					dimBlock.x=GS.x;
					dimGrid.x=1;
				}
			}
			else
			{
				if(BldBCatPlane[0])
				{
					dimBlock.x=512;
					dimGrid.x=(GS.x-2)/512;
					if((GS.x-2)%512!=0)dimGrid.x++;
				}
				else
				{
					dimBlock.x=512;
					dimGrid.x=GS.x/512;
					if(GS.x%512!=0)dimGrid.x++;
				}
				
			}
			dimGrid.y=GS.z;
			if(BldBCatPlane[2])
			{
				dimGrid.y=GS.z-2;
			}
			printf("dimBlock [%d,%d,%d]\n",dimBlock.x,dimBlock.y,dimBlock.z);
			printf("dimGrid [%d,%d,%d]\n",dimGrid.x,dimGrid.y,dimGrid.z);
			KGOAtomsCoulBCXZOnCUDA<<<dimGrid, dimBlock>>>(cuda_rq,Nq, GS,eps, cuda_pot,cBldBCatPlane);
			cudaThreadSynchronize();
			
			cudaMemcpy(Potential, cuda_pot, GS.w*GS.z*sizeof(float), cudaMemcpyDeviceToHost);
			
			cudaFree(cuda_pot);
		}
		cudaFree(cuda_rq);
	}
	else
	{
		int4 GS;
		GS.x=GridSize[0];
		GS.y=GridSize[1];
		GS.z=GridSize[2];
		GS.w=GridSize[0]*GridSize[1];
		
		int GrPnt;
		int iq;
		int i,j,k;
		int iStart,jStart,kStart;
		int iEnd,jEnd,kEnd;
		
		float r1;
		
		if(BldBCatPlane[2])
		{
			for(i=0;i<GS.x;i++)
				for(j=0;j<GS.y;j++)//XY
			{
				k=0;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
					
				k=GS.z-1;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
			}
		}
		if(BldBCatPlane[0])
		{
			jStart=0;
			kStart=0;
			jEnd=GS.y;
			kEnd=GS.z;
			if(BldBCatPlane[2])
			{
				kStart=1;
				kEnd=GS.z-1;
			}
			for(j=jStart;j<jEnd;j++)//YZ
				for(k=kStart;k<kEnd;k++)
			{
				i=0;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
					
				i=GS.x-1;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
			}
		}
		if(BldBCatPlane[1])
		{
			iStart=0;
			kStart=0;
			iEnd=GS.x;
			kEnd=GS.z;
			if(BldBCatPlane[2])
			{
				kStart=1;
				kEnd=GS.z-1;
			}
			if(BldBCatPlane[0])
			{
				iStart=1;
				iEnd=GS.x-1;
			}
			for(i=iStart;i<iEnd;i++)//ZX
				for(k=kStart;k<kEnd;k++)
			{
				j=0;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
					
				j=GS.y-1;
				GrPnt=i+j*GS.x+k*GS.w;
					
				for(iq=0;iq < Nq;iq++)
				{
					r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
					Potential[GrPnt]+=rq[iq].w/(eps*r1);
				}
			}
		}
	}
	return EXIT_SUCCESS;
}

#endif
