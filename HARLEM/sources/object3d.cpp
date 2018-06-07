/*! \file object3d.cpp
 
    Classes to define general displayed 3D objects

    \author Igor Kurnikov  
    \date 1999-2003
*/

#define OBJECT3D_CPP

#include "hastl.h"
#include "halinalg.h"
#include "object3d.h" 

#include "hamolview.h"//<mikola 30July06
#include "canvas3d.h"//<mikola 30July06
#include "math.h"//<mikola 13Sep07
#include "hasurface.h"//<mikola 26March08

Object3D::Object3D(const int new_obj_type, const char* new_name)
{
	obj_type = new_obj_type;
	SetObjName( new_name );

 	connect_flag= true;
	displayed_flag = true;
	transparency = 0.0;
}

Object3D::Object3D(const Object3D& obj_ref)
{
	obj_type = obj_ref.obj_type;
	SetObjName( obj_ref.GetObjName() );

	connect_flag= obj_ref.connect_flag;
	displayed_flag = true;
	transparency = 0.0;
}

Object3D::~Object3D()
{

}

int
Object3D::RotateX(double theta, const Vec3D& cnt )
{
	double cost = cos(theta);  
	double sint = sin(theta);
 
	HaMat_double Tr(3,3);

	Tr(1,1) = 1.0; Tr(1,2) = 0.0;   Tr(1,3) = 0.0;
    Tr(2,1) = 0.0; Tr(2,2) = cost;  Tr(2,3) = sint;
	Tr(3,1) = 0.0; Tr(3,2) = -sint; Tr(3,3) = cost;

	int ires = RotateObj(Tr,cnt);
	
	return ires;
}

int
Object3D::RotateY(double theta, const Vec3D& cnt )
{
	double cost = cos(theta);  
	double sint = sin(theta);
 
	HaMat_double Tr(3,3);

	Tr(1,1) = cost;  Tr(1,2) = 0.0; Tr(1,3) = sint;
	Tr(2,1) = 0.0;   Tr(2,2) = 1.0; Tr(2,3) = 0.0;
	Tr(3,1) = -sint; Tr(3,2) = 0.0; Tr(3,3) = cost;

	int ires = RotateObj(Tr,cnt);
	
	return ires;
}

int
Object3D::RotateZ(double theta, const Vec3D& cnt )
{
	double cost = cos(theta);  
	double sint = sin(theta);
 
	HaMat_double Tr(3,3);

	Tr(1,1) = cost; Tr(1,2) = -sint; Tr(1,3) = 0.0;
	Tr(2,1) = sint; Tr(2,2) = cost;  Tr(2,3) = 0.0;
	Tr(3,1) = 0.0;  Tr(3,2) = 0.0;   Tr(3,3) = 1.0;
	
	int ires = RotateObj(Tr,cnt);
	
	return ires;
}


int 
Object3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt)
{
    return 1;
}

int 
Object3D::Translate(const Vec3D& tr_vec )
{
	return 1;
}

int 
Object3D::SetTransparency(double trasp_new)
{
	return FALSE;
}
//////////////////////////////////////////////////////////////////////
BoxObj3D::BoxObj3D(const char* new_name,int r,int g,int b,float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1)
  :Object3D(OBJ3D_BOX,new_name)
{
  Color=HaColor::RegisterColor(r,g,b);
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
}
BoxObj3D::~BoxObj3D()
{
}
void BoxObj3D::SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1)
{
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
}
int BoxObj3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt )
{
  return 1;
}
int BoxObj3D::Translate( const Vec3D& tr_vec)
{
  return 1;
}
bool BoxObj3D::GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const
{
	MinX_v = x0*1.2;
	MinY_v = y0*1.2;
	MinZ_v = z0*1.2;
	MaxX_v = x1*1.2;
	MaxY_v = y1*1.2;
	MaxZ_v = z1*1.2;
	return true;
}
int BoxObj3D::DrawObj(HaMolView* mview)
{
  int ixadd=  mview->pCanv->XRange()/2 ;
  int iyadd=  mview->pCanv->YRange()/2 ;
  int izadd=  mview->ZOffset() ;
  double Scale=mview->Scale;
  double x_tr, y_tr, z_tr;

  //>this is a cheat, to trick the depth colors, otherwise box would not be pritty
  double TrickTheDepthColors;
  
  double fmax = 0.0;
  double fdist;
  mview->GetTransfCoord( x0, y0, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y0, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y1, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y1, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y0, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y0, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y1, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y1, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;

  if(fmax < 6.0) fmax = 6.0;
  double tmpDScale = 1.0/(2.0*sqrt(fmax));

  TrickTheDepthColors=tmpDScale/mview->DScale;

  if (TrickTheDepthColors>1.0)
	  TrickTheDepthColors=1.0;
  /*else
  {
	  mview->Scale=mview->Scale*TrickTheDepthColors;
	  Scale=mview->Scale;
  }*/
  //<end TrickTheDepthColors calculation
  int r1[3],r2[3],r3[3],r4[3],r5[3],r6[3],r7[3],r8[3];
  
  
  
  mview->GetTransfCoord( x0, y0, z0, x_tr, y_tr, z_tr);
  r1[0] = (int)(x_tr * Scale) + ixadd;
  r1[1] = (int)(y_tr * Scale) + iyadd;
  r1[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, z0, x_tr, y_tr, z_tr);
  r2[0] = (int)(x_tr * Scale) + ixadd;
  r2[1] = (int)(y_tr * Scale) + iyadd;
  r2[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, z0, x_tr, y_tr, z_tr);
  r3[0] = (int)(x_tr * Scale) + ixadd;
  r3[1] = (int)(y_tr * Scale) + iyadd;
  r3[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, z0, x_tr, y_tr, z_tr);
  r4[0] = (int)(x_tr * Scale) + ixadd;
  r4[1] = (int)(y_tr * Scale) + iyadd;
  r4[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y0, z1, x_tr, y_tr, z_tr);
  r5[0] = (int)(x_tr * Scale) + ixadd;
  r5[1] = (int)(y_tr * Scale) + iyadd;
  r5[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, z1, x_tr, y_tr, z_tr);
  r6[0] = (int)(x_tr * Scale) + ixadd;
  r6[1] = (int)(y_tr * Scale) + iyadd;
  r6[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, z1, x_tr, y_tr, z_tr);
  r7[0] = (int)(x_tr * Scale) + ixadd;
  r7[1] = (int)(y_tr * Scale) + iyadd;
  r7[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, z1, x_tr, y_tr, z_tr);
  r8[0] = (int)(x_tr * Scale) + ixadd;
  r8[1] = (int)(y_tr * Scale) + iyadd;
  r8[2] = (int)(z_tr * Scale) + izadd;
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r4[0],r4[1],r4[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r4[0],r4[1],r4[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r6[0],r6[1],r6[2],r8[0],r8[1],r8[2], Color,Color);
  mview->pCanv->ClipTwinVector(r7[0],r7[1],r7[2],r8[0],r8[1],r8[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r5[0],r5[1],r5[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r4[0],r4[1],r4[2],r8[0],r8[1],r8[2], Color,Color);

  return 1;
}
//////////////////////////////////////////////////////////////////////
MatrixObj3D::MatrixObj3D(const char* new_name,int r,int g,int b,float the_x0,float the_y0,float the_z0,float the_dx1,float the_dy1,float the_dz1,float the_dx2,float the_dy2,float the_dz2)
  :Object3D(OBJ3D_MATRIX,new_name)
{
  Color=HaColor::RegisterColor(r,g,b);
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  dx1=the_dx1;
  dy1=the_dy1;
  dz1=the_dz1;
  dx2=the_dx2;
  dy2=the_dy2;
  dz2=the_dz2;
  Ni=0;
  Nj=0;
  Color=HaColor::RegisterColor(r,g,b);
}
MatrixObj3D::~MatrixObj3D()
{
}
void MatrixObj3D::SetMatrixGuiders(float the_x0,float the_y0,float the_z0,float the_dx1,float the_dy1,float the_dz1,float the_dx2,float the_dy2,float the_dz2)
{
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  dx1=the_dx1;
  dy1=the_dy1;
  dz1=the_dz1;
  dx2=the_dx2;
  dy2=the_dy2;
  dz2=the_dz2;
  
}
void MatrixObj3D::SetColors(int Nx,int Ny,float * fmap,double elpot_high_val,double elpot_low_val)
{
  Ni=Nx;
  Nj=Ny;
  ColMat.newsize(Ni,Nj);
  
  int ic;
  int ncol = 128; // the number of colors in the scale

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

		
  double fspan = elpot_high_val - elpot_low_val;
	
  if( fspan < 1E-10) 
  {
    ErrorInMod("ElectrostMod::ColorMolSurfElPot()",
               "The difference between high and low boundary potential values is too small");
    return;
  }

  int i,j;

  for(i=0;i<Ni;i++)
  {
    for(j=0;j<Nj;j++)
    {
  
      double fval = fmap[i+j*Ni];
    
      double didx_col = 1.0 + ((double) (ncol-1) ) * ( (fval - elpot_low_val)/fspan);
  
      int idx_col = (int) didx_col;
      if( (didx_col - (double)idx_col ) > 0.5 )
        idx_col++;
  
      if( idx_col < 1) idx_col = 1;
      if( idx_col > ncol ) idx_col = ncol;
  
      col = col_num(idx_col);
      ColMat(i+1,j+1) = col;
    }
  }
}
int MatrixObj3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt )
{
  return 1;
}
int MatrixObj3D::Translate( const Vec3D& tr_vec)
{
  return 1;
}
int MatrixObj3D::DrawObj(HaMolView* mview)
{
  int ixadd=  mview->pCanv->XRange()/2 ;
  int iyadd=  mview->pCanv->YRange()/2 ;
  int izadd=  mview->ZOffset() ;
  double Scale=mview->Scale;
  
  //int r1[3],r2[3],r3[3],r4[3];
  float dx12=0.5*dx1,dy12=0.5*dy1,dz12=0.5*dz1;
  float dx22=0.5*dx2,dy22=0.5*dy2,dz22=0.5*dz2;
  double x_tr, y_tr, z_tr;
  double x1, y1, z1;
  int i,j;
  static Poly p;
  p.count=4;
  float nx=dy1*dz2-dy2*dz1,ny=dz1*dy2-dy1*dz2,nz=dx1*dy2-dx2*dy1;
  float zn1=nx*mview->Rot(3,1)  + ny*mview->Rot(3,2)  + nz*mview->Rot(3,3);
  /*Vec3D nrm,r1,r2;
  r1.SetX(dx1);
  r1.SetY(dy1);
  r1.SetZ(dz1);
  r2.SetX(dx2);
  r2.SetY(dy2);
  r2.SetZ(dz2);
  r1.normalize();
  r2.normalize();
  Vec3D::VecProduct(nrm,r1,r2);
  nrm.normalize();
  zn1=nrm.GetZ()*mview->Rot(3,1)  + nrm.GetY()*mview->Rot(3,2)  + nrm.GetZ()*mview->Rot(3,3);
  float fColourDepth=(int)ColourDepth;
  if(zn1<0.0)zn1=-zn1;
  if(zn1>0.99)zn1=0.99;*/
  for(i=0;i<Ni;i++)
  {
    for(j=0;j<Nj;j++)
    {
      x1=x0+(float)i*dx1+(float)j*dx2;
      y1=y0+(float)i*dy1+(float)j*dy2;
      z1=z0+(float)i*dz1+(float)j*dz2;
      mview->GetTransfCoord( x1-dx12-dx22, y1-dy12-dy22, z1-dz12-dz22, x_tr, y_tr, z_tr);
      p.v[0].x = (int)(x_tr * Scale) + ixadd;
      p.v[0].y = (int)(y_tr * Scale) + iyadd;
      p.v[0].z = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( x1+dx12-dx22, y1+dy12-dy22, z1+dz12-dz22, x_tr, y_tr, z_tr);
      p.v[1].x = (int)(x_tr * Scale) + ixadd;
      p.v[1].y = (int)(y_tr * Scale) + iyadd;
      p.v[1].z = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( x1-dx12+dx22, y1-dy12+dy22, z1-dz12+dz22, x_tr, y_tr, z_tr);
      p.v[3].x = (int)(x_tr * Scale) + ixadd;
      p.v[3].y = (int)(y_tr * Scale) + iyadd;
      p.v[3].z = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( x1+dx12+dx22, y1+dy12+dy22, z1+dz12+dz22, x_tr, y_tr, z_tr);
      p.v[2].x = (int)(x_tr * Scale) + ixadd;
      p.v[2].y = (int)(y_tr * Scale) + iyadd;
      p.v[2].z = (int)(z_tr * Scale) + izadd;
      
      /*if( have_norms) 
      {
        p.v[0].inten = ColMat(i+1,j+1) + ColourDepth*zn1;
        p.v[1].inten = ColMat(i+1,j+1) + ColourDepth*zn2; 
        p.v[2].inten = ColMat(i+1,j+1) + ColourDepth*zn3;
      }*/
//       p.v[0].inten = ColMat(i+1,j+1) + (int)(fColourDepth*zn1);
//       p.v[1].inten = ColMat(i+1,j+1) + (int)(fColourDepth*zn1);
//       p.v[2].inten = ColMat(i+1,j+1) + (int)(fColourDepth*zn1);
//       p.v[3].inten = ColMat(i+1,j+1) + (int)(fColourDepth*zn1);
      p.v[0].inten = ColMat(i+1,j+1)+31;
      p.v[1].inten = ColMat(i+1,j+1)+31;
      p.v[2].inten = ColMat(i+1,j+1)+31;
      p.v[3].inten = ColMat(i+1,j+1)+31;

      mview->pCanv->ClipPolygon( &p,0.0 );
      //Color=ColMat(i+1,j+1);
      /*mview->pCanv->ClipTwinVector( p.v[0].x,p.v[0].y,p.v[0].z, p.v[1].x,p.v[1].y,p.v[1].z, Color,Color);
      mview->pCanv->ClipTwinVector( p.v[0].x,p.v[0].y,p.v[0].z, p.v[2].x,p.v[2].y,p.v[2].z, Color,Color);
      mview->pCanv->ClipTwinVector( p.v[1].x,p.v[1].y,p.v[1].z, p.v[2].x,p.v[2].y,p.v[2].z, Color,Color);
      mview->pCanv->ClipTwinVector( p.v[0].x,p.v[0].y,p.v[0].z, p.v[3].x,p.v[3].y,p.v[3].z, Color,Color);
      mview->pCanv->ClipTwinVector( p.v[0].x,p.v[0].y,p.v[0].z, p.v[3].x,p.v[3].y,p.v[3].z, Color,Color);
      mview->pCanv->ClipTwinVector( p.v[3].x,p.v[3].y,p.v[3].z, p.v[2].x,p.v[2].y,p.v[2].z, Color,Color);*/
      
    }
  }
  
  
  return 1;
}
//////////////////////////////////////////////////////////////////////
PlaneViewOfHaField3D::PlaneViewOfHaField3D(HaField3D *_field, const char* new_name, int OwnerOfData)
	:Object3D(OBJ3D_PlaneViewerOfHaField3D,new_name)
{
	//int r=255, g=255, b=255;
	//Color=HaColor::RegisterColor(r,g,b);
	Plane=PlaneXY;
	field=_field;
	Level=field->GetNz()/2;
	HideZeroValues=false;
	OwnerOfData=OwnerOfData;
	
	float fValueMin,fValueMax;
	field->GetMinMaxValue(&fValueMin,&fValueMax);
	ValueMin=fValueMin;ValueMax=fValueMax;
	SetColorsOfPlane();
}
PlaneViewOfHaField3D::~PlaneViewOfHaField3D()
{
	if(OwnerOfData==1)
	{
		delete field;
	}
}
void PlaneViewOfHaField3D::SetMinMax(double m_Min,double m_Max)
{
	ValueMin=m_Min;ValueMax=m_Max;
	SetColorsOfPlane();
}
void PlaneViewOfHaField3D::SetLevel(int NewLevel)
{
	if(Plane==PlaneXY)
	{
		if((NewLevel>=0)&&(NewLevel<=field->GetNz()-1))
			Level=NewLevel;
	}
	else if(Plane==PlaneYZ)
	{
		if((NewLevel>=0)&&(NewLevel<=field->GetNx()-1))
			Level=NewLevel;
	}
	else
	{
		if((NewLevel>=0)&&(NewLevel<=field->GetNy()-1))
			Level=NewLevel;
	}
	//printf("Plane=%d Level=%d NewLevel=%d [%d %d %d]\n",Plane,Level,NewLevel,field->GetNx(),field->GetNy(),field->GetNz());
	SetColorsOfPlane();
}
void PlaneViewOfHaField3D::SetHideZeroValues(bool newHideZeroValues)
{
	HideZeroValues=newHideZeroValues;
	SetColorsOfPlane();
}
void PlaneViewOfHaField3D::SetPlaneXY()
{
	Plane=PlaneXY;
	SetLevel(field->GetNz()/2);
}
void PlaneViewOfHaField3D::SetPlaneYZ()
{
	Plane=PlaneYZ;
	SetLevel(field->GetNx()/2);
}
void PlaneViewOfHaField3D::SetPlaneZX()
{
	Plane=PlaneZX;
	SetLevel(field->GetNy()/2);
}
void PlaneViewOfHaField3D::SetPlane(int newPlane)
{
	if(newPlane==PlaneXY)
		SetPlaneXY();
	if(newPlane==PlaneYZ)
		SetPlaneYZ();
	if(newPlane==PlaneZX)
		SetPlaneZX();
}
void PlaneViewOfHaField3D::SetColorsOfPlane()
{
	x0=field->GetXmin();
	y0=field->GetYmin();
	z0=field->GetZmin();
	
	switch(Plane)
	{
		case PlaneXY:
		{
			Ni=field->GetNx();
			Nj=field->GetNy();
			
			dx1=field->GetXstep();
			dy1=0.0;
			dz1=0.0;
			dx2=0.0;
			dy2=field->GetYstep();
			dz2=0.0;
			z0+=field->GetZstep()*Level;
			
			break;
		}
		case PlaneYZ:
		{
			Ni=field->GetNy();
			Nj=field->GetNz();
			
			dx1=0.0;
			dy1=field->GetYstep();
			dz1=0.0;
			dx2=0.0;
			dy2=0.0;
			dz2=field->GetZstep();
			x0+=field->GetXstep()*Level;
			
			break;
		}
		case PlaneZX:
		{
			Ni=field->GetNz();
			Nj=field->GetNx();
			
			dx1=0.0;
			dy1=0.0;
			dz1=field->GetZstep();
			dx2=field->GetXstep();
			dy2=0.0;
			dz2=0.0;
			y0+=field->GetYstep()*Level;
			
			break;
		}
	}
	//printf("Plane=%d Level=%d r0=[%f %f %f] e0=[%f %f %f] e1=[%f %f %f]\n",Plane,Level,x0,y0,z0,dx1,dy1,dz1,dx2,dy2,dz2);
	ColMat.newsize(Ni,Nj);
	
  //int Nx,int Ny,float * fmap,double elpot_high_val,double elpot_low_val
	int ic;
	int ncol = 128; // the number of colors in the scale

	int nmid = ncol/2 + 1;
	int nspan2 = ncol/2;

	HaVec_int col_num(ncol);

	int r,g,b;

	int col;

	//prepare color shades interpolated between red -> white -> blue
	for(ic = 1; ic <= ncol; ic++) 
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

		
	double fspan = ValueMax - ValueMin;
	
	if( fspan < 1E-10) 
	{
		ErrorInMod("ElectrostMod::ColorMolSurfElPot()",
							 "The difference between high and low boundary potential values is too small");
		ValueMax = 1.05*ValueMax;
		ValueMin = 0.95*ValueMin;
		fspan = ValueMax - ValueMin;
	}

	int i,j;

	float *Map=field->GetFieldPtr();
	int Nx=field->GetNx();
	int Nxy=field->GetNx()*field->GetNy();
	
	double fval;
	for(i=0;i<Ni;i++)
	{
		for(j=0;j<Nj;j++)
		{
			//switch(Plane)
			//{
			//	case PlaneXY:
			if(Plane==PlaneXY)
					fval = Map[i+j*Nx+Level*Nxy];
			else if(Plane==PlaneYZ)
				fval = Map[Level+i*Nx+j*Nxy];
			else
				fval = Map[j+Level*Nx+i*Nxy];
			double didx_col = 1.0 + ((double) (ncol-1) ) * ( (fval - ValueMin)/fspan);
			
			int idx_col = (int) didx_col;
			if( (didx_col - (double)idx_col ) > 0.5 )
				idx_col++;
			
			if( idx_col < 1) idx_col = 1;
			if( idx_col > ncol ) idx_col = ncol;
			
			col = col_num(idx_col);
			ColMat(i+1,j+1) = col;
			if(HideZeroValues&&fval==0.0)ColMat(i+1,j+1) = 0;
		}
	}
	
}
int PlaneViewOfHaField3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt )
{
	return 1;
}
int PlaneViewOfHaField3D::Translate( const Vec3D& tr_vec)
{
	return 1;
}
int PlaneViewOfHaField3D::DrawObj(HaMolView* mview)
{
	int ixadd=  mview->pCanv->XRange()/2 ;
	int iyadd=  mview->pCanv->YRange()/2 ;
	int izadd=  mview->ZOffset() ;
	double Scale=mview->Scale;
  
  //int r1[3],r2[3],r3[3],r4[3];
	float dx12=0.5*dx1,dy12=0.5*dy1,dz12=0.5*dz1;
	float dx22=0.5*dx2,dy22=0.5*dy2,dz22=0.5*dz2;
	double x_tr, y_tr, z_tr;
	double x1, y1, z1;
	int i,j;
	static Poly p;
	p.count=4;
	float nx=dy1*dz2-dy2*dz1,ny=dz1*dy2-dy1*dz2,nz=dx1*dy2-dx2*dy1;
	float zn1=nx*mview->Rot(3,1)  + ny*mview->Rot(3,2)  + nz*mview->Rot(3,3);
	
	for(i=0;i<Ni;i++)
	{
		for(j=0;j<Nj;j++)
		{
			if(ColMat(i+1,j+1)>0)
			{
				x1=x0+(float)i*dx1+(float)j*dx2;
				y1=y0+(float)i*dy1+(float)j*dy2;
				z1=z0+(float)i*dz1+(float)j*dz2;
				mview->GetTransfCoord( x1-dx12-dx22, y1-dy12-dy22, z1-dz12-dz22, x_tr, y_tr, z_tr);
				p.v[0].x = (int)(x_tr * Scale) + ixadd;
				p.v[0].y = (int)(y_tr * Scale) + iyadd;
				p.v[0].z = (int)(z_tr * Scale) + izadd;
				mview->GetTransfCoord( x1+dx12-dx22, y1+dy12-dy22, z1+dz12-dz22, x_tr, y_tr, z_tr);
				p.v[1].x = (int)(x_tr * Scale) + ixadd;
				p.v[1].y = (int)(y_tr * Scale) + iyadd;
				p.v[1].z = (int)(z_tr * Scale) + izadd;
				mview->GetTransfCoord( x1-dx12+dx22, y1-dy12+dy22, z1-dz12+dz22, x_tr, y_tr, z_tr);
				p.v[3].x = (int)(x_tr * Scale) + ixadd;
				p.v[3].y = (int)(y_tr * Scale) + iyadd;
				p.v[3].z = (int)(z_tr * Scale) + izadd;
				mview->GetTransfCoord( x1+dx12+dx22, y1+dy12+dy22, z1+dz12+dz22, x_tr, y_tr, z_tr);
				p.v[2].x = (int)(x_tr * Scale) + ixadd;
				p.v[2].y = (int)(y_tr * Scale) + iyadd;
				p.v[2].z = (int)(z_tr * Scale) + izadd;
				
				p.v[0].inten = ColMat(i+1,j+1)+31;
				p.v[1].inten = ColMat(i+1,j+1)+31;
				p.v[2].inten = ColMat(i+1,j+1)+31;
				p.v[3].inten = ColMat(i+1,j+1)+31;
				
				mview->pCanv->ClipPolygon( &p,0.0 );
			}
		}
	}
  
  
	return 1;
}
bool PlaneViewOfHaField3D::GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const
{
	MinX_v = field->GetXmin()*1.2;
	MinY_v = field->GetYmin()*1.2;
	MinZ_v = field->GetZmin()*1.2;
	MaxX_v = field->GetXmax()*1.2;
	MaxY_v = field->GetYmax()*1.2;
	MaxZ_v = field->GetZmax()*1.2;
	return true;
}
//////////////////////////////////////////////////////////////////////
MembraneZObj3D::MembraneZObj3D(const char* new_name,int r,int g,int b,float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1,float the_mem_z1,float the_mem_z2)
  :Object3D(OBJ3D_MEMBRANEZ,new_name)
{
  Color=HaColor::RegisterColor(r,g,b);
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
  mem_z1=the_mem_z1;
  mem_z2=the_mem_z2;
  hstep=20.0;
}
MembraneZObj3D::~MembraneZObj3D()
{
}
void MembraneZObj3D::SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1)
{
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
  hstep=(x1-x0)/5.0;
}
void MembraneZObj3D::SetMembraneZ(float the_mem_z1,float the_mem_z2)
{
  mem_z1=the_mem_z1;
  mem_z2=the_mem_z2;
}
int MembraneZObj3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt )
{
  return 1;
}
int MembraneZObj3D::Translate( const Vec3D& tr_vec)
{
  return 1;
}
int MembraneZObj3D::DrawObj(HaMolView* mview)
{
  int ixadd=  mview->pCanv->XRange()/2 ;
  int iyadd=  mview->pCanv->YRange()/2 ;
  int izadd=  mview->ZOffset() ;
  double Scale=mview->Scale;
  
  float fx0,fy0,fx1,fy1,fz;
  int r1[3],r2[3],r3[3],r4[3],r5[3],r6[3],r7[3],r8[3];
  
  double x_tr, y_tr, z_tr;
  
  //big box
  mview->GetTransfCoord( x0, y0, mem_z1, x_tr, y_tr, z_tr);
  r1[0] = (int)(x_tr * Scale) + ixadd;
  r1[1] = (int)(y_tr * Scale) + iyadd;
  r1[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, mem_z1, x_tr, y_tr, z_tr);
  r2[0] = (int)(x_tr * Scale) + ixadd;
  r2[1] = (int)(y_tr * Scale) + iyadd;
  r2[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, mem_z1, x_tr, y_tr, z_tr);
  r3[0] = (int)(x_tr * Scale) + ixadd;
  r3[1] = (int)(y_tr * Scale) + iyadd;
  r3[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, mem_z1, x_tr, y_tr, z_tr);
  r4[0] = (int)(x_tr * Scale) + ixadd;
  r4[1] = (int)(y_tr * Scale) + iyadd;
  r4[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y0, mem_z2, x_tr, y_tr, z_tr);
  r5[0] = (int)(x_tr * Scale) + ixadd;
  r5[1] = (int)(y_tr * Scale) + iyadd;
  r5[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, mem_z2, x_tr, y_tr, z_tr);
  r6[0] = (int)(x_tr * Scale) + ixadd;
  r6[1] = (int)(y_tr * Scale) + iyadd;
  r6[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, mem_z2, x_tr, y_tr, z_tr);
  r7[0] = (int)(x_tr * Scale) + ixadd;
  r7[1] = (int)(y_tr * Scale) + iyadd;
  r7[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, mem_z2, x_tr, y_tr, z_tr);
  r8[0] = (int)(x_tr * Scale) + ixadd;
  r8[1] = (int)(y_tr * Scale) + iyadd;
  r8[2] = (int)(z_tr * Scale) + izadd;
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r4[0],r4[1],r4[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r4[0],r4[1],r4[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r6[0],r6[1],r6[2],r8[0],r8[1],r8[2], Color,Color);
  mview->pCanv->ClipTwinVector(r7[0],r7[1],r7[2],r8[0],r8[1],r8[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r5[0],r5[1],r5[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r4[0],r4[1],r4[2],r8[0],r8[1],r8[2], Color,Color);
  //net on top
  int i,N,j;
  
  N=(int)((x1-x0)/hstep);
  fx0,fy0,fx1,fy1;
  fy0=y0;
  fy1=y1;
  
  fz=mem_z1;
  for(j=0;j<2;j++)
  {
    for(i=1;i<=N;i++)
    {
      fx0=hstep*i+x0;
      mview->GetTransfCoord( fx0, fy0, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( fx0, fy1, fz, x_tr, y_tr, z_tr);
      r2[0] = (int)(x_tr * Scale) + ixadd;
      r2[1] = (int)(y_tr * Scale) + iyadd;
      r2[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
      
      mview->GetTransfCoord( fx0, fy0, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( fx0, fy1, fz, x_tr, y_tr, z_tr);
      r2[0] = (int)(x_tr * Scale) + ixadd;
      r2[1] = (int)(y_tr * Scale) + iyadd;
      r2[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
    }
    fz=mem_z2;
  }
  N=(int)((y1-y0)/hstep);
  fx0=x0;
  fx1=x1;
  fz=mem_z1;
  for(j=0;j<2;j++)
  {
    for(i=1;i<=N;i++)
    {
      fy0=hstep*i+y0;
      mview->GetTransfCoord( fx0, fy0, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( fx1, fy0, fz, x_tr, y_tr, z_tr);
      r2[0] = (int)(x_tr * Scale) + ixadd;
      r2[1] = (int)(y_tr * Scale) + iyadd;
      r2[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
      
      mview->GetTransfCoord( fx0, fy0, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( fx1, fy0, fz, x_tr, y_tr, z_tr);
      r2[0] = (int)(x_tr * Scale) + ixadd;
      r2[1] = (int)(y_tr * Scale) + iyadd;
      r2[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
    }
    fz=mem_z2;
  }
  return 1;
}
//////////////////////////////////////////////////////////////////////
TubeObj3D::TubeObj3D(const char* new_name,int r,int g,int b,float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1)
  :Object3D(OBJ3D_MEMBRANE_TUBE,new_name)
{
  Color=HaColor::RegisterColor(r,g,b);
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
  hstep=5.0;
  phistep=PI/6.0;
  Style=0;
}
TubeObj3D::~TubeObj3D()
{
}
void TubeObj3D::SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1)
{
  x0=the_x0;
  y0=the_y0;
  z0=the_z0;
  x1=the_x1;
  y1=the_y1;
  z1=the_z1;
  hstep=(x1-x0)/5.0;
  phistep=PI/6.0;
}
void TubeObj3D::SetTube3d(float the_mem_z1,float the_mem_z2,float the_mem_x,float the_mem_y,float the_R1,float the_R2)
{
  mem_z1=the_mem_z1;
  mem_z2=the_mem_z2;
  mem_x=the_mem_x;
  mem_y=the_mem_y;
  R1=the_R1;
  R2=the_R2;
}
int TubeObj3D::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt )
{
  return 1;
}
int TubeObj3D::Translate( const Vec3D& tr_vec)
{
  return 1;
}
bool TubeObj3D::GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const
{
	MinX_v = x0*1.2;
	MinY_v = y0*1.2;
	MinZ_v = z0*1.2;
	MaxX_v = x1*1.2;
	MaxY_v = y1*1.2;
	MaxZ_v = z1*1.2;
	return true;
}
int TubeObj3D::DrawObj(HaMolView* mview)
{
  int ixadd=  mview->pCanv->XRange()/2 ;
  int iyadd=  mview->pCanv->YRange()/2 ;
  int izadd=  mview->ZOffset() ;
  double Scale=mview->Scale;
  
  float fx0,fy0,fx1,fy1;
  int r1[3],r2[3],r3[3],r4[3],r5[3],r6[3],r7[3],r8[3];
  
  double x_tr, y_tr, z_tr;

  //>this is a cheat, to trick the depth colors, otherwise box would not be pritty
  /*double TrickTheDepthColors;
  
  double fmax = 0.0;
  double fdist;
  mview->GetTransfCoord( x0, y0, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y0, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y1, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y1, z0, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y0, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y0, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x0, y1, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;
  mview->GetTransfCoord( x1, y1, z1, x_tr, y_tr, z_tr);
  fdist = x_tr*x_tr + y_tr*y_tr + z_tr*z_tr;
  if( fdist > fmax )fmax = fdist;

  if(fmax < 6.0) fmax = 6.0;
  double tmpDScale = 1.0/(2.0*sqrt(fmax));

  TrickTheDepthColors=tmpDScale/mview->DScale;

  if (TrickTheDepthColors>1.0)TrickTheDepthColors=1.0;*/
  //<end TrickTheDepthColors calculation

  float fx,fy,fz;
  
  //big box
  mview->GetTransfCoord( x0, y0, mem_z1, x_tr, y_tr, z_tr);
  r1[0] = (int)(x_tr * Scale) + ixadd;
  r1[1] = (int)(y_tr * Scale) + iyadd;
  r1[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, mem_z1, x_tr, y_tr, z_tr);
  r2[0] = (int)(x_tr * Scale) + ixadd;
  r2[1] = (int)(y_tr * Scale) + iyadd;
  r2[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, mem_z1, x_tr, y_tr, z_tr);
  r3[0] = (int)(x_tr * Scale) + ixadd;
  r3[1] = (int)(y_tr * Scale) + iyadd;
  r3[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, mem_z1, x_tr, y_tr, z_tr);
  r4[0] = (int)(x_tr * Scale) + ixadd;
  r4[1] = (int)(y_tr * Scale) + iyadd;
  r4[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y0, mem_z2, x_tr, y_tr, z_tr);
  r5[0] = (int)(x_tr * Scale) + ixadd;
  r5[1] = (int)(y_tr * Scale) + iyadd;
  r5[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y0, mem_z2, x_tr, y_tr, z_tr);
  r6[0] = (int)(x_tr * Scale) + ixadd;
  r6[1] = (int)(y_tr * Scale) + iyadd;
  r6[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x0, y1, mem_z2, x_tr, y_tr, z_tr);
  r7[0] = (int)(x_tr * Scale) + ixadd;
  r7[1] = (int)(y_tr * Scale) + iyadd;
  r7[2] = (int)(z_tr * Scale) + izadd;
  mview->GetTransfCoord( x1, y1, mem_z2, x_tr, y_tr, z_tr);
  r8[0] = (int)(x_tr * Scale) + ixadd;
  r8[1] = (int)(y_tr * Scale) + iyadd;
  r8[2] = (int)(z_tr * Scale) + izadd;
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r2[0],r2[1],r2[2], Color,Color);
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r4[0],r4[1],r4[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r4[0],r4[1],r4[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r6[0],r6[1],r6[2],r8[0],r8[1],r8[2], Color,Color);
  mview->pCanv->ClipTwinVector(r7[0],r7[1],r7[2],r8[0],r8[1],r8[2], Color,Color);
  
  mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r5[0],r5[1],r5[2], Color,Color);
  mview->pCanv->ClipTwinVector(r2[0],r2[1],r2[2],r6[0],r6[1],r6[2], Color,Color);
  mview->pCanv->ClipTwinVector(r3[0],r3[1],r3[2],r7[0],r7[1],r7[2], Color,Color);
  mview->pCanv->ClipTwinVector(r4[0],r4[1],r4[2],r8[0],r8[1],r8[2], Color,Color);
  //net on top
  float phi=0.0,dphi=PI/32.0;
  int j;
  
  
  //circles
  fz=mem_z1;
  for(j=0;j<2;j++)
  {
    fx=R1*cos(phi)+mem_x;
    fy=R1*sin(phi)+mem_y;
    mview->GetTransfCoord( fx, fy, fz, x_tr, y_tr, z_tr);
    r3[0] = (int)(x_tr * Scale) + ixadd;
    r3[1] = (int)(y_tr * Scale) + iyadd;
    r3[2] = (int)(z_tr * Scale) + izadd;
    r5[0] = r3[0];
    r5[1] = r3[1];
    r5[2] = r3[2];
    for(phi=dphi;phi<2*PI;phi+=dphi)
    {
      fx=R1*cos(phi)+mem_x;
      fy=R1*sin(phi)+mem_y;
      
      mview->GetTransfCoord( fx, fy, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
      r3[0] = r1[0];
      r3[1] = r1[1];
      r3[2] = r1[2];
    }
    mview->pCanv->ClipTwinVector(r5[0],r5[1],r5[2],r3[0],r3[1],r3[2], Color,Color);
    fz=mem_z2;
  }
  //lines beside circles
  dphi=PI/16.0;
  if(R1!=0.0)
  {
	for(phi=0.0;phi<2*PI;phi+=dphi)
	{
		fx=R1*cos(phi)+mem_x;
		fy=R1*sin(phi)+mem_y;
    
		mview->GetTransfCoord( fx, fy, mem_z1, x_tr, y_tr, z_tr);
		r1[0] = (int)(x_tr * Scale) + ixadd;
		r1[1] = (int)(y_tr * Scale) + iyadd;
		r1[2] = (int)(z_tr * Scale) + izadd;
		mview->GetTransfCoord( fx, fy, mem_z2, x_tr, y_tr, z_tr);
		r3[0] = (int)(x_tr * Scale) + ixadd;
		r3[1] = (int)(y_tr * Scale) + iyadd;
		r3[2] = (int)(z_tr * Scale) + izadd;
		mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
	}
  }
  //radial lines
  float phi1,phi2,phi3,phi4,fR;
  phi1=atan((y1-mem_y)/(x1-mem_x));
  phi4=2*PI+atan((y0-mem_y)/(x1-mem_x));
  phi2=PI+atan((y1-mem_y)/(x0-mem_x));
  phi3=PI+atan((y0-mem_y)/(x0-mem_x));
  fz=mem_z1;
  for(j=0;j<2;j++)
  {
    for(phi=0.0;phi<2.0*PI;phi+=dphi)
    {
      if(phi>=phi4||phi<phi1)
      {
        fx0=R1*cos(phi)+mem_x;
        fy0=R1*sin(phi)+mem_y;
        fR=R2;//(y1-mem_y)/sin(phi);
        fx1=x1;
        fy1=(x1-mem_x)*tan(phi)+mem_y;
      }
      if(phi>=phi1&&phi<phi2)
      {
        fx0=R1*cos(phi)+mem_x;
        fy0=R1*sin(phi)+mem_y;
        fR=R2;//(y1-mem_y)/sin(phi);
        fx1=(y1-mem_y)/tan(phi)+mem_x;
        fy1=y1;
      }
      if(phi>=phi2&&phi<phi3)
      {
        fx0=R1*cos(phi)+mem_x;
        fy0=R1*sin(phi)+mem_y;
        fR=R2;//(y1-mem_y)/sin(phi);
        fx1=x0;
        fy1=(x0-mem_x)*tan(phi)+mem_y;
      }
      if(phi>=phi3&&phi<phi4)
      {
        fx0=R1*cos(phi)+mem_x;
        fy0=R1*sin(phi)+mem_y;
        fR=R2;//(y1-mem_y)/sin(phi);
        fx1=(y0-mem_y)/tan(phi)+mem_x;
        fy1=y0;
      }
      mview->GetTransfCoord( fx0, fy0, fz, x_tr, y_tr, z_tr);
      r1[0] = (int)(x_tr * Scale) + ixadd;
      r1[1] = (int)(y_tr * Scale) + iyadd;
      r1[2] = (int)(z_tr * Scale) + izadd;
      mview->GetTransfCoord( fx1, fy1, fz, x_tr, y_tr, z_tr);
      r3[0] = (int)(x_tr * Scale) + ixadd;
      r3[1] = (int)(y_tr * Scale) + iyadd;
      r3[2] = (int)(z_tr * Scale) + izadd;
      mview->pCanv->ClipTwinVector(r1[0],r1[1],r1[2],r3[0],r3[1],r3[2], Color,Color);
    }
    fz=mem_z2;
  }
  return 1;
}
