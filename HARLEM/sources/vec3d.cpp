/*! vec3d.cpp

  Classes to define a vector in 3D space 

  \author Igor Kurnikov  

   \date 2000-2003
*/

#define VEC3D_CPP

#include "haconst.h"
#include <memory>
#include <float.h>
#include <math.h>
#include "haio.h"
#include "vec3d.h"
#include "halinalg.h"
#include "haatgroup.h"

Vec3D::Vec3D() 
{
	pos[0] = 0.0;
	pos[1] = 0.0;
	pos[2] = 0.0;
}
Vec3D::Vec3D(double x,double y,double z)
{
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
}

Vec3D::~Vec3D()
{
	
}

Vec3D& Vec3D::operator=(const Vec3D& ref_vec)
{
	this->SetX( ref_vec.GetX());
	this->SetY( ref_vec.GetY());
	this->SetZ( ref_vec.GetZ());

	return(*this);
}

void Vec3D::SetCoordFrom(const Vec3D& ref_vec)
{
	this->SetX( ref_vec.GetX());
	this->SetY( ref_vec.GetY());
	this->SetZ( ref_vec.GetZ());
}

Vec3D operator*(const double factor, const Vec3D& v)
{
	Vec3D v_out;
	v_out = v;
	v_out.Scale(factor);
	return v_out;
}

Vec3D operator+(const Vec3D& v1,const Vec3D& v2)
{
	Vec3D v_out;
	Vec3D::sum(v_out,v1,v2);
	return v_out;
}

Vec3D operator-(const Vec3D& v1,const Vec3D& v2)
{
	Vec3D v_out;
	Vec3D::diff(v_out,v1,v2);
	return v_out;
}

Vec3D& Vec3D::operator +=(const Vec3D& v)
{
	pos[0] += v.GetX();
	pos[1] += v.GetY();
	pos[2] += v.GetZ();
	return(*this);
}


void Vec3D::SetX(const double x_new)
// Assumes X_new in Bohr
{
	pos[0]=x_new;
}

void Vec3D::SetY(const double y_new)
{
	pos[1]=y_new;
}

void Vec3D::SetZ(const double z_new)
{
	pos[2]=z_new;
}

void Vec3D::SetX_Bohr( double x_new )
// Assumes X_new in Bohr
{
	pos[0]= x_new * BOHR_TO_ANG;
}

void Vec3D::SetY_Bohr( double y_new )
{
	pos[1]= y_new * BOHR_TO_ANG;
}

void Vec3D::SetZ_Bohr( double z_new)
{
	pos[2]= z_new * BOHR_TO_ANG;
}

void Vec3D::SetX_Ang( double x_new)
// Assumes x_new is in Angstroms
{
	pos[0]=x_new;
}

void Vec3D::SetY_Ang( double y_new)
{
	pos[1]=y_new;
}

void Vec3D::SetZ_Ang( double z_new)
{
	pos[2]=z_new;
}

int Vec3D::SetFromStr(const char* str)
{
	istrstream is(str);
	int i;

	for( i = 0; i < 3; i++)
	{
		is >> pos[i];
		if(is.fail())
		{
			ErrorInMod(" Vec3D::SetFromStr()", 
				       " Failed to read a vector element ");
			return FALSE;
		}
	}
	return TRUE;
}

double Vec3D::length() const 
{ 
   return sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]); 
}

double Vec3D::length2() const 
{ 
   return (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]); 
}


int Vec3D::IsClose(const Vec3D& pt, double crit ) const
{
	if( fabs(pt[0] - (*this)[0]) > crit ) return FALSE;
	if( fabs(pt[1] - (*this)[1]) > crit ) return FALSE;
	if( fabs(pt[2] - (*this)[2]) > crit ) return FALSE;
	return TRUE;
}


void Vec3D::Scale(double scale)
{
	pos[0]*= scale;
	pos[1]*= scale;
	pos[2]*= scale;
}

int Vec3D::normalize()
{
	double len = this->length();
	if(len < DBL_EPSILON) return FALSE;

	this->Scale(1.0/len);
	return TRUE;
}

int Vec3D::RotatePt(Vec3D& still_pt, Vec3D& n, double ca, double sa)
{
    double dx = this->GetX() - still_pt.GetX();
	double dy = this->GetY() - still_pt.GetY();
	double dz = this->GetZ() - still_pt.GetZ();


	double atx = (n[0]*n[0] + (1-n[0]*n[0])*ca)*dx+
	             (n[0]*n[1]*(1-ca)  -  n[2]*sa)*dy+ 
				 (n[0]*n[2]*(1-ca)  +  n[1]*sa)*dz;

	double aty = (n[0]*n[1]*(1-ca)  +  n[2]*sa)*dx+
		         (n[1]*n[1] + (1-n[1]*n[1])*ca)*dy+
                 (n[1]*n[2]*(1-ca)  -  n[0]*sa)*dz;

	double atz = (n[0]*n[2]*(1-ca) - n[1]*sa)*dx+
		         (n[1]*n[2]*(1-ca) + n[0]*sa)*dy+
                 (n[2]*n[2] + (1-n[2]*n[2])*ca)*dz;	

	this->SetX(still_pt.GetX() + atx);
    this->SetY(still_pt.GetY() + aty);
	this->SetZ(still_pt.GetZ() + atz);

	return FALSE;
}

int Vec3D::Rotate(Vec3D& n, double ca, double sa)
{
	double atx = (n[0]*n[0] + (1-n[0]*n[0])*ca)*GetX() +
	             (n[0]*n[1]*(1-ca)  -  n[2]*sa)*GetY() + 
				 (n[0]*n[2]*(1-ca)  +  n[1]*sa)*GetZ();

	double aty = (n[0]*n[1]*(1-ca)  +  n[2]*sa)*GetX()+
		         (n[1]*n[1] + (1-n[1]*n[1])*ca)*GetY()+
                 (n[1]*n[2]*(1-ca)  -  n[0]*sa)*GetZ();

	double atz = (n[0]*n[2]*(1-ca) - n[1]*sa)*GetX()+
		         (n[1]*n[2]*(1-ca) + n[0]*sa)*GetY()+
                 (n[2]*n[2] + (1-n[2]*n[2])*ca)*GetZ();	

	this->SetX(atx);
    this->SetY(aty);
	this->SetZ(atz);

	return FALSE;
}

double Vec3D::CalcDistance(const Vec3D* atm1, const Vec3D* atm2, const CoordUnits units )
{
    double dx,dy,dz;
    double dist2;

    dx = atm1->GetX() - atm2->GetX();
    dy = atm1->GetY() - atm2->GetY();
    dz = atm1->GetZ() - atm2->GetZ();

	dist2= dx*dx + dy*dy + dz*dz;
 	if(units == ANGSTROM_U)
		return(sqrt(dist2));
	else
		return(sqrt(dist2)/BOHR_TO_ANG);

}

double Vec3D::CalcDistanceSq(const Vec3D* atm1, const Vec3D* atm2, const CoordUnits units )
{
    double dx,dy,dz;
    double dist2;

    dx = atm1->GetX() - atm2->GetX();
    dy = atm1->GetY() - atm2->GetY();
    dz = atm1->GetZ() - atm2->GetZ();

	dist2= dx*dx + dy*dy + dz*dz;
 	if(units == ANGSTROM_U)
		return(dist2);
	else
		return(dist2/(BOHR_TO_ANG*BOHR_TO_ANG));

}

double Vec3D::CalcAngle(const Vec3D* atm1, const Vec3D* atm2, const Vec3D* atm3 )
{
    register double ulen2,vlen2;
    register double ux,uy,uz;
    register double vx,vy,vz;
    register double temp;

	if(atm1 == NULL || atm2 == NULL || atm3 == NULL )
	{
		ErrorInMod(" Vec3D::CalcAngle() ",
			       " One of the atoms of the valence angle does not exist \n");
		return 0.0;
	}

    ux = atm1->GetX() - atm2->GetX();
    uy = atm1->GetY() - atm2->GetY();
    uz = atm1->GetZ() - atm2->GetZ();
    if( !ux && !uy && !uz )
        return( 0.0 );
    ulen2 = ux*ux + uy*uy + uz*uz;

    vx = atm3->GetX() - atm2->GetX();
    vy = atm3->GetY() - atm2->GetY();
    vz = atm3->GetZ() - atm2->GetZ();
    if( !vx && !vy && !vz )
        return( 0.0 );
    vlen2 = vx*vx + vy*vy + vz*vz;

    temp = (ux*vx + uy*vy + uz*vz)/sqrt(ulen2*vlen2);
    return( acos(temp) );
}


double Vec3D::CalcTorsion(const Vec3D* atm1, const Vec3D* atm2,
	   				   const Vec3D* atm3, const Vec3D* atm4 )
{
    double ax, ay, az;
    double bx, by, bz;
    double cx, cy, cz;
    double c12,c13,c23;
    double s12,s23;

    double cossq,sgn,om;
    double cosom,sinom;
    double len;

	if(atm1 == NULL || atm2 == NULL || atm3 == NULL || atm4 == NULL)
	{
		ErrorInMod(" Vec3D::CalcTorsion() ",
			       " One of the atoms of the torsion angle does not exist \n");
		return 0.0;
	}

    ax = atm2->GetX() - atm1->GetX();
    ay = atm2->GetY() - atm1->GetY();
    az = atm2->GetZ() - atm1->GetZ();
    if( !ax && !ay && !az )
        return( 0.0 );

    bx = atm3->GetX() - atm2->GetX();
    by = atm3->GetY() - atm2->GetY();
    bz = atm3->GetZ() - atm2->GetZ();
    if( !bx && !by && !bz )
        return( 0.0 );

    cx = atm4->GetX() - atm3->GetX();
    cy = atm4->GetY() - atm3->GetY();
    cz = atm4->GetZ() - atm3->GetZ();
    if( !cx && !cy && !cz )
        return( 0.0 );

    az = -az;  bz = -bz;  cz = -cz;

    len = sqrt(ax*ax + ay*ay + az*az);
    ax /= len;  ay /= len;  az /= len;
    len = sqrt(bx*bx + by*by + bz*bz);
    bx /= len;  by /= len;  bz /= len;
    len = sqrt(cx*cx + cy*cy + cz*cz);
    cx /= len;  cy /= len;  cz /= len;

    c12 = ax*bx + ay*by + az*bz;
    c13 = ax*cx + ay*cy + az*cz;
    c23 = bx*cx + by*cy + bz*cz;

    s12 = sqrt(1.0-c12*c12);
    s23 = sqrt(1.0-c23*c23);

    cosom = (c12*c23-c13)/(s12*s23);
    cossq = cosom*cosom;

    if( cossq >= 1.0 )
    {
		if( cosom < 0.0 )
        {
			return( PI );
        }
		else
			return( 0.0 );
    }

    sinom = sqrt(1.0-cossq);
    om = atan2(sinom,cosom);

    sgn =  ax*((by*cz)-(bz*cy));
    sgn += ay*((bz*cx)-(bx*cz));
    sgn += az*((bx*cy)-(by*cx));

    return( (sgn<0)? om : -om );
}

double Vec3D::CalcDihedral(const Vec3D* atm1,const Vec3D* atm2,const Vec3D* atm3,const Vec3D* atm4)
{
    return( PI - CalcTorsion(atm1,atm2,atm3,atm4) );
}

int Vec3D::VecProduct(Vec3D& vprod, const Vec3D& vec1, const Vec3D& vec2)
{
      vprod.SetX( vec1.GetY() * vec2.GetZ() - vec1.GetZ() * vec2.GetY());
	  vprod.SetY( vec1.GetZ() * vec2.GetX() - vec1.GetX() * vec2.GetZ());
	  vprod.SetZ( vec1.GetX() * vec2.GetY() - vec1.GetY() * vec2.GetX());
	  return TRUE;
}

double Vec3D::DotProduct(const Vec3D& vec1, const Vec3D& vec2)
{
	double dot = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];

	return dot;
}

int Vec3D::diff(Vec3D& c, const Vec3D& vec1, const Vec3D& vec2)
{
	c[0] = vec1[0] - vec2[0];
	c[1] = vec1[1] - vec2[1];
	c[2] = vec1[2] - vec2[2];

	return TRUE;
}

int Vec3D::sum(Vec3D& c, const Vec3D& vec1, const Vec3D& vec2)
{
	c[0] = vec1[0] + vec2[0];
	c[1] = vec1[1] + vec2[1];
	c[2] = vec1[2] + vec2[2];

	return TRUE;
}

void Vec3D::SetZeros()
{
	this->pos[0] = 0.0;
	this->pos[1] = 0.0;
	this->pos[2] = 0.0;
}


int Vec3D::SetAtomPos( Vec3D* pptr, const Vec3D* pptr1, const Vec3D* pptr2, const Vec3D* pptr3,
					   double dist, double val_angle, double dih_angle)
{
	  if( pptr == NULL || pptr1 == NULL || pptr2 == NULL || pptr3 == NULL)
	  {
		  ErrorInMod("Vec3D::SetAtomPos()", "Some of the points are invalid");
		  return FALSE;
	  }
	  int i;
	  Vec3D u1, u2, u3, u4;

	  for(i=0; i < 3; i++)
	  {
		u1[i]= (*pptr1)[i] - (*pptr2)[i];
		u2[i]= (*pptr3)[i] - (*pptr2)[i];
	  }
	  double len1 = u1.length();
	  double len2 = u2.length();

	  if(len1 < DBL_EPSILON )
	  {
		  ErrorInMod("Vec3D::SetAtomPos()", "points 1 and 2 coincide");
		  return FALSE;
	  }

	  if(len2 < DBL_EPSILON )
	  {
		  ErrorInMod("Vec3D::SetAtomPos()", "points 2 and 3 coincide");
		  return FALSE;
	  }

	  u1.Scale(1.0/len1);
	  u2.Scale(1.0/len2);
		
	  VecProduct(u3,u1,u2);
	  double len3 = u3.length();
	  if(len3 < DBL_EPSILON )
	  {
		  ErrorInMod("Vec3D::SetAtomPos", "points 1 and 2 and 4 are on one line ");
		  return FALSE;
	  }
	  u3.Scale(1.0/len3);

	  VecProduct(u4,u3,u1);
	  double len4 = u4.length();
	  u4.Scale(1.0/len4);

      double cosa=cos(val_angle);
      double sina=sin(val_angle);
      double cosb=cos(dih_angle);
      double sinb=sin(dih_angle);

      for(i = 0; i < 3; i++)
        (*pptr)[i]=(*pptr1)[i]+dist*(-u1[i]*cosa+u4[i]*sina*cosb + u3[i]*sina*sinb);

      return TRUE;
}


Quaternion::Quaternion()
{
	W = 1.0;
	X = 0.0;
	Y = 0.0;
	Z = 0.0;
}

Quaternion::Quaternion(const double w, const double x, const double y, const double z)
{
	W = w;
	X = x;
	Y = y;
	Z = z;
}

Quaternion::~Quaternion()
{
	
}


Quaternion operator * (const Quaternion &a, const Quaternion &b)
{
	double w,x,y,z;
	
	w = a.W*b.W - (a.X*b.X + a.Y*b.Y + a.Z*b.Z);
	
	x = a.W*b.X + b.W*a.X + a.Y*b.Z - a.Z*b.Y;
	y = a.W*b.Y + b.W*a.Y + a.Z*b.X - a.X*b.Z;
	z = a.W*b.Z + b.W*a.Z + a.X*b.Y - a.Y*b.X;
	
	return Quaternion(w,x,y,z); 
}

const Quaternion& Quaternion::operator *= (const Quaternion &q)
{
	double w = W*q.W - (X*q.X + Y*q.Y + Z*q.Z);
	
	double x = W*q.X + q.W*X + Y*q.Z - Z*q.Y;
	double y = W*q.Y + q.W*Y + Z*q.X - X*q.Z;
	double z = W*q.Z + q.W*Z + X*q.Y - Y*q.X;
	
	W = w;
	X = x;
	Y = y;
	Z = z;
	return *this;
}

const Quaternion& Quaternion::operator ~ ()
{
	X = -X;
	Y = -Y;
	Z = -Z;
	return *this;
}

const Quaternion& Quaternion::operator - ()
{
	double norme = sqrt(W*W + X*X + Y*Y + Z*Z);
	if (norme == 0.0)
		norme = 1.0;
	
	double recip = 1.0 / norme;
	
	W =  W * recip;
	X = -X * recip;
	Y = -Y * recip;
	Z = -Z * recip;
	
	return *this;
}

const Quaternion& Quaternion::Normalize()
{
	double norme = sqrt(W*W + X*X + Y*Y + Z*Z);
	if (norme == 0.0)
    {
		W = 1.0; 
		X = Y = Z = 0.0;
	}
	else
    {
		double recip = 1.0/norme;
		
		W *= recip;
		X *= recip;
		Y *= recip;
		Z *= recip;
	}
	return *this;
}

int Quaternion::GetQuaternion(double& w, double& x, double& y, double& z)
{
	w= W;
	x= X;
	y= Y;
	z= Z;
	return TRUE;
}

int Quaternion::SetQuaternion(double& w, double& x, double& y, double& z)
{
	W= w;
	X= x;
	Y= y;
	Z= z;
	return TRUE;
}


const Quaternion& Quaternion::QuaternionFromAxis(const double Angle, double x, double y, double z)
{
	double omega, s, c;
    
	s = sqrt(x*x + y*y + z*z);
	
	if (fabs(s) > FLT_EPSILON)
    {
		c = 1.0/s;
		
		x *= c;
		y *= c;
		z *= c;
		
		omega = -0.5f * Angle;
		s = (double)sin(omega);
		
		X = s*x;
		Y = s*y;
		Z = s*z;
		W = (double)cos(omega);
	}
	else
    {
		X = Y = 0.0f;
		Z = 0.0f;
		W = 1.0f;
	}
	Normalize();
	return *this;
}


int Quaternion::QuaternionToRotMat(const Quaternion& q, HaMat_double& rmat)
{
	rmat(1,1) = 1.0 - 2*q.Y*q.Y - 2*q.Z*q.Z;
	rmat(1,2) = 2*q.X*q.Y + 2*q.W*q.Z;
	rmat(1,3) = 2*q.X*q.Z - 2*q.W*q.Y;
	
	rmat(2,1) = 2*q.X*q.Y - 2*q.W*q.Z;
	rmat(2,2) = 1.0 - 2*q.X*q.X - 2*q.Z*q.Z;
	rmat(2,3) = 2*q.Y*q.Z + 2*q.W*q.X;
	
	rmat(3,1) = 2*q.X*q.Z + 2*q.W*q.Y;
	rmat(3,2) = 2*q.Y*q.Z - 2*q.W*q.X;
	rmat(3,3) = 1.0 - 2*q.X*q.X - 2*q.Y*q.Y;
//	   	PrintLog("Q2R %2.3f %2.3f %2.3f\n",rmat(1,1), rmat(1,2), rmat(1,3));
//		PrintLog("Q2R rot_mat %2.3f %2.3f %2.3f\n",rmat(2,1), rmat(2,2), rmat(2,3));
// 		PrintLog("Q2R rot_mat %2.3f %2.3f %2.3f\n",rmat(3,1), rmat(3,2), rmat(3,3));

	return TRUE;
}

int Quaternion::RotMatToQuaternion(const HaMat_double& rmat, Quaternion& q)
{
	double	t = rmat(1,1) + rmat(2,2) + rmat(3,3);
//	PrintLog("T %2.3f\n", t);
	double s, s1;
	
	
	if (t > 0.0f)
	{
		s = sqrt(t+1.0) ;
		q.W = 0.5 * s;
		s1 = 0.5/ s;
		q.X=  (rmat(2,3) - rmat(3,2)) * s1;
		q.Y= (rmat(3,1) - rmat(1,3)) * s1;
		q.Z= (rmat(1,2) - rmat(2,1)) * s1;
		q.Normalize();	
//		PrintLog("1!!QQQ\n");
//		q.PrintOn();
//		PrintLog("11!!QQQ\n");
		return TRUE;
	}
	else
	{
		if (rmat(2,2) > rmat(1,1) )
		{
			if (rmat(3,3) > rmat(2,2))
			{
				s = sqrt(1.0 + rmat(3,3) - rmat(1,1) - rmat(2,2));
				q.Z= 0.5 * s;	
				s1 = 0.5/s;
				q.X= (rmat(1,3) + rmat(3,1)) * s1;
				q.Y= (rmat(3,2) + rmat(2,3)) * s1;
				
				q.W= (rmat(1,2) - rmat(2,1) ) * s1;
				q.Normalize();
//				PrintLog("4!!QQQ\n");
				return TRUE;
			}
			else
			{
				s = sqrt(1.0 + rmat(2,2) - rmat(1,1) - rmat(3,3));
				q.Y= 0.5 * s;
				s1 = 0.5 /s;
				q.X= (rmat(2,1) + rmat(1,2)) * s1;
				q.Z= (rmat(3,2) + rmat(2,3)) * s1;
				q.W= (rmat(3,1) - rmat(1,3)) * s1;
				q.Normalize();
//				PrintLog("3!!QQQ\n");
				return TRUE;
			}
			
		}
		else
		{
			if( rmat(3,3) > rmat(1,1) )
			{
				s = sqrt(1.0 + rmat(3,3) - rmat(1,1) - rmat(2,2));
				q.Z= 0.5 * s;	
				s1 = 0.5/s;
				q.X= (rmat(1,3) + rmat(3,1)) * s1;
				q.Y= (rmat(3,2) + rmat(2,3)) * s1;
				
				q.W= (rmat(1,2) - rmat(2,1) ) * s1;
				q.Normalize();
//				PrintLog("4!!QQQ\n");
				return TRUE;
			}
			else
			{
				s = sqrt(1.0 + rmat(1,1) - rmat(2,2) - rmat(3,3));
				s1 = 0.5/s;
				q.X= 0.5 * s;
				q.Y= (rmat(2,1) + rmat(1,2)) * s1;
				q.Z= (rmat(1,3) + rmat(3,1)) * s1;
				q.W= (rmat(2,3) - rmat(3,2)) * s1;
				q.Normalize();	
//				PrintLog("2!!QQQ\n");
				return TRUE;
			}
		}
		
	}
}
	

void Quaternion::Slerp(const Quaternion &a,const Quaternion &b, const double t)
{
	double omega, cosom, sinom, sclp, sclq;
	
	cosom = a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W;
	
	
	if ((1.0f+cosom) > FLT_EPSILON)
	{
		if ((1.0f-cosom) > FLT_EPSILON)
		{
			omega = acos(cosom);
			sinom = sin(omega);
			sclp = sin((1.0f-t)*omega) / sinom;
			sclq = sin(t*omega) / sinom;
		}
		else
		{
			sclp = 1.0f - t;
			sclq = t;
		}
		
		X = sclp*a.X + sclq*b.X;
		Y = sclp*a.Y + sclq*b.Y;
		Z = sclp*a.Z + sclq*b.Z;
		W = sclp*a.W + sclq*b.W;
		
	}
	else
	{
		X =-a.Y;
		Y = a.X;
		Z =-a.W;
		W = a.Z;
		
		sclp = sin((1.0f-t) * PI * 0.5);
		sclq = sin(t * PI * 0.5);
		
		X = sclp*a.X + sclq*b.X;
		Y = sclp*a.Y + sclq*b.Y;
		Z = sclp*a.Z + sclq*b.Z;
		
	}
}

const Quaternion& Quaternion::exp()
{                               
	double Mul;
	double Length = sqrt(X*X + Y*Y + Z*Z);
	
	if (Length > 1.0e-4)
		Mul = sin(Length)/Length;
	else
		Mul = 1.0;
	
	W = cos(Length);
	
	X *= Mul;
	Y *= Mul;
	Z *= Mul; 
	
	return *this;
}

const Quaternion& Quaternion::log()
{
	double Length;
	
	Length = sqrt(X*X + Y*Y + Z*Z);
	Length = atan(Length/W);
	
	W = 0.0;
	
	X *= Length;
	Y *= Length;
	Z *= Length;
	
	return *this;
}

void Quaternion::PrintOn()
{
	PrintLog("Quaternion W= %2.3f X= %2.3f Y= %2.3f Z= %2.3f Length= %2.3f\n", W, X, Y, Z, sqrt(W*W+ X*X + Y*Y+ Z*Z));
}



int 
Rot3D::RotMatToEuler(const HaMat_double& rmat, double& phi, double& cos_theta, double& psi)
{
	phi = 0.0;
	cos_theta = 1.0;
	psi = 0.0;

	if( fabs(rmat(3,3)) > 1.0 )
	{
		ErrorInMod("rot_mat_to_euler()", " abs(rmat(3,3)) > 1 "); 
		return FALSE;
	}

	cos_theta = rmat(3,3);

	double sin_theta = sqrt(1.0 - cos_theta* cos_theta);
	double sin_psi,cos_psi,sin_phi,cos_phi;
	
	if( sin_theta > 1.0E-8)
	{
		sin_psi = rmat(1,3)/sin_theta;
		cos_psi = rmat(2,3)/sin_theta;
		sin_phi = rmat(3,1)/sin_theta;
		cos_phi = -rmat(3,2)/sin_theta;
	}
	else
	{
		cos_psi =  rmat(1,1);
		sin_psi = -rmat(2,1);
		sin_phi = 0.0;
		cos_phi = 1.0;
	}

	if( fabs(sin_psi) > 1.0 ||  fabs(cos_psi) > 1.0 || 
		fabs(sin_phi) > 1.0 ||  fabs(cos_phi) > 1.0 )
	{
		ErrorInMod("rot_mat_to_euler()", 
		"abs(sin(psi)) > 1 or abs(cos(psi)) > 1 or abs(sin(phi)) > 1 or abs(cos(phi)) > 1");
		return FALSE;
	}

	psi = acos(cos_psi);
	if(sin_psi < 0.0) psi = -psi;

	phi = acos(cos_phi);
	if(sin_phi < 0.0) phi = -phi;

	return TRUE;
}

int 
Rot3D::EulerToRotMat(const double phi, const double cos_theta, const double psi, 
		                        HaMat_double& rmat)
{
	if( fabs(cos_theta) > 1.0 )
	{
		ErrorInMod("euler_to_rot_mat()",
			       " abs( cos_theta) > 1.0 ");
		return FALSE;
	}
	double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	double sin_psi = sin(psi);
	double cos_psi = cos(psi);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);

	rmat(1,1) = cos_phi * cos_psi - sin_phi * cos_theta * sin_psi;
	rmat(1,2) = sin_phi * cos_psi + cos_phi * cos_theta * sin_psi;
	rmat(1,3) = sin_theta * sin_psi;

	rmat(2,1) = -cos_phi * sin_psi - sin_phi * cos_theta * cos_psi;
	rmat(2,2) = -sin_phi * sin_psi + cos_phi * cos_theta * cos_psi;
	rmat(2,3) = sin_theta * cos_psi;

	rmat(3,1) = sin_phi * sin_theta;
	rmat(3,2) = -cos_phi * sin_theta;
    rmat(3,3) = cos_theta;

	return TRUE;
}


int
Rot3D::IncrEulerAng( double& phi, double& cos_theta, double& psi, 
					      double delt_phi, double delt_cos_theta, double delt_psi)
{
    int ires = NormalizeEulerAng(phi, cos_theta, psi); 
    if( !ires) return FALSE;

	phi += delt_phi;
	if( phi > PI  ) phi = phi - 2*PI;
	if( phi < -PI ) phi = phi + 2*PI;

	psi += delt_psi;
	if( psi > PI  ) psi = psi - 2*PI;
	if( psi < -PI ) psi = psi + 2*PI;

	cos_theta += delt_cos_theta;

	int change_psi = 0;

	if( cos_theta > 1.0) 
	{
		cos_theta = 2.0 - cos_theta;
		change_psi = 1;
	}
	if( cos_theta < -1.0) 
	{
		cos_theta = -2.0 - cos_theta;
		change_psi = 1;
	}

	if( change_psi )
	{
//		double psi_old = psi;

//		psi = phi;
//		phi = psi_old;	

		if( psi > 0 ) 
		{
//			psi = PI - psi;
			psi = psi - PI;
		}
		else
		{
//			psi = -PI - psi;
			psi = psi + PI;
		}

		if( phi > 0 ) 
		{
//			phi = PI - phi;
			phi = phi - PI;
		}
		else
		{
//			phi = -PI - phi;
			phi = phi + PI;
		}

	}
	return TRUE;
}


int
Rot3D::NormalizeEulerAng( double& phi, double& cos_theta, double& psi)
{
	if( phi > PI || phi < -PI) 
	{
		int n_2pi = (int) (( phi + PI)/(2*PI));
		phi = phi - 2*PI*n_2pi;
	}

	if( psi > PI || psi < -PI) 
	{
		int n_2pi = (int) (( psi + PI)/(2*PI));
		psi = psi - 2*PI*n_2pi;
	}
	
	if( cos_theta > 1.0 || cos_theta < -1.0 )
	{
		PrintLog("Error in Rot3D::NormalizeEulerAng() \n");
		PrintLog(" Cos(Theta) = %9.3f is out of limits (-1,1) \n", cos_theta);
		return FALSE;
	}

	return TRUE;
}

int
Rot3D::RotMatToQuat(const HaMat_double& rmat, Quaternion& q)
{
	q.RotMatToQuaternion(rmat,q);
	return TRUE;
}

int 
Rot3D::QuatToRotMat(const Quaternion& q, HaMat_double& rmat)
{
	q.QuaternionToRotMat(q,rmat);
	return TRUE;
}


int 
Vec3D::CalcHlxParams(const Vec3D& c0_1, const Vec3D& v1_1, const Vec3D& v2_1, const Vec3D& v3_1,
                     const Vec3D& c0_2, const Vec3D& v1_2, const Vec3D& v2_2, const Vec3D& v3_2,
					 double& shift, double& slide, double& rise,
					 double& tilt, double& roll,double& twist,
					 int idir)
{
	double dir = 1.0;
	if(idir == -1)dir = -1.0;

	Vec3D n;    // Average normal
	Vec3D q;    // Average reference point
	Vec3D d;    // Average v1
	Vec3D f;    // base pair stack averages vector to major groove (f)
	Vec3D v,pl,pu;
				
	Vec3D::sum(n,v3_1,v3_2);
	n.normalize();
	n.Scale(dir);
				
	Vec3D::sum(q,c0_1,c0_2);
	q.Scale(0.5);
		
	Vec3D::sum(d,v1_1,v1_2);
				
	double dot = Vec3D::DotProduct(n,d);
		
	int j;

	for(j=0; j < 3; j++)
		d[j] -= dot* n[j]; // subtract component parallel to n
				
	d.normalize();
				
	Vec3D::VecProduct(f,n,d);
	Vec3D::diff(v,q,c0_1);
				
	double dot1 = Vec3D::DotProduct(v,n);
	double dot2 = dir * Vec3D::DotProduct(v3_1,n);
	double dl = dot1/dot2;
				
	for(j=0; j < 3; j++)
		pl[j] = c0_1[j] + dir*v3_1[j]*dl;
				
	Vec3D::diff(v,c0_2,q);
				
	dot1 = Vec3D::DotProduct(v,n);
	dot2 = dir * Vec3D::DotProduct(v3_2,n);
	double du = dot1/dot2;
				
	for(j=0; j < 3; j++)
		pu[j] = c0_2[j] - dir*v3_2[j]*du;
				
	rise = (dl + du);
	Vec3D::diff(v,pu,pl);
				
	shift = Vec3D::DotProduct(d,v);
	slide = Vec3D::DotProduct(f,v);
				
	Vec3D t,r;
	Vec3D::VecProduct(t,v3_2,d);
	t.Scale(dir);
	double tlen = t.length();
	dot = Vec3D::DotProduct(f,t)/tlen;
				
	if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
				
	double cln = RAD_TO_DEG* acos(dot);
	Vec3D::VecProduct(v,f,t);
	if(Vec3D::DotProduct(v,d) < 0.0) cln = -cln;
	Vec3D::VecProduct(r,d,t);
	double rlen = r.length();
	dot = dir* Vec3D::DotProduct(v3_2,r)/rlen;
				
	if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
				
	double tip = RAD_TO_DEG* acos(dot);
	Vec3D::VecProduct(v,r,v3_2);
	v.Scale(dir);
	if( Vec3D::DotProduct(v,t) < 0.0) tip = -tip;
				
	tilt = 2*cln;
	roll = 2*tip;
	twist = 0.0;
				
	int ii;
	for(ii=0; ii < 2; ii++)
	{
		double sa = sin(-DEG_TO_RAD*cln);
		double ca = cos(DEG_TO_RAD*cln);
		if(ii == 1) sa = sin(DEG_TO_RAD*cln);
					
		Vec3D fp = f;
		fp.Rotate(d,ca,sa);
		Vec3D vh,vn;
		if(ii == 0 ) { vh = v2_1; vn = v3_1; }
        if(ii == 1 ) { vh = v2_2; vn = v3_2; } 
		dot = -1.0 * Vec3D::DotProduct(fp,vh);
		if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
		double wdg = RAD_TO_DEG* acos(dot);
					
		Vec3D::VecProduct(v,fp,vh);
		v.Scale(-1.0);
		dot = dir* Vec3D::DotProduct(v,vn);
		if( ii == 0 && dot > 0.0) wdg = -wdg;
		if( ii == 1 && dot < 0.0) wdg = -wdg;
		twist += wdg;
	}
				
	twist -= 360.0*floor( twist/360.0);
	if(fabs(twist) > 180.0) twist = (twist > 0)? twist - 360.0: twist + 360.0;

	return TRUE;
}

extern "C"
{
	void fitstr_(int* nat,double* x1, double* x2, double* eps,
		int* irot, double* r,double* cc,double* w,int* ia);
}


double PointContainer::GetSuperimposeMat( HaVec_double& ref_crd, PointContainer& g2,
                                                  HaMat_double& rot_mat, HaVec_double& transl_vec )
{
	double eps;
	Vec3DValArray ref_crd_vec3d;
	int npt = g2.GetNumPt();
	ref_crd_vec3d.resize(npt);
	ref_crd_vec3d.SetCrdFromArray(ref_crd);
	int ires = PointContainer::GetSuperimposeMat(ref_crd_vec3d, g2, rot_mat, transl_vec, eps );
	if( ires ) return eps;
	return -1.0;
}

int PointContainer::GetSuperimposeMat( PointContainer& grp1,  PointContainer& grp2, 
		                               HaMat_double& rot_mat,  HaVec_double& transl_vec, double& eps)
{
//! calls fortran
//!	subroutine fitstr(nat,xfmol,xsmol,eps,irot,rot_mat,trans_mat,work1,work2)
//! nat - number of atoms
//! xfmol(3*nat), xsmol(3*nat) - atomic coord. arrays (x,y,z etc)
//! work1(2*NR)- working array , work2(NR+3*NR+3*NR) - working array
//! eps = 0.5*W*(X-XP)**2  - in matched orientation
//! irot = 1 - in a successful run, 0 - in a failure
//! rot_mat(3,3) - rotation matrix
//! trans_mat(3) -translational vector
//! such that to restore xfmol:
//! xfmol(i,j) ~= \sum_k R(i,k) xsmol(k) +cc(i)

//	PrintLog(" PointContainer::GetSuperimposeMat() pt 1 \n");

	rot_mat.newsize(3,3);
	rot_mat = 0.0;
	rot_mat(1,1) = 1.0; rot_mat(2,2) = 1.0; rot_mat(3,3) = 1.0; 
	transl_vec.newsize(3);
	transl_vec = 0.0;
	eps = 0.0;

	int n = grp1.GetNumPt();
	if( n != grp2.GetNumPt() )
	{
		ErrorInMod("PointContainer::GetSuperimposeMat()",
			       "Different number of points in two collections");
		return FALSE;
	}

	PointIterator* p_itr1 = grp1.GetPointIteratorPtr();
	PointIterator* p_itr2 = grp2.GetPointIteratorPtr();

	HaVec_double xfmol(3*n+6,0.0), xsmol(3*n+6,0.0); // fitstr_() can add up two extra axxiliary points therefore allocate spare for 3 more coordinates
    
	Vec3D* ppt = NULL;
	int i=0;
	for (ppt= p_itr1->GetFirstPt(); ppt; ppt= p_itr1->GetNextPt())
	{
		xfmol[3*i  ] = ppt->GetX();  //coord in Ang
		xfmol[3*i+1] = ppt->GetY();
		xfmol[3*i+2] = ppt->GetZ();
		i++;
	}
	delete p_itr1;

	i=0;
	for (ppt= p_itr2->GetFirstPt(); ppt; ppt= p_itr2->GetNextPt())
	{
		xsmol[3*i  ] = ppt->GetX();  //coord in Ang
		xsmol[3*i+1] = ppt->GetY();
		xsmol[3*i+2] = ppt->GetZ();
		i++;
	}
	delete p_itr2;

// create variables needed for fortran subroutine and call it	

	HaVec_double work1(7*n+14);
    HaVec_int iwork2(2*n+4);

//array.begin() - address of the first element
    int irot;
	int nat = n;
	fitstr_(&nat,xfmol.v(),xsmol.v(), &eps, &irot,
		rot_mat.v(),transl_vec.v(),work1.begin(),iwork2.begin());

	eps = sqrt(2.0*eps/nat);

	if(irot == 0 ) return FALSE;

    return TRUE;
}

double PointContainer::CalcRMSD(PointContainer& g1, PointContainer& g2, int transform )
{
	int nat = g1.GetNumPt();
	if( nat != g2.GetNumPt() || nat == 0 )
	{
		PrintLog( "Error in PointContainer::CalcRMSD() \n");
		PrintLog( "Different number of points in two collections \n");
		return -1.0;
	}

	double eps = -1.0;
	if( transform ) 
	{
		HaMat_double rot_mat(3,3);
		HaVec_double trans_vec(3);
		double eps;
		int ires = GetSuperimposeMat( g1, g2, rot_mat, trans_vec, eps);
		if(!ires) eps = -1.0;
	}
	else
	{
		std::auto_ptr<PointIterator> p_auto1( g1.GetPointIteratorPtr() );
		PointIterator* p_itr1 = p_auto1.get();
		std::auto_ptr<PointIterator> p_auto2( g2.GetPointIteratorPtr() );
		PointIterator* p_itr2 = p_auto2.get();

		Vec3D* ppt1 = p_itr1->GetFirstPt();
		Vec3D* ppt2 = p_itr2->GetFirstPt();

		eps = 0.0;
		for( ; ppt1 != NULL; ppt1 = p_itr1->GetNextPt(),ppt2 = p_itr2->GetNextPt() )
		{
			eps += Vec3D::CalcDistanceSq( ppt1, ppt2 );
		}
		eps = sqrt( eps/nat );
	}
	return eps;
}

int PointContainer::Transform( const HaMat_double& rot_mat, const HaVec_double& trans_vec)
{
	Vec3D* pp;
	std::auto_ptr<PointIterator> pitr(GetPointIteratorPtr());
	
	for (pp = pitr->GetFirstPt();pp; pp = pitr->GetNextPt() )
	{
		double x=  pp->GetX();
       	double y = pp->GetY();
		double z = pp->GetZ();
		
		double xnew = x*rot_mat(1,1)+y*rot_mat(1,2)+z*rot_mat(1,3)+trans_vec(1);
		double ynew = x*rot_mat(2,1)+y*rot_mat(2,2)+z*rot_mat(2,3)+trans_vec(2);
		double znew = x*rot_mat(3,1)+y*rot_mat(3,2)+z*rot_mat(3,3)+trans_vec(3);
		
		pp->SetX(xnew);
		pp->SetY(ynew);
		pp->SetZ(znew);
	}
	return TRUE;
}

int PointContainer::SaveCrdToArray(HaVec_double& crd_arr) 
{
	int npt = this->GetNumPt();
	crd_arr.resize(3*npt);
	
	Vec3D* pp;
	std::auto_ptr<PointIterator> pitr(GetPointIteratorPtr());
	
	int i = 0;
	for (pp = pitr->GetFirstPt();pp; pp = pitr->GetNextPt() )
	{
		crd_arr[3*i  ] = pp->GetX();
		crd_arr[3*i+1] = pp->GetY();
		crd_arr[3*i+2] = pp->GetZ();
		i++;
	}
	return TRUE;
}

HaVec_double PointContainer::GetCrdArray() 
{
	HaVec_double crd_arr;
	SaveCrdToArray(crd_arr);
	return crd_arr;
}

int PointContainer::SetCrdFromArray( const HaVec_double& crd_arr)
{
	int npt = this->GetNumPt();
	if( crd_arr.size() != npt*3 )
	{
		PrintLog(" Error in PointContainer::SetCrdFromArray() \n ");
		PrintLog(" Size of coordinate array = %d  doesn't correspond to number of points = %d \n ",
			       crd_arr.size(), npt );
		return FALSE;
	}
	
	Vec3D* pp;
	std::auto_ptr<PointIterator> pitr(GetPointIteratorPtr());
	int i = 0;
	for (pp = pitr->GetFirstPt();pp; pp = pitr->GetNextPt() )
	{
		pp->SetX( crd_arr[3*i  ] );
		pp->SetY( crd_arr[3*i+1] );
		pp->SetZ( crd_arr[3*i+2] );
		i++;
	}
	return TRUE;
}

int PointContainer::IsWithinRadius(Vec3D* pptr, double dlimit2 ) const
{
    double dx,dy,dz;
    double dist2;

	const Vec3D* pp;

	PointIteratorGen_const pitr(*this);
	
	for( pp = pitr.GetFirstPt(); pp ; pp = pitr.GetNextPt() )
	{    
		dx = pptr->GetX() - pp->GetX();
		dist2 = dx*dx;
		if( dist2 > dlimit2 ) continue;
		dy = pptr->GetY() - pp->GetY();
		dist2 += dy*dy;
		if( dist2 > dlimit2 ) continue;
		dz = pptr->GetZ() - pp->GetZ();
		dist2 += dz*dz;
		if( dist2 > dlimit2 ) continue;

		return( TRUE );
	}
    return( FALSE );
}


bool PointContainer::GetMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
					           double& MaxX_v, double& MaxY_v, double& MaxZ_v) const
{
	PointIteratorGen_const pitr(*this);
	const Vec3D* pptr;

	MinX_v = 10.0E6;
	MinY_v = 10.0E6;
	MinZ_v = 10.0E6;
	MaxX_v = -10.0E6;
	MaxY_v = -10.0E6;
	MaxZ_v = -10.0E6;

	int n_sel_at = 0;
	for(pptr = pitr.GetFirstPt(); pptr; pptr = pitr.GetNextPt())
	{
		double x = pptr->GetX();
		double y = pptr->GetY();
		double z = pptr->GetZ();
		
		if( x < MinX_v) MinX_v = x;
		if( y < MinY_v) MinY_v = y;
		if( z < MinZ_v) MinZ_v = z;
		if( x > MaxX_v) MaxX_v = x;
		if( y > MaxY_v) MaxY_v = y;
		if( z > MaxZ_v) MaxZ_v = z;
		
		n_sel_at++;
	}
	
	if(n_sel_at == 0)
	{
		MinX_v = 0.0;
		MinY_v = 0.0;
		MinZ_v = 0.0;
		MaxX_v = 0.0;
		MaxY_v = 0.0;
		MaxZ_v = 0.0;
		return false;
	}
	return true;
}

bool PointContainer::GetAverageCoord(double& avx, double& avy, double& avz) const
{
	avx = 0.0;
	avy = 0.0;
	avz = 0.0;

	PointIteratorGen_const pitr(*this);
	
	const Vec3D* pptr;
	int npt = 0;

	for(pptr= pitr.GetFirstPt(); pptr; pptr = pitr.GetNextPt() )
	{
		avx += pptr->GetX();
		avy += pptr->GetY();
		avz += pptr->GetZ();
		npt++;
	}

	if(npt == 0) return false;

	avx = avx/npt;
    avy = avy/npt;
    avz = avz/npt;

	return true;
}

Vec3D PointContainer::GetAverageCoord() const
{
	Vec3D avg_pt;
	GetAverageCoord( avg_pt[0],avg_pt[1],avg_pt[2]);
	return avg_pt;
}

int PointContainer::FindCoordMatch( PointContainer& g1, PointContainer& g2, PtrPtrMap& pt_pt_map)
{
	pt_pt_map.clear();

	if( g1.GetNumPt() == 0 || g2.GetNumPt() == 0) return FALSE;

	Vec3D pt_min_1, pt_max_1, pt_min_2, pt_max_2;

	g1.GetMinMaxCrd( pt_min_1[0], pt_min_1[1], pt_min_1[2],
				     pt_max_1[0], pt_max_1[1], pt_max_1[2]);

	g2.GetMinMaxCrd( pt_min_2[0], pt_min_2[1], pt_min_2[2],
				     pt_max_2[0], pt_max_2[1], pt_max_2[2]);


	int i;

	for(i=0; i < 3; i++)
	{
		if( pt_min_2[i] < pt_min_1[i]) pt_min_1[i] = pt_min_2[i];
		if( pt_max_2[i] > pt_max_1[i]) pt_max_1[i] = pt_max_2[i];
		pt_min_1[i] = pt_min_1[i] - 0.5;
		pt_max_1[i] = pt_max_1[i] + 0.5;
	}	

	BoxPartition part_1, part_2;

	part_1.SetDimensions( 10, 10, 10);
	part_2.SetDimensions( 10, 10, 10);

	part_1.SetBoundaries(pt_min_1[0],pt_min_1[1],pt_min_1[2],
		                 pt_max_1[0],pt_max_1[1],pt_max_1[2]);

	part_2.SetBoundaries(pt_min_1[0],pt_min_1[1],pt_min_1[2],
		                 pt_max_1[0],pt_max_1[1],pt_max_1[2]);

	part_1.DistributePointsToCells(g1);
	part_2.DistributePointsToCells(g2);

	int nc = part_1.size();

	double crit = 0.001;

	for( i = 0; i < nc; i++)
	{
        VecPtr& pts_1 = part_1[i];
		VecPtr& pts_2 = part_2[i];

		int n1 = pts_1.size();
		int n2 = pts_2.size();

		if( n1 == 0) continue;
		if( n2 == 0) continue;

		int j1,j2;

		for( j1=0; j1 < n1; j1++ )
		{
			Vec3D* pt1 = (Vec3D*) pts_1[j1];

			for( j2=0; j2 < n2; j2++)
			{
				Vec3D* pt2 = (Vec3D*) pts_2[j2];
				if( fabs( (*pt1)[0] - (*pt2)[0]) > crit) continue;
				if( fabs( (*pt1)[1] - (*pt2)[1]) > crit) continue;
				if( fabs( (*pt1)[2] - (*pt2)[2]) > crit) continue;

				pt_pt_map[pt1] = pt2;
				pt_pt_map[pt2] = pt1;
			}
		}
	}
	return TRUE;
}



Vec3DValArrayIterator::Vec3DValArrayIterator(Vec3DValArray* pt_array)
{
	cur_arr = pt_array;
	cur_idx = 0;
}

Vec3DValArrayIterator::~Vec3DValArrayIterator()
{

}

Vec3D*
Vec3DValArrayIterator::GetFirstPt()
{
	cur_idx = 0;
	if(cur_arr == NULL) return NULL;
	if(cur_arr->size() == 0) return NULL;
	return &((*cur_arr)[0]);
}

Vec3D*
Vec3DValArrayIterator::GetNextPt()
{
	if(cur_arr == NULL) return NULL;
	cur_idx++;
	if(cur_arr->size() == 0 || cur_idx >= cur_arr->size()) return NULL;
	return &((*cur_arr)[cur_idx]);
}

Vec3DValArrayIterator_const::Vec3DValArrayIterator_const(const Vec3DValArray* pt_array)
{
	cur_arr = pt_array;
	cur_idx = 0;
}

Vec3DValArrayIterator_const::~Vec3DValArrayIterator_const()
{

}

const Vec3D*
Vec3DValArrayIterator_const::GetFirstPt()
{
	cur_idx = 0;
	if(cur_arr == NULL) return NULL;
	if(cur_arr->size() == 0) return NULL;
	return &((*cur_arr)[0]);
}

const Vec3D*
Vec3DValArrayIterator_const::GetNextPt()
{
	if(cur_arr == NULL) return NULL;
	cur_idx++;
	if(cur_arr->size() == 0 || cur_idx >= cur_arr->size()) return NULL;
	return &((*cur_arr)[cur_idx]);
}

Vec3DValArray::Vec3DValArray()
{

}

Vec3DValArray::Vec3DValArray(int n): vector<Vec3D>(n)
{

}


Vec3DValArray::~Vec3DValArray()
{

}

PointIterator* 
Vec3DValArray::GetPointIteratorPtr() 
{ 
	return new Vec3DValArrayIterator(this); 
}

PointIterator_const* 
Vec3DValArray::GetPointIteratorPtr() const
{ 
	return new Vec3DValArrayIterator_const(this); 
}

int
Vec3DValArray::GetNumPt() const
{ 
	return size(); 
}


BoxPartition::BoxPartition()
{
	SetDimensions(21,21,21);
}

int
BoxPartition::SetDimensions(int nx_new, int ny_new, int nz_new)
{
	if( nx_new <= 0 || ny_new <= 0 || nz_new <= 0)
	{
		nx = 1;
		ny = 1;
		nz = 1;
		resize(1);
		clear();
		return FALSE;
	}

	nx = nx_new;
	ny = ny_new;
	nz = nz_new;

	clear();
	resize(nx*ny*nz);
	return TRUE;
}

int
BoxPartition::GetPointCellIdx(const Vec3D* ppt)
{
	return this->GetPointCellIdx( ppt->GetX(), ppt->GetY(), ppt->GetZ());
}

int
BoxPartition::GetPointCellIdx(double x, double y, double z)
{
	int ix = GetIdx_x(x);
	int iy = GetIdx_y(y);
	int iz = GetIdx_z(z);

	if( ix < 0 || iy < 0 || iz < 0) return -1;

	return ( ix + nx*iy + nx*ny*iz); 
}


int
BoxPartition::GetLinCellIdx( int ix, int iy, int iz)
{
	return ( ix + nx*iy + nx*ny*iz); 
}


int
BoxPartition::SetBoundaries(double xmin_new, double ymin_new, double zmin_new, 
		                    double xmax_new, double ymax_new, double zmax_new)
{
	xmin = xmin_new;
	ymin = ymin_new;
	zmin = zmin_new;

	dx = (xmax_new - xmin_new)/nx;
	dy = (ymax_new - ymin_new)/ny;
	dz = (zmax_new - zmin_new)/nz;
	
	return True;
}

int BoxPartition::DistributePointsToCells( PointContainer& pt_coll)
{
	PointIteratorGen pitr(pt_coll);

	Vec3D* pt;

	for( pt = pitr.GetFirstPt(); pt; pt = pitr.GetNextPt())
	{
		this->AddPoint(pt);
	}
	return TRUE;
}

int BoxPartition::AddPoint( Vec3D* pt)
{
	if( pt == NULL) return FALSE;
	int idx = GetPointCellIdx(pt->GetX(), pt->GetY(), pt->GetZ());
	if( idx < 0) return FALSE;
	(this->at(idx)).push_back(pt);
	return TRUE;
}

int 
BoxPartition::GetIdx_x( double x_coord ) 
{ 
	int ix = (int)(( x_coord - xmin )/ dx); 
	if( ix < 0 ) ix = -1; 
	if( ix >= nx  ) ix = -1 ;
	return ix;
}

int 
BoxPartition::GetIdx_y( double y_coord ) 
{ 
	int iy = (int)(( y_coord - ymin )/dy); 
	if( iy < 0 ) iy = -1; 
	if( iy >= ny) iy = -1;
	return iy;
}

int 
BoxPartition::GetIdx_z( double z_coord ) 
{ 
	int iz = (int) (( z_coord - zmin )/dz); 
	if( iz < 0 ) iz = -1; 
	if( iz >= nz) iz = -1 ;
	return iz;
}

int 
BoxPartition::SetRegionRad(double dist)
{
	if( fabs(dx) < DBL_EPSILON ||  fabs(dx) < DBL_EPSILON || fabs(dy) < DBL_EPSILON )
	{
		ErrorInMod("BoxPartition::SetRegionRad()",
			       " dx, dy or dz equal 0.0 ");
		ixrad = 0;
		iyrad = 0;
		izrad = 0;
		return False;

	}

	if( dist < DBL_EPSILON)
	{
		ixrad = 0;
		iyrad = 0;
		izrad = 0;
	}
	else
	{
		ixrad = (int) (dist/dx) +1;
		iyrad = (int) (dist/dy) +1;
		izrad = (int) (dist/dz) +1;
	}
	return true;
}

int
BoxPartition::GetNeighbors(Vec3D& pt, AtomGroup& neighbors)
{
	neighbors.clear();

	int nxd = GetIdx_x( pt.GetX());
	int nyd = GetIdx_y( pt.GetY());
	int nzd = GetIdx_z( pt.GetZ());

	int iax_min = MaxFun(nxd - ixrad, 0);
	int iay_min = MaxFun(nyd - iyrad, 0);
	int iaz_min = MaxFun(nzd - izrad, 0);
	int iax_max = MinFun(nxd + ixrad, nx - 1);
	int iay_max = MinFun(nyd + iyrad, ny - 1);
	int iaz_max = MinFun(nzd + izrad, nz - 1);

	int nn = 0;
	for(int iax = iax_min; iax <= iax_max; iax++)
	{
		for(int iay = iay_min; iay <= iay_max; iay++)
		{
			for(int iaz = iaz_min; iaz <= iaz_max; iaz++)
			{
				 int idx = GetLinCellIdx(iax,iay,iaz);
				 if( (*this)[idx].empty())
					 continue;
				 VecPtr::iterator aitr;
				 for(aitr = (*this)[idx].begin(); aitr != (*this)[idx].end(); aitr++)
				 {
					neighbors.push_back((HaAtom*)(*aitr));
					nn++;
				 }
			}
		}
	}
	return nn;
}

BoxRegionPointIterator
BoxPartition::GetRegionPointIterator()
{
	BoxRegionPointIterator iter;
	return iter;
}

BoxRegionPointIterator::BoxRegionPointIterator()
{
	iax_min = -1;
	iay_min = -1;
	iaz_min = -1;
	iax_max = -1;
	iay_max = -1;
	iaz_max = -1;

	partition = NULL;
}

BoxRegionPointIterator::BoxRegionPointIterator(const BoxRegionPointIterator& iter_ref)
{
	iax_min = iter_ref.iax_min;
	iay_min = iter_ref.iay_min;
	iaz_min = iter_ref.iaz_min;
	iax_max = iter_ref.iax_max;
	iay_max = iter_ref.iay_max;
	iaz_max = iter_ref.iax_max;

	partition = iter_ref.partition;
}

BoxRegionPointIterator::~BoxRegionPointIterator()
{

}

Vec3D*
BoxRegionPointIterator::GetFirstPt()
{
	return NULL;
}

Vec3D*
BoxRegionPointIterator::GetNextPt()
{
	return NULL;
}

int
BoxRegionPointIterator::GetNumPt()
{
	return 0;
}



