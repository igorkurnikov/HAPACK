/*TEX
%
% QCHEM - IPACK
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% 
\subsection{Implementation of \name{HalfSphere}} 
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "rec.h"
#include "tquadruple.h"
#include "four_stor.h"
#include "sphere.h"

ARRAY<Mat>              Sphere::TRANSFORM;
ARRAY<ARRAY<Momentum> > Sphere::SYMINFO;

void Sphere::transform_clear()
{
  int i;
  for(i = 0; i < TRANSFORM.size(); i++)
    TRANSFORM[i].reset(0,0);
  TRANSFORM.reset(0);

  for(i = 0; i < SYMINFO.size(); i++)
    SYMINFO[i].reset(0);
  SYMINFO.reset(0);
}

void transform_init();

HalfSphere::HalfSphere(RecursionStorage& rec, int dbg,Allocator& tmpl) 
     : alloc(tmpl)
{  
  if (Sphere::TRANSFORM.size() == 0)
    Sphere::transform_init();
  
  Momentum m;
  l1      = rec.l1;
  l2      = rec.l2;
  max_mu  = rec.max_mu;
  max_aux = rec.max_aux;
  debug   = dbg;
  no_l1 =  2*l1 +1;
  no_l2 =  m.number(l2);
  no_mu = rec.no_mu;

  fields.reset(no_l1*no_l2*(max_aux+1)*no_mu);
  fields.set(0);
 
  if (debug)
    cout << " +++ Entries in Half-Sphere: " << no_l1 <<  " " 
      << no_l2 << " " << no_l1*no_l2*(max_aux+1)*no_mu << endl;

  if (l1 >= Sphere::TRANSFORM.size())
    error(" ++++ Sphere -- Transformation Matrix not initialized. ");
  Mat& T = Sphere::TRANSFORM[l1];
  
  for(int l=0; l < no_l1; l++)
  {
    if (debug) 
      cout << " +++ Creating Function: " << l << endl;
    
    for(int aux = 0; aux <= max_aux; aux++)
    for(int mut = 0; mut <= max_mu; mut++)
    for(Momentum mu(mut); mu(); ++mu)
    for(Momentum m2(l2); m2(); ++m2)
    {
      Momentum m1(l1);
      IntegArray* source = rec.find(m1,m2,aux,mu);
      if (!source) break;
      
      IntegArray *array = new IntegArray(*source,alloc);
      array -> zero();
      
      for(; m1(); ++m1) // now loop over m1
      {
	source = rec.find(m1,m2,aux,mu);
	assert(source != 0);
	double cc = T(l,m1.local_index());
	if (fabs(cc) > 0)
	{
	  if (debug)
	    cout << "            Term: l1/m1: " << l << " " << m1 
	         << " C: " << cc << endl;
	  array -> add_to(cc,*source);
	}
      }
      
      int idx = compute_index(l,m2,aux,mu);
      if (debug)
	cout << " +++ Storing: " << l << " " << m2 << " Aux: " << aux << " " 
	   << mu << endl;
      fields[idx] = array;
    }
  }
}
/*TEX
\subsection{Implementation of \name{Sphere}}
*/  
Sphere::Sphere(RecursionStorage& rec, int dbg,Allocator& tmpl) : alloc(tmpl)
{  
  debug = dbg;
  
  if (TRANSFORM.size() == 0)
    transform_init();

  l1      = rec.l1;
  l2      = rec.l2;
  max_mu  = rec.max_mu;
  max_aux = rec.max_aux;
  debug   = dbg;
  
  no_l1 =  2*l1 +1;
  no_l2 =  2*l2 +1;
  no_mu =  rec.no_mu;  

  fields.reset(no_l1*no_l2*(max_aux+1)*no_mu);
  fields.set(0);
 
  if (debug)
    cout << " +++ Entries in Sphere: " << no_l1 <<  " " 
      << no_l2 << " " << no_l1*no_l2*(max_aux+1)*no_mu << endl;

  HalfSphere half(rec,dbg,tmpl);
    
  if (l2 >= TRANSFORM.size())
    error(" ++++ Sphere -- Transformation Matrix not initialized. ");
  Mat& T = TRANSFORM[l2];

  for(int l2new=0; l2new < no_l2; l2new++)
  {
    if (debug) 
      cout << " +++ Sphere: Creating Function: " << l2new << endl;
    
    for(int l1new = 0 ; l1new < no_l1; l1new++)
    for(int aux = 0; aux <= max_aux; aux++)
    for(int mut = 0; mut <= max_mu; mut++)
    for(Momentum mu(mut); mu(); ++mu)
    {
      Momentum m2(l2);
      IntegArray* source = half.find(l1new,m2,aux,mu);
      if (!source) break;
       
      IntegArray *array = new IntegArray(*source,alloc);
      array -> zero();
      
      for(; m2(); ++m2) // now loop over m1
      {
	source = half.find(l1new,m2,aux,mu);
	assert(source != 0);
	double cc = T(l2new,m2.local_index());
	if (fabs(cc) > 0)  
	{
	  if (debug)
	    cout << "            Term: l2/m2: " << l2new 
	         << " " << m2 << " C: " << cc << endl;
	  array -> add_to(cc,*source);
	}
      }
      int idx = compute_index(l1new,l2new,aux,mu);
      if (debug) 
	cout << " +++ Storing: " << l1new << " " << l2new 
	  << " Aux: " << aux << " " << mu << endl;
      fields[idx] = array;
    }
  }
}

HalfSphere::~HalfSphere() 
{
  for(int i=0; i< fields.size(); i++)
    if (fields[i]) delete fields[i];
  fields.reset(0);
}


/*TEX
\subsection{Transformation Coefficients}

The following function initializes the transformation coeffcients for up to g-waves. 
*/ 

void Sphere::transform_init()
{
  TRANSFORM.reset(6);
  SYMINFO.reset(6);
  
  Mat& ms = TRANSFORM[0];
  ms.reset(1,1);
  ms(0,0) = 1.0;

  ARRAY<Momentum>& ss = Sphere::SYMINFO[0];
  ss.reset(1);
  ss[0]   = Momentum(0,0,0);
  

  Mat& mp = TRANSFORM[1];
  mp.reset(3,3);
  mp.set(0);
  mp(0,0) = 1.0;
  mp(1,1) = 1.0;
  mp(2,2) = 1.0;

  ARRAY<Momentum>& sp = SYMINFO[1];
  sp.reset(3);
  for(Momentum p(1); p(); ++p)
    sp[p.local_index()] = p;
  
  // Spherical D-WAVE functions
  //   xy, xz, yz, 2 zz - xx - yy, xx - yy

  Mat& md = TRANSFORM[2];
  md.reset(5,6);
  md.set(0);

  double ff = 1/sqrt(3.0);
  
//  md(0, Momentum(1,1,0).local_index()) =  1;  // xy
//  md(1, Momentum(1,0,1).local_index()) =  1;  // xz
//  md(2, Momentum(0,1,1).local_index()) =  1;  // yz

//  md(3, Momentum(0,0,2).local_index()) =  2;  // 2.0* zz - xx - yy
//  md(3, Momentum(2,0,0).local_index()) = -1;
//  md(3, Momentum(0,2,0).local_index()) = -1;

//  md(4, Momentum(2,0,0).local_index()) =  1;  // xx - yy
//  md(4, Momentum(0,2,0).local_index()) = -1;

// Igor Kurnikov modifications to coincide with Gaussian sequence of spherical orbitals

  md(0, Momentum(0,0,2).local_index()) =  2; // 2.0* zz - xx - yy 
  md(0, Momentum(2,0,0).local_index()) = -1;
  md(0, Momentum(0,2,0).local_index()) = -1;

  md(1, Momentum(1,0,1).local_index()) =  1;  // xz
  md(2, Momentum(0,1,1).local_index()) =  1;  // yz

  md(3, Momentum(2,0,0).local_index()) =  1;  // xx - yy
  md(3, Momentum(0,2,0).local_index()) = -1;
 
  md(4, Momentum(1,1,0).local_index()) =  1;   // xy
	  
  ARRAY<Momentum>& sd = SYMINFO[2];
  sd.reset(md.size1());
  int l;
  for(l = 0; l < md.size1(); l++)
    for(Momentum m(2); m(); ++m)
      if (fabs(md(l,m.local_index())) > 0)
      {
	sd[l] = m;
	break;
      }
	  
  // F-WAVE Functions

  Mat& mf = TRANSFORM[3];
  mf.reset(7,10);
  mf.set(0);

  // new function 0: xyz
    
  mf(0, Momentum(1,1,1).local_index()) = 1;
  
  // new function 1: fx = 2 xxx - 3 xyy - 3 xzz     
  
  mf(1, Momentum(3,0,0).local_index()) =  2;
  mf(1, Momentum(1,2,0).local_index()) = -3;
  mf(1, Momentum(1,0,2).local_index()) = -3;
      
  // new function 2: fy = 2 yyy - 3 xxy - 3 yzz     

  mf(2, Momentum(0,3,0).local_index()) =  2;
  mf(2, Momentum(2,1,0).local_index()) = -3;
  mf(2, Momentum(0,1,2).local_index()) = -3;
  
  // new function 3: fz = 2 zzz - 3 xxz - 3 yyz     

  mf(3, Momentum(0,0,3).local_index()) =  2;
  mf(3, Momentum(2,0,1).local_index()) = -3;
  mf(3, Momentum(0,2,1).local_index()) = -3;

  // new function 4: f'x = xyy - xzz
    
  mf(4, Momentum(1,2,0).local_index()) =  1;
  mf(4, Momentum(1,0,2).local_index()) = -1;

  // new function 5: f'y = yzz - xxy

  mf(5, Momentum(0,1,2).local_index()) =  1;
  mf(5, Momentum(2,1,0).local_index()) = -1;
    
  // new function 6: f'z = xxz - yyz

  mf(6, Momentum(2,0,1).local_index()) =  1;
  mf(6, Momentum(0,2,1).local_index()) = -1;

  ARRAY<Momentum>& sf = SYMINFO[3];
  sf.reset(mf.size1());
  for(l = 0; l < mf.size1(); l++)
    for(Momentum m(3); m(); ++m)
      if (fabs(mf(l,m.local_index())) > 0)
      {
	sf[l] = m;
	break;
      }

  // G-WAVE FUNCTIONS

  Mat& mg = TRANSFORM[4];
  mg.reset(9,15);
  mg.set(0);

  int no;

  // new function g0: 2 xxxx + 2 yyyy + 2 zzzz - 6 xxyy - 6 xxzz - 6 yyzz

  no = 0;
  mg(no, Momentum("xxxx").local_index()) =  2;
  mg(no, Momentum("yyyy").local_index()) =  2;
  mg(no, Momentum("zzzz").local_index()) =  2;
  mg(no, Momentum("xxyy").local_index()) = -6;
  mg(no, Momentum("xxzz").local_index()) = -6;
  mg(no, Momentum("yyzz").local_index()) = -6;

  // new function 1: g1 = 6 xxyz - yyyz - yzzz

  no = 1;
  mg(no, Momentum("xxyz").local_index()) =  6;
  mg(no, Momentum("yyyz").local_index()) = -1;
  mg(no, Momentum("yzzz").local_index()) = -1;

  // new function 2: g2 = 6 xyyz - xxxz - xzzz

  no = 2;
  mg(no, Momentum("xyyz").local_index()) =  6;
  mg(no, Momentum("xxxz").local_index()) = -1;
  mg(no, Momentum("xzzz").local_index()) = -1;

  // new function 3: g3 = 6 xyzz - xxxy - xyyy

  no = 3;
  mg(no, Momentum("xyzz").local_index()) =  6;
  mg(no, Momentum("xxxy").local_index()) = -1;
  mg(no, Momentum("xyyy").local_index()) = -1;

  // new function 4: g4 = 6 xxzz - xxxx + yyyy - 6 yyzz

  no = 4;
  mg(no, Momentum("xxzz").local_index()) =  6;
  mg(no, Momentum("xxxx").local_index()) = -1;
  mg(no, Momentum("yyyy").local_index()) =  1;
  mg(no, Momentum("yyzz").local_index()) = -6;

  // new function 5: g5 =  xxxz - xzzz

  no = 5;
  mg(no, Momentum("xxxz").local_index()) =  1;
  mg(no, Momentum("xzzz").local_index()) = -1;

  // new function 6: g6 =  yyyz - yzzzz

  no = 6;
  mg(no, Momentum("yyyz").local_index()) =  1;
  mg(no, Momentum("yzzz").local_index()) = -1;

  // new function 7: g7 =  xxxy - xyyy

  no = 7;
  mg(no, Momentum("xxxy").local_index()) =  1;
  mg(no, Momentum("xyyy").local_index()) = -1;

  // new function 8: g8 = 2 zzzz - xxxx - yyyy + 12 xxyy - 6 yyzz - 6 xxzz

  no = 8;
  mg(no, Momentum("zzzz").local_index()) =  2;
  mg(no, Momentum("xxxx").local_index()) = -1;
  mg(no, Momentum("yyyy").local_index()) = -1;
  mg(no, Momentum("xxyy").local_index()) = 12;
  mg(no, Momentum("yyzz").local_index()) = -6;
  mg(no, Momentum("xxzz").local_index()) = -6;

  ARRAY<Momentum>& sg = SYMINFO[4];
  sg.reset(mg.size1());
  for(l = 0; l < mg.size1(); l++)
    for(Momentum m(4); m(); ++m)
      if (fabs(mg(l,m.local_index())) > 0)
      {
	sg[l] = m;
	break;
      }


  // H-WAVE FUNCTIONS

  Mat& mh = TRANSFORM[5];
  mh.reset(11,30);
  mh.set(0);

  // new function h0:  xxxyz + xyyyz - 2 xyzzz        SYM: OOO 

  no = 0;
  mh(no, Momentum("xxxyz").local_index()) =   1;
  mh(no, Momentum("xyyyz").local_index()) =   1;
  mh(no, Momentum("xyzzz").local_index()) =  -2;

  // new function h1: xxxyz - xyyyz 

  no = 1;
  mh(no, Momentum("xxxyz").local_index()) =   1;
  mh(no, Momentum("xyyyz").local_index()) =  -1;

  // new function h2 = xxxxc - yyyyz - 2 xxzzz + 2 yyzzz  SYM: EEO  

  no = 2;
  mh(no, Momentum("xxxxz").local_index()) =  1;
  mh(no, Momentum("yyyyz").local_index()) = -1;
  mh(no, Momentum("xxzzz").local_index()) = -2;
  mh(no, Momentum("yyzzz").local_index()) =  1;

  // new function 3: g3 = xxxxz - 6 xxyyz + yyyyz

  no = 3;
  mh(no, Momentum("xxxxz").local_index()) =  1;
  mh(no, Momentum("yyyyz").local_index()) =  1;
  mh(no, Momentum("xxyyz").local_index()) = -6;

  // new function 4: g4 = 15 xxyyz - 5 xxzzz - 5 yyzzz + zzzz

  no = 4;
  mh(no, Momentum("xxyyz").local_index()) = 15;
  mh(no, Momentum("xxzzz").local_index()) = -5;
  mh(no, Momentum("yyzzz").local_index()) = -5;
  mh(no, Momentum("zzzzz").local_index()) =  1;

  // new function h2 = yyyyx - zzzzx - 2 yyxxx + 2 zzxxx  SYM: EEO  

  no = 5;
  mh(no, Momentum("yyyyx").local_index()) =  1;
  mh(no, Momentum("zzzzx").local_index()) = -1;
  mh(no, Momentum("yyxxx").local_index()) = -2;
  mh(no, Momentum("zzxxx").local_index()) =  1;

  // new function 3: g3 = yyyyx - 6 yyzzx + zzzzx

  no = 6;
  mh(no, Momentum("yyyyx").local_index()) =  1;
  mh(no, Momentum("zzzzx").local_index()) =  1;
  mh(no, Momentum("yyzzx").local_index()) = -6;

  // new function 4: g4 = 15 yyzzx - 5 yyxxx - 5 zzxxx + xxxx

  no = 7;
  mh(no, Momentum("yyzzx").local_index()) = 15;
  mh(no, Momentum("yyxxx").local_index()) = -5;
  mh(no, Momentum("zzxxx").local_index()) = -5;
  mh(no, Momentum("xxxxx").local_index()) =  1;

  // new function h2 = zzzzy - xxxxy - 2 zzyyy + 2 xxyyy  SYM: EEO  

  no = 8;
  mh(no, Momentum("zzzzy").local_index()) =  1;
  mh(no, Momentum("xxxxy").local_index()) = -1;
  mh(no, Momentum("zzyyy").local_index()) = -2;
  mh(no, Momentum("xxyyy").local_index()) =  1;

  // new function 3: g3 = azzzy - 6 zzxxy + xxxxy

  no = 9;
  mh(no, Momentum("zzzzy").local_index()) =  1;
  mh(no, Momentum("xxxxy").local_index()) =  1;
  mh(no, Momentum("zzxxy").local_index()) = -6;

  // new function 4: g4 = 15 zzxxy - 5 zzyyy - 5 xxyyy + yyyy

  no = 10;
  mh(no, Momentum("zzxxy").local_index()) = 15;
  mh(no, Momentum("zzyyy").local_index()) = -5;
  mh(no, Momentum("xxyyy").local_index()) = -5;
  mh(no, Momentum("yyyyy").local_index()) =  1;


  ARRAY<Momentum>& sh = SYMINFO[5];
  sh.reset(mh.size1());
  for(l = 0; l < mh.size1(); l++)
    for(Momentum m(5); m(); ++m)
      if (fabs(mh(l,m.local_index())) > 0)
      {
	sh[l] = m;
	break;
      }
}

Sphere::~Sphere() 
{
  for(int i=0; i< fields.size(); i++)
    if (fields[i]) delete fields[i];
  fields.reset(0);
}
