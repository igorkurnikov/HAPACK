/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: twoel.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "basis.h"
#include "function.h"
#include "integ_array.h"
#include "quad.h"
#include "vrr.h"
#include "hrr1.h"
#include "hrr2.h"
#include <numer.h>
#if !defined(_MSC_VER)
#include <sys/times.h>
#include <sys/resource.h>
#endif
#include "symop.h"
#include "integ_file.h"
#include "twoel.h"
#include "sphere.h"

static TaskCount twoeltasks;

TwoElectronIntegrals::TwoElectronIntegrals(
		     InternalBasis& bb1,InternalBasis& bb2,
		     InternalBasis& bb3,InternalBasis& bb4,
		     int ss12,int ss34,int ssp,IntegFile& oobuf) :
     b1(bb1),b2(bb2),b3(bb3),b4(bb4),same12(ss12),same34(ss34),
     same_pair(ssp),obuf(oobuf)
{
  timer -> start(t_2p_total);
  debug = 0;
  logic(); 
  timer -> stop(t_2p_total);
}

TwoElectronIntegrals::~TwoElectronIntegrals()
{
  for(int i =0 ; i< tconf.size(); i++)
    delete tconf[i];
  tconf.reset(0); // clear everything 

  cset12.reset(0,0);
  hset12.reset(0,0);
  cset34.reset(0,0);
  hset34.reset(0,0);
  
}

/*TEX
\subsection{Generating the Coefficients}  

Run through all allowed primitives and determine the combinations of symmetry operations which are required. 
*/
void TwoElectronIntegrals::generate_coefs()
{
  int n,m,l1,l2;
    
  cset12.reset(8,8);
  hset12.reset(8,8);
  for(n = 0; n < cset12.size1(); n++) 
  for(m = 0; m < cset12.size1(); m++) 
  {
    cset12(n,m).reset(b1.max_momentum(),b2.max_momentum());
    hset12(n,m).reset(b1.max_momentum(),b2.max_momentum());    
    cset12(n,m).set(0);
    hset12(n,m).set(0);
  }

  int i;
  for(i=0; i < primitive.size(); i++)
  {
    SymOp s1(primitive[i](0)); // obtain the indices of the symops.
    SymOp s2(primitive[i](1)); 
    Array2D<Coefficient_Set*>& cset  = cset12(s1,s2);
    Array2D<HRR_Set*>&         hset  = hset12(s1,s2);
    
    for(l1=0; l1 <  b1.max_momentum(); l1++)
      if (b1.noprim[l1] > 0)
      {
	//
        // the second angular momentum is smaller than the first 
	// if both functions come form the smae basis set.
        //
       int maxl2 = l1 + 1;
	if (!same12)
	  maxl2 = b2.max_momentum();

	for(l2=0; l2 < maxl2; l2++) 
	  if (b2.noprim[l2] > 0 && !cset(l1,l2))
	  {
	    cset(l1,l2) = new Coefficient_Set(b1.radials[l1],l1,s1,
					      b2.radials[l2],l2,s2,
					      false,TWOEL);
	    hset(l1,l2) = new HRR_Set(b1.fcenter[l1],s1,b2.fcenter[l2],s2,
				      false);
	  }
      }  
  }

    
  cset34.reset(8,8);
  hset34.reset(8,8);
  for(n = 0; n < cset34.size1(); n++) 
  for(m = 0; m < cset34.size1(); m++) 
  {
    cset34(n,m).reset(b3.max_momentum(),b4.max_momentum());
    hset34(n,m).reset(b3.max_momentum(),b4.max_momentum());    
    cset34(n,m).set(0);
    hset34(n,m).set(0);
  }

  for(i=0; i < primitive.size(); i++)
  {
    SymOp s1(primitive[i](2)); // obtain the indices of the symops.
    SymOp s2(primitive[i](3)); 
    Array2D<Coefficient_Set*>& cset  = cset34(s1,s2);
    Array2D<HRR_Set*>&        hset  = hset34(s1,s2);
    
    for(l1=0; l1 <  b3.max_momentum(); l1++)
      if (b3.noprim[l1] > 0)
      {
	//
        // the second angular momentum is smaller than the first 
	// if both functions come form the smae basis set.
        //
       int maxl2 = l1 + 1;
	if (!same34)
	  maxl2 = b4.max_momentum();

	for(l2=0; l2 < maxl2; l2++) 
	  if (b4.noprim[l2] > 0 && !cset(l1,l2))
	  {
	    cset(l1,l2) = new Coefficient_Set(b3.radials[l1],l1,s1,
					      b4.radials[l2],l2,s2,
                                              false,TWOEL);
	    hset(l1,l2) = new HRR_Set(b3.fcenter[l1],s1,b4.fcenter[l2],s2,
				      false);
	  }
      }  
  }
}


/*TEX

*/

void TwoElectronIntegrals::logic()
{
  generate_primitives();
  generate_coefs();
   
  int l1,l2,l3,l4;
  
  for(l1=0; l1 <  b1.max_momentum(); l1++)
  if (b1.noprim[l1] > 0)
  {
    //
    // the second angular momentum is smaller than the first 
    // if both functions come form the smae basis set.
    //
    int maxl2 = l1 + 1;
    if (!same12)
      maxl2 = b2.max_momentum();
    for(l2=0; l2 < maxl2; l2++) 
    if (b2.noprim[l2] > 0)
    {
      //
      // the third angular momentum is smaller than the first
      // if the basis pair 12 equals the basis pair 34
      //
      int maxl3 = b3.max_momentum();
      if (same_pair)
	maxl3 = l1 + 1;
      
      for(l3=0; l3 < maxl3; l3++)   
      if (b3.noprim[l3] > 0)
      {
	// 
	// the last index is doubly constrained:
	// if b3 == b4, the second index l4 < l3
        // if b1 == b3 and b2 == b4 the pair on the left must be larger
        // than that on the right.   
        // 
	int maxl4 = b4.max_momentum();
	if (same34)
	  maxl4 = l3 + 1;
	if (same_pair && l1 == l3)
	  maxl4 = l2 + 1;
	for(l4=0; l4 < maxl4; l4++)   
	if (b4.noprim[l4] > 0)
	{
	  // evaluate the integrals for the given shell
	  evaluate(l1,l2,l3,l4);
	}
      }
    }
  }
  clear_coefs();
}


/*TEX
Step 3: For the choosen quadruplet of angular momenta we now
compute the integrals. First we form the quadruple sets:
*/

void TwoElectronIntegrals::evaluate(int l1,int l2,int l3,int l4)
{
    
  if(!twoeltasks.local())
  {
    twoeltasks++;
    return;
  }
  twoeltasks++;
  // cout << " +++ TwoElectronIntegrals:: evaluate " << l1 << " " 
  // << l2 << " " << l3 << " " << l4 << " on node: " 
  // << pm -> id() << endl;
    
  generate_targets(l1,l2,l3,l4); 
  tconf.reset(target.size());
  int i;
  for(i = 0; i < tconf.size(); i++)
  {
    tconf[i] = new FourStorage(b1.radial_functions(l1).nofun(),
			       b2.radial_functions(l2).nofun(),
			       b3.radial_functions(l3).nofun(),
			       b4.radial_functions(l4).nofun(),1);
    assert(tconf[i]);
    if (!flaglist(Spherical))
      tconf[i] -> reset(l1,l2,l3,l4,T0);
    else 
      tconf[i] -> reset(l1,l2,l3,l4,T4);
  }

  for(int ip=0; ip < primitive.size(); ip++)
  {
    assert(cset12(primitive[ip](0),primitive[ip](1))(l1,l2) != 0);
    assert(cset34(primitive[ip](2),primitive[ip](3))(l3,l4) != 0);

    Quadruple_Set qset(*cset12(primitive[ip](0),primitive[ip](1))(l1,l2),
		       *cset34(primitive[ip](2),primitive[ip](3))(l3,l4),
		       false);

    // getrusage(0,&RUsage);	  
    // cout << " Maximum Resident Set Size " << RUsage.ru_maxrss << endl;
    
    VRR vrr(qset,l1,l2,l3,l4);
    vrr.clean_aux();
    
    if (l2 > 0 || l4 > 0) 
    {
      assert(hset12(primitive[ip](0),primitive[ip](1))(l1,l2) != 0);
      HRR1 hrr1(vrr,*hset12(primitive[ip](0),primitive[ip](1))(l1,l2),
		l1,l2,l3,l4);
      assert(hset34(primitive[ip](2),primitive[ip](3))(l3,l4) != 0);
      HRR2 hrr2(hrr1,*hset34(primitive[ip](2),primitive[ip](3))(l3,l4));
      distribute(l1,l2,l3,l4,primitive[ip],hrr2);
    }
    else
      distribute(l1,l2,l3,l4,primitive[ip],vrr);
    vrr.clean();
  }

  timer -> start(t_2p_output);
  for(int it=0; it < tconf.size(); it++)
  {
    Quadruple<SymDesignator>& symt = target[it];
    // cout << " +++ Output for Target: " << symt << endl;

    if (!flaglist(Spherical))
    {
      for(Momentum m1(l1);m1(); ++m1)
      for(Momentum m2(l2);m2(); ++m2)
      for(Momentum m3(l3);m3(); ++m3)
      for(Momentum m4(l4);m4(); ++m4)
      {
	Storage* tar = tconf[it] -> find(m1,m2,m3,m4);
	assert(tar != 0);

	tar -> output(b1,symt[0],b2,symt[1],b3,symt[2],b4,symt[3],
		      l1,m1.local_index(),l2,m2.local_index(),
		      l3,m3.local_index(),l4,m4.local_index(),
		      obuf,0,same12,same34,same_pair);
      }
    }
    else
    {
      for(int m1 = 0;m1 < 2*l1+1; ++m1)
      for(int m2 = 0;m2 < 2*l2+1; ++m2)
      for(int m3 = 0;m3 < 2*l3+1; ++m3)
      for(int m4 = 0;m4 < 2*l4+1; ++m4)
      {
	Storage* tar = tconf[it] -> find(m1,m2,m3,m4);
	assert(tar != 0);

	tar -> output(b1,symt[0],b2,symt[1],b3,symt[2],b4,symt[3],
		      l1,m1,l2,m2,l3,m3,l4,m4,
		      obuf,0,same12,same34,same_pair);
      }
      
    }
    delete tconf[it];
    tconf[it] = 0;
  }
  
  timer -> stop(t_2p_output);  

}
/*TEX
\subsection{Distribution}
*/
void TwoElectronIntegrals::distribute(int l1,int l2,int l3,int l4,
				      Quadruple<SymOp>& primitive,HRR2& hrr)
{
  timer -> start(t_2p_disth);
  Momentum zero;
  generate_sources(primitive);
  
  if (!flaglist(Spherical))
  {    
    for(Momentum m3(l3); m3();++m3)
    for(Momentum m4(l4); m4();++m4)
    {
      Storage* ss = hrr.find(m3,m4,0,zero);      
      assert(ss != 0);
            
      int bcount = 0;
      for(Momentum m1(l1); m1();++m1)
      for(Momentum m2(l2); m2();++m2)
      {	
	Storage field(*ss,bcount++,tconf[0] -> allocator() );
	assert(ss -> size1()  == field.size1());
	assert(ss -> size2()  == field.size2());
	assert(ss -> size3()  == field.size3());
	assert(ss -> size4()  == field.size4());
	assert(1 == field.blocks());
	
	for(int src = 0; src < sources.size(); src++)
	{
	  Quadruple<Momentum> m(m1,m2,m3,m4);	  
	  distribute(field,m,sources[src],operation[src]);	
	}
      }
    }
  }
  else
  {
    Sphere sphere(hrr,false,hrr.allocator());
    for(int m3 = 0; m3 < 2*l3+1;++m3)
    for(int m4 = 0; m4 < 2*l4+1;++m4)
    {
      IntegArray* ss = sphere.find(m3,m4,0,zero);
      assert(ss != 0);
      
      // now run over the momenta which are packed into the IntegArray
      int block = 0;
      for(int m1 = 0; m1 < 2*l1+1; ++m1)
      for(int m2 = 0; m2 < 2*l2+1; ++m2)
      {
	if (debug)
	  cout  << " TWOEL::outblock + sphere " << m1 << m2  
  	        << " starts at: " << block << endl;

	// Storage field(*ss,block++);
	Storage field(*ss,block++, tconf[0] -> allocator());
	for(int src = 0; src < sources.size(); src++)
	{
	  Quadruple<int> m(m1,m2,m3,m4);
	  // determine the effective momenta:
    
	  Quadruple<Momentum> meff(sphere.syminfo(l1,m[0]),
				   sphere.syminfo(l2,m[1]),
				   sphere.syminfo(l3,m[2]),
				   sphere.syminfo(l4,m[3]));
	  distribute(field,m,meff,sources[src],operation[src]);	
	}
      }
    }    
  }
  timer -> stop(t_2p_disth);
}

void TwoElectronIntegrals::distribute(int l1,int l2,int l3,int l4,
				      Quadruple<SymOp>& primitive,
				      VRR& vrr)
{
  timer -> start(t_2p_distv);
  generate_sources(primitive);
  assert(l2 == 0 && l4 == 0);

  Momentum zero;

  if (!flaglist(Spherical))
  {
    for(Momentum m1(l1); m1();++m1)
    for(Momentum m3(l3); m3();++m3)
    {
      Storage* ss = vrr.find(m1,m3,0,zero);
      assert(ss != 0);
            
      for(int src = 0; src < sources.size(); src++)
      {
	Quadruple<Momentum> m(m1,zero,m3,zero);
	distribute(*ss,m,sources[src],operation[src]);	
      }
    }
  }
  else
  {
    Sphere sphere(vrr,false,vrr.allocator());
    int zero = 0;
    
    for(int m1=0; m1 < 2*l1+1;m1++)
      for(int m3=0; m3 < 2*l3+1;m3++)
      {
	Storage* ss = sphere.find(m1,m3,0,zero);      
	assert(ss != 0);
	for(int src = 0; src < sources.size(); src++)
	{
	  Quadruple<int> m(m1,zero,m3,zero);

	  // determine the effective momenta:
    
	  Quadruple<Momentum> meff(sphere.syminfo(l1,m[0]),
				   sphere.syminfo(l2,m[1]),
				   sphere.syminfo(l3,m[2]),
				   sphere.syminfo(l4,m[3]));
	  
	  distribute(*ss,m,meff,sources[src],operation[src]);	
	}
      }
  }
  timer -> stop(t_2p_distv);
}


/*TEX
\subsection{Distribution of an Individual Block of Integrals} 

*/

void TwoElectronIntegrals::distribute(Storage& ss,Quadruple<int>& m,
				      Quadruple<Momentum>& meff,
				      Quadruple<SymOp>& source,SymOp& oper)
{
  Momentum zero;
  
  for(int itar = 0; itar < tconf.size(); itar++)
  {
    Storage *tar = tconf[itar] -> find(m[0],m[1],m[2],m[3]);
    if (tar == 0)
    {
      Storage* tt = new Storage(ss, tconf[itar] -> allocator());
      tt -> zero();
      
      assert(tt != 0);
      tconf[itar] -> insert(m[0],m[1],m[2],m[3],tt);
      tar = tt;
    }    
    double coef = compute_factor(target[itar],meff,source,oper);
    tar -> add_to(coef,ss);  
  }
}

void TwoElectronIntegrals::distribute(Storage& ss,Quadruple<Momentum>& m,
				      Quadruple<SymOp>& source,SymOp& oper)
{  
  Momentum zero;
  for(int itar = 0; itar < tconf.size(); itar++)
  {
    Storage *tar = tconf[itar] -> find(m[0],m[1],m[2],m[3]);
    if (tar == 0)
    {
      //Storage* tt = new IntegArray(ss.size1(),ss.size2(),
      //			   ss.size3(),ss.size4(),ss.blocks());
      Storage* tt = new Storage(ss, tconf[itar] -> allocator());
      assert(tt != 0);
      assert(tt -> size1() == ss.size1());
      assert(tt -> size2() == ss.size2());
      assert(tt -> size3() == ss.size3());
      assert(tt -> size4() == ss.size4());
      assert(tt -> blocks() == ss.blocks());
      tt -> zero();
      
      tconf[itar] -> insert(m[0],m[1],m[2],m[3],tt);
      tar = tt;
    }
    double coef = compute_factor(target[itar],m,source,oper);
    // cout << "distribute: " << source << " Op: " << oper << " to " << 
    //  target[itar] << " with: " << coef << endl;

    tar -> add_to(coef,ss);  
  }
}

/*TEX
\subsection{Cleanup of Coefficient Sets}
 
*/
void TwoElectronIntegrals::clear_coefs()
{
  int l1,l2;
  for(int n=0; n < cset12.size1(); n++)
  for(int m=0; m < cset12.size2(); m++)
  {
    for(l1=0; l1 <  b1.max_momentum(); l1++)
      if (b1.noprim[l1] > 0)
      {
	int maxl2 = l1 + 1;
	if (!same12)
	  maxl2 = b2.max_momentum();
	for(l2=0; l2 < maxl2; l2++) 
	  if (b2.noprim[l2] > 0)
	  {
	    delete cset12(n,m)(l1,l2);
	    delete hset12(n,m)(l1,l2);
	  }
      }

    cset12(n,m).reset(0,0);
    hset12(n,m).reset(0,0);

    for(l1=0; l1 <  b3.max_momentum(); l1++)
      if (b3.noprim[l1] > 0)
      {
	int maxl2 = l1 + 1;
	if (!same34)
	  maxl2 = b4.max_momentum();
	for(l2=0; l2 < maxl2; l2++) // s-functions only 
	  if (b4.noprim[l2] > 0)
	  {
	    delete cset34(n,m)(l1,l2);
	    delete hset34(n,m)(l1,l2); 
	  }  
      }

    cset34(n,m).reset(0,0);
    hset34(n,m).reset(0,0);

  }

}
