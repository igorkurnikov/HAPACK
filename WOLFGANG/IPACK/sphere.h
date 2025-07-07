/*TEX
\section{Transformation to Spherical Gaussians}

\subsection{Concepts}

In this section we discuss the transformation from cartesian to
shperical gaussian functions. This transformation can generally be
represented as:
\be
\phi_l (x) = \sum{\vec{m}} \alpha_{l\vec{m}} \phi_{\vec{m}} (x) 
\ee
where $0 \le l < 2L +1$ is the linear index of the spherical function and $\vec{m}$
designates the set of cartesian angular momenta of total angular
momentum $L$.  Accordingly a set of integrals dependent on two angular
momenta is transformed in two steps:
\be
   i'(l,\vec{m}) = \sum{\vec{n}} \alpha_{l\vec{n}} i(\vec{n},\vec{m})
\ee
and
\be
   i''(l,l') = \sum{\vec{n}} \alpha_{l'\vec{n}} i(l,\vec{n})
\ee
Since this transformation depends solely on the angular characteristic
of the integrals it can be immedeately extended to entire
\name{IntegArrays}, which is the formulation of the  two classes
presented below.
 
\subsection{Class Header}  

*/  

#ifndef SPHERE_H
#define SPHERE_H

#include "tquadruple.h"
#include "four_stor.h"
#include "storage.h"

class RecursionStorage;

class HalfSphere
{
  int debug;
  int l1;
  int l2;
  int no_l1;
  int no_l2;
  int max_aux;
  int max_mu;
  int no_mu;

  Allocator& alloc;
  ARRAY<IntegArray*> fields; 

int compute_index(const int l1, const Momentum& l2,
		  const int aux, const Momentum& mu) const
{
  assert( ((l1*no_l2+l2.local_index())*(max_aux+1) + aux)*no_mu + mu.index()
    < fields.size());
  return ((l1*no_l2+l2.local_index())*(max_aux+1) + aux)*no_mu + mu.index();
}
  
public:
            HalfSphere(RecursionStorage& rec,int dbx,Allocator& tmpl);
           ~HalfSphere();

            // return a pointer to the specified IntegArray
IntegArray* find(const int l1, const Momentum& l2,const int aux, const Momentum& mu)
            { return fields[compute_index(l1,l2,aux,mu)]; }   
}; 

class Sphere
{
  friend class HalfSphere;
  Allocator& alloc;  
  ARRAY<IntegArray*> fields; 
  int debug;
  int l1;
  int l2;
  int no_l1;
  int no_l2;
  int max_aux;
  int max_mu;
  int no_mu;

int compute_index(const int l1, const int l2,const int aux, 
		  const Momentum& mu) const
		  { return ((l1*no_l2+l2)*(max_aux+1)+aux)*no_mu +mu.index(); }
    static ARRAY<ARRAY<Momentum> > SYMINFO;
    static ARRAY<Mat>              TRANSFORM;
public:
            Sphere(RecursionStorage& rec,int dbx,Allocator& tmpl);
           ~Sphere();
static void transform_init();
static void transform_clear();

static Momentum&  
            syminfo(int l,int no) 
	    { 
	      if (Sphere::TRANSFORM.size() == 0)
	      Sphere::transform_init();
	      return SYMINFO[l][no]; 
	    } 
            // return a pointer to the specified IntegArray
IntegArray* find(const int l1,const int l2,const int aux, const Momentum& mu)
            { return fields[compute_index(l1,l2,aux,mu)]; }   
}; 

#endif


