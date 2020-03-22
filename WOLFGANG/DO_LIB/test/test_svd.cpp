#include <math.h>
#include <errorip.h>
#include <assert.h>
#include <iomanip.h>
#include <numer.h>
#include <diagonalize.h>


int svd_test(D_Matrix& a,D_Matrix& u,D_Vector& w,D_Matrix& v)
{
  int i,j;
  cout << " +++ Test of Singular Value Decompostion " << endl;

  cout << " +++ Input Matrix 1 : " << endl;

  u = a;
  u.print(cout,5," %10.3f");
  svdcmp(u,w,v);

  cout << " +++ Diagonal Vector: " << endl;
  w.print(cout,5," %10.3f");
  cout << " +++ Output Matrix  : " << endl;
  v.print(cout,5," %10.3f");

  // Test

  assert(w.size() == a.size2());
  D_Matrix diag(a.size2(),a.size2(),a.offset2(),a.offset2());
  diag.set(0);
  for(i = a.offset2(); i < a.offset2() + a.size2(); i++)
    diag(i,i) = w(i);

  // diag.print(cout,5," %10.3f");

  D_Matrix tmp1(a.size2(),a.size2(),a.offset2(),a.offset2());
  tmp1.set(0);
  tmp1.multiply_transpose2(diag,v);
  
  // tmp1.print(cout,5," %10.3f");

  D_Matrix tmp2(a.size1(),a.size2(),a.offset1(),a.offset2());
  tmp2.set(0);
  tmp2.multiply(u,tmp1);
  
  // tmp2.print(cout,5," %10.3f");
    
  int err = 0;
  for(i = a.offset1(); i < a.offset1() + a.size1(); i++)
  for(j = a.offset2(); j < a.offset2() + a.size2(); j++)
    if (fabs(a(i,j) - tmp2(i,j)) > 1E-10)
    {
      cout << " --- Error in SVD: " << setw(3) << i << " "
	   << setw(3) << j << " A(i,j) = " << a(i,j) << " " 
	   << " PROD(i,j) = " << tmp2(i,j) << endl;
      err++;
    }
  
  if  (!err)
    cout << " +++ SVD Sucessful. " << endl;
  return err;
}

double solve(D_Matrix& m,D_Matrix& o,D_Vector& eval,D_Matrix& evec)
{
  int real_size = m.size1();
  
  D_Matrix u(real_size,real_size);
  orthonormalize(o,u);

  D_Matrix utrans(u);
  utrans.transpose();

  // form Heff in M by two transforms with U

  D_Matrix h1(real_size,real_size);
  D_Matrix h2(real_size,real_size);
  h1.set(0);
  h1.multiply(m,u);
  h2.set(0);
  h2.multiply(u,h1);

  D_Vector snew(real_size);  
  D_Matrix tempvec(real_size,real_size);  

  
  ::diagonalize(h2,tempvec,eval);

  // the eigenvector in the old basis
  for(int k = 0; k < real_size; k++)
  {
    int i;
    for(i=0; i< real_size; i++)
      snew[i] = tempvec(i,k);

    double new_energy = eval[k]; 

    D_Vector test(real_size);
  
    test.set(0);
    test.multiply(h2,snew);
    for(i=0; i< real_size; i++)
      if (fabs(test[i] - new_energy*snew[i]) > 1E-8)
      {
	char buf[100];
	sprintf(buf,"Vec: %3i Index: %3i LHS: %10.6f RHS: %10.6f DIFF: %10.3e",
		k,i,test[i],new_energy*snew[i],test[i]-new_energy*snew[i]);
	cout << " +++ Diag1: " << buf << endl;
      }

    cout << "U:" << u << endl;
    D_Vector s(real_size);
    s.set(0);
    s.multiply(u,snew);    
    test.set(0);
    test.multiply(m,s);
    D_Vector tr(real_size);
    tr.set(0);
    tr.multiply(o,s);
    
    double normfac = s*tr;
    tr *= new_energy;
    
    for(i =0 ; i < real_size; i++)
      if (fabs(test[i]-tr[i]) > 1E-8)
      {
	char buf[100];
	sprintf(buf,"Vec: %3i Index: %3i LHS: %10.6f RHS: %10.6f DIFF: %10.3e",
		k,i,test[i],tr[i],test[i]-tr[i]);
	cout << " +++ Diag2: " << buf << endl;
      }

    for(i =0 ; i< real_size;i++)
      evec(i,k) = s(i)/sqrt(normfac);
  }
  return 0;
  
}


/*
\subsection{Constrained Minimization}

Given is the hamiltonian matrix \name{H} in a non-orthogonal basis $\phi_i$ of size N with overlap 
matrix \name{S}. We seek to minimize 
$$\frac{\langle  Psi | H | Psi \rangle}{\langle  Psi | Psi \rangle}$$

subject to the constraint that the solution is orthogonal to a set of K vectors $\psi_k$.
Let $sx_{ki} = \langle \psi_k | \phi_i \rangle$ the external overlap matrix. 

The problem is solved in three steps. 

\begin{itemize}
\item[(1)] We embed the $N \times K$ matrix $sx$ (with $K \leq N$) in a $N \times N$ matrix $X$
           by filling the lower k rows with zeroes. We then use singular value decomposition (SVD)
	   to construct the nullspace of this matrix, i.e. the set of solution vectors
                    
           $$         X \vec{a}_i = 0 $$

	   this set (N-K) of vectors defines a basis 

	   $$         \phi'_i = \sum a_{ik} \phi_k  $$

	   of the Hilbert space which satisfies  the constraints. 

\item[(2)] We express the Hamiltonian and the overlap matrices in the new basis of size N-K
           and use SVD to orthogonalize this matrix. This procedure defines $M \leq N-K$ 
	   new basis vectors

	   $$         \phi''_i = \sum b_{ik} \phi'_k  $$

	   wich are both orthonormal and satisfy the constraints. 

\item[(3)] We express the Hamiltonian in the new basis and find the lowest eigenvalue. 
           Let the corrspounf eigenvector be $\vec{u}$.  We compute the coefficients for this 
	   eigenvector in terms of the original basis as:

	   $$         r_l = \sum u_i b_{ik} a_{kl}   $$
                      
	   */

double constrained_diag(D_Matrix& h,D_Matrix& s,D_Matrix& sx,D_Vector& v_res)
{
  int debug = 1;

  int N =  h.size1();
  int K = sx.size1();
  assert( h.size2() ==  N);
  assert( s.size1() ==  N);
  assert( s.size2() ==  N);
  assert(sx.size2() ==  N);

  if (K > N)
  {
    cout << " Constrained-Diag Error: More constraints than degrees of freedom. " << endl;
    error("Fatal Error in constrained diagonalization. ");
  }

  // step 1

  int i,j,k,l;
  D_Matrix X(N,N);
  X.set(0);

  for(i = 0; i < K; i++)
  for(j = 0; j < N; j++)
    X(i,j) = sx(i,j);

  D_Matrix u(X);
  D_Matrix v;
  D_Vector w;
  
  svdcmp(u,w,v);

  int count = 0;
  for(i = 0; i < N; i++)
    count += fabs(w(i)) < 1E-10;
  
  if (count == 0)
  {
    cout << " Constrained-Diag Error: No effective degrees of freedom. " << endl;
    error("Fatal Error in constrained diagonalization. ");
  }
  

  D_Matrix a(count,N);
  count = 0;
  for(i = 0; i < N; i++)
    if (fabs(w(i)) < 1E-10) // found a solution
    {
      D_Vector sol(N);
      for(j = 0; j < N; j++)
      {
	a(count,j) = v(j,i);
	sol[j]     = v(j,i);
      }
      count++;
      
      D_Vector tmp(N);
      tmp.set(0);
      tmp.multiply(X,sol);
      
      double err = tmp*tmp;
      if (fabs(err) > 1E-10)
      {
	cout << " --- Error: The following vector is not a solution in step 1 " << endl;
	sol.print(cout,5," %10.3f");	
	error("Fatal Error in constrained diagonalization. ");
      }
    }

  int NN = count;
  D_Matrix hp(NN,NN);
  D_Matrix sp(NN,NN);

  for(i = 0; i < NN; i++)
  for(j = 0; j < NN; j++)
  {
    double hsum = 0;
    double ssum = 0;
    for(k = 0; k < N; k++)
    for(l = 0; l < N; l++)
    {
      hsum += a(i,k)*h(k,l)*a(j,l);
      ssum += a(i,k)*s(k,l)*a(j,l);
    }
    hp(i,j) = hsum;
    sp(i,j) = ssum;
  }

  // step 2

  u = sp;
  svdcmp(u,w,v);
  
  int M = 0;
  for(i = 0; i < NN; i++)
    if (fabs(w(i)) > 1E-10) // found a solution
      M++;

  D_Matrix b(M,NN);
  count = 0; 
  for(i = 0; i < NN; i++)
    if (fabs(w(i)) > 1E-10) // found a solution
    {
      for(j = 0; j < NN; j++)
	b(count,j) = u(j,i) / sqrt(w(i));
      count++;
    }

  D_Matrix hpp(M,M);

  for(i = 0; i < M; i++)
    for(j = 0; j < M; j++)
    {
      double hsum = 0;
      double ssum = 0;
      for(k = 0; k < NN; k++)
      for(l = 0; l < NN; l++)
      {
	ssum += b(i,k)*sp(k,l)*b(j,l);
	hsum += b(i,k)*hp(k,l)*b(j,l);
      }

      if (fabs(ssum - (i==j)) > 1E-10) 
      {
	cout << " Error in Overlap : " << i << " " << j << " " << ssum << endl;
	error("Fatal Error in constrained diagonalization. ");
      }
      hpp(i,j) = hsum;
    }
  
  // step 3:

  D_Matrix evec(M,M);
  D_Vector eval(M);

  diagonalize(hpp,evec,eval);
  
  double energy = eval(0);
  
  // now checking:

  D_Vector vpp(M);
  for(i = 0; i < M; i++)
    vpp[i] = evec(i,0);

  D_Vector tpp1(M);
  D_Vector tpp2(vpp);
  tpp1.set(0);
  tpp1.multiply(hpp,vpp);
  tpp2 *= energy;
  
  tpp1 -= tpp2;
  double err_pp = tpp1*tpp1;
  if (fabs(err_pp) > 1E-10)
  {
    cout << " +++ Error in Matrix Diag: " << energy << endl;
    hpp.print(cout,5," %10.6f");    
    vpp.print(cout,5," %10.6f");
    tpp1.print(cout,5," %10.6f");
    tpp2.print(cout,5," %10.6f");
    error("Fatal Error in constrained diagonalization. ");
  }
  
  // now transform back into the non-orthogonal basis

  D_Vector vp(NN);
  vp.set(0);
  for(i = 0; i < NN; i++)
  for(j = 0; j <  M; j++)
    vp[i] += b(j,i) * vpp[j];
 
  D_Vector tp1(NN);
  D_Vector tp2(NN);
  tp1.set(0);
  tp2.set(0);
  tp1.multiply(hp,vp);
  tp2.multiply(sp,vp);
  tp2 *= energy;

  tp1 -= tp2;
  double err_p = tp1*tp1;
  if (fabs(err_p) > 1E-10)
  {
    cout << " +++ Error in Nonorthogonal Basis Diag: " << energy << endl;

    cout << "HP:  " << endl;    hp.print(cout,5," %10.6f");    
    cout << "SP:  " << endl;    sp.print(cout,5," %10.6f");    
    cout << "B :  " << endl;     b.print(cout,5," %10.6f");    

    cout << "VP:  " << endl;    vp.print(cout,5," %10.6f");
    cout << "VPP: " << endl;   vpp.print(cout,5," %10.6f");


    D_Matrix vvv(NN,NN);
    D_Vector vve(NN);
    D_Vector vv (NN);
    solve(hp,sp,vve,vvv);
    cout << "   Computed EV: " << vve(0) << endl;
    for(i = 0; i < NN; i++)
      vv[i] = vvv(i,0);

    cout << "VVV: " << endl;   vv.print(cout,5," %10.6f");

    D_Vector tx1(NN);
    D_Vector tx2(NN);
    tx1.set(0);
    tx2.set(0);
    tx1.multiply(hp,vv);
    tx2.multiply(sp,vv);
    tx2 *= vve(0);
    tx1 -= tx2;
    cout << " TEST2: " << tx1*tx1 << " " << err_p <<  endl;
    error("Fatal Error in constrained diagonalization. ");

  }

  // now transform into the original basis

  v_res.reset(N);
  v_res.set(0);
  for(i = 0; i <  N; i++)
  for(j = 0; j < NN; j++)
    v_res[i] += vp[j] * a(j,i);
 
  D_Vector t1(N);
  D_Vector t2(N);
  t1.set(0);
  t2.set(0);
  t1.multiply(h,v_res);
  t2.multiply(s,v_res);
  t2 *= energy;

  t1 -= t2;
  double err = t1*t1;
  if (fabs(err) > 1E-10)
  {
    cout << " +++ Error in Nonorthogonal Unconstrained Basis Diag: " << energy << endl;
    h.print(cout,5," %10.6f");    
    s.print(cout,5," %10.6f");    
    v_res.print(cout,5," %10.6f");
    t1.print(cout,5," %10.6f");
  }

  return energy;
  
}

main()
{
  int i,j;
  int m = 5;
  int n = 3;
  D_Matrix a(m,n);
  D_Matrix v;
  D_Vector w;

  for(i = a.offset1(); i < a.offset1() + a.size1(); i++)
  for(j = a.offset2(); j < a.offset2() + a.size2(); j++)
    a(i,j) = 2*(i+1)*(j+3) + 3*i*(i+1)/(j+2);

  D_Matrix u;
  svd_test(a,u,w,v);

  // orthogonal space
  // given is a set of M equations in N variables (M < N). There are
  // N-M linearly independent solutions to this set of M equations.

  m = 0;
  n = 5;

  a.reset(n,n);
  a.set(0);
  for(i = 0; i < m; i++)
  for(j = 0; j < n; j++)
    a(i,j) = 2*(i+1)*(j+3) + 3*i*(i+1)/(j+2);
  
  svd_test(a,u,w,v);

  int count = 0;
  for(i = 0; i < n; i++)
    if (fabs(w(i)) < 1E-10) // found a solution
    {
      cout << " +++ Solution Vector: " << setw(3) << i << " W = " << w(i)
	<< endl;
      D_Vector sol(n);
      for(j = 0; j < n; j++)
	sol[j] = v(j,i);
      
      D_Vector tmp(n);
      tmp.set(0);
      tmp.multiply(a,sol);

      double err = tmp*tmp;
      if (fabs(err) > 1E-10)
      {
	cout << " --- Error: The following vector is not a solution: " << endl;
	tmp.print(cout,5," %10.3f");	
      }
      else
	count++;
    }
  
  cout << " +++ There were " << count << " linearly independent solututions. "
       << endl;

  /*
    Let S_{ij} be the overlap matrix of the given basis and S = U W V^+ its 
    SVD decomposition. We represent the basis

                   c_i = V W^{-1/2} U^T c_i'

    and the overlap matrix in the new basis c_i' becomes:

                U W^{-1/2} V^T (U W V^+) V W^{-1/2} U^T
                U W^{-1/2} V^T  U W^{1/2} U^T

  */
  cout << " +++ Test to construct an orthornomal basis " << endl;

  n = 2;
  m = 1;
  D_Matrix h (n,n);
  D_Matrix sx(m,n);
  a.reset(n,n);
  a.set(0);
  for(i = 0; i < n; i++)
  for(j = 0; j <= i; j++)
  {
    a(i,j) = (0.02*(i+1)*(j+3) + 0.053*i*(i+1)/(j+2) + 0.4*(i+1)*(i==j));
    a(j,i) = a(i,j);
    h(i,j) = 4*(i+3)*(j+2) + 3*i*(j+1)/(i+2);
    h(j,i) = h(i,j);
  }

  for(i = 0; i < m; i++)
  for(j = 0; j < n; j++)
    sx(i,j) = 0.1221*(n+1)*(n+2) + 3*m;

  D_Matrix evec(n,n);
  D_Vector eval(n);

  diagonalize(a,evec,eval);
  cout << "EigenValues: " << eval << endl;

  svd_test(a,u,w,v);

  count = 0;
  for(i = 0; i < n; i++)
    if (fabs(w(i)) > 1E-10) // found a solution
      count++;

  D_Matrix new_basis(count,n);
  count = 0; 
  for(i = 0; i < n; i++)
    if (fabs(w(i)) > 1E-10) // found a solution
    {
      cout << " +++ Solution Vector: " << setw(3) << i << " W = " << w(i)
	<< endl;
      for(j = 0; j < n; j++)
	new_basis(count,j) = u(j,i);
      count++;
    }
  
  for(i = 0; i < count; i++)
    for(j = 0; j < i; j++)
    {
      double ovl1 = 0;
      double ovl2 = 0;
      for(int k = 0; k < n; k++)
      for(int l = 0; l < n; l++)
	ovl1 += new_basis(i,k)*a(k,l)*new_basis(j,l);
      cout << " Overlap : " << i << " " << j << " " << ovl1 << endl;
    }


  cout << " **************** Cecking Constrained Diag ************************ " << endl;
  cout << " **************** Cecking Constrained Diag ************************ " << endl;
  cout << " **************** Cecking Constrained Diag ************************ " << endl;
  cout << " **************** Cecking Constrained Diag ************************ " << endl;

  D_Vector vres;
  double ee = constrained_diag(h,a,sx,vres);
    
}

