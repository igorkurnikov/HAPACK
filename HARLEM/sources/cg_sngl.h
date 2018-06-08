//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
#ifndef CG_SNGL_H
#define CG_SNGL_H


template < class GMatrix, class GVector, class Preconditioner, class Real >
int 
CG_solv1(const GMatrix &A, GVector &x,  GVector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
// Solution for one vector 
{
  Real resid;
  GVector p, z, q;
  Real alpha, beta, rho, rho_1;

  Real normb = norm(b);
  GVector r = b - A*x;

  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return FALSE;
  }

  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    rho = dot(r, z);
    
    if (i == 1)
      p = z;
    else {
      beta = rho / rho_1;
      p = z + beta * p;
    }
    
    q = A*p;
    alpha = rho / dot(p, q);
    
    x += alpha * p;
    r -= alpha * q;
    
    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }

    char buf[132]; 
 	sprintf(buf,"resid=norm(r)/ normb at iteration %d  is %16.9f \n",i,resid );
	PrintLog("%s",buf);

    rho_1 = rho;
  }
  
  tol = resid;
  return 1;
}


template < class GMatrix, class GVector, class Preconditioner, class Real >
int 
CG_sngl(const GMatrix &A, GVector &x,  GVector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
// Solution for one vector 
{
  Real resid;
  GVector p, z, q;
  Real alpha, beta, rho, rho_1;

  Real normb = norm(b);
  GVector r = b - A*x;

  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return FALSE;
  }

  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    rho = dot(r, z);
    
    if (i == 1)
      p = z;
    else {
      beta = rho / rho_1;
      p = z + beta * p;
    }
    
    q = A*p;
    alpha = rho / dot(p, q);
    
    x += alpha * p;
    r -= alpha * q;
    
    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }

    char buf[132]; 
 	sprintf(buf,"resid=norm(r)/ normb at iteration %d  is %16.9f \n",i,resid );
	PrintLog("%s",buf);

    rho_1 = rho;
  }
  
  tol = resid;
  return 1;
}

template < class GMatrix, class GVector, class Preconditioner, class Real >
int 
CG_mult(GMatrix &A, GVector &x,  GVector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
// Solution for one vector 
{
  Real resid;
  GVector p, z, q;
  Real alpha, beta, rho, rho_1;

  Real normb = norm2(b); 
  GVector tmp = A*x; 
  GVector r = b - tmp;  // r = b - Ax 

  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = norm2(r) / normb) < tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    rho = dot2(r, z);
    
    if (i == 1)
      p = z;
    else 
    {
      beta = rho / rho_1;
      tmp = beta*p;
      p = z + tmp;        // p = z + (rho / rho_1)* p
    }
    
    q = A*p;
    alpha = rho / dot2(p, q);
    
    tmp = alpha * p; 
    x += tmp;  
    tmp = alpha * q;
    r -= tmp;
    
    if ((resid = norm2(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }

//    PrintLog("resid=norm(r)/ normb at iteration %d \n is %12.6e \n",i,(double)resid);

    rho_1 = rho;
  }
  
  tol = resid;
  return 1;
}


#endif /* !CG_SNGL_H */
