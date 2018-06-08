/*
 ************************************************************************
 *
 *			  Numerical Math Package
 *
 * The present package implements various algorithms of Numerical Math
 *
 * $Id: math_num.h,v 1.1.1.1 2008/04/08 19:44:29 igor Exp $
 *
 ************************************************************************
 */

//#ifndef __GNUC__
//#pragma once
//#endif

#if !defined(MATH_NUM_H)
#define MATH_NUM_H 1

#ifdef __GNUC__
#if ! defined(TWINE)
//    #pragma interface
#endif
#endif

#include <float.h>
#include <math.h>
#include "haio.h"
#include "builtin.h"


//------------------------------------------------------------------------
//	Assertions, aborts, and related stderr printing

				// Print an error message on stderr and
				// abort, arguments like those of printf()
//void _error(const char * message,...);

				// Print a message on stderr
				// It looks and acts like printf(), only
				// prints on stderr (which is unbuffered...)
void message(const char * message,...);

#ifndef assert
#define assert(ex) \
        (void)((ex) ? 1 : \
              (fprintf(stderr,"Failed assertion " #ex " at line %d of `%s'.\n", \
               __LINE__, __FILE__), 0))
#define assertval(ex) assert(ex)
#endif

#define assure(expr,message)				\
	if	(expr) {;}				\
	else fprintf(stderr,"%s\n at line %d of '%s'.",message,__LINE__, __FILE__)

/*
 *------------------------------------------------------------------------
 *				Some constants
 * Compile and run the program epsilon.c to determine the values below for
 * your computer
 */

//#define EPSILON		2.22045e-16	// DBL_EPSILON
//#define SQRT_EPSILON	1.49012e-08

/*
 *------------------------------------------------------------------------
 *		Brent's minimum and zero finders for 
 *		  a function of a single argument
 */

			// A function of a single argument
class UnivariateFunctor
{
  public:
      UnivariateFunctor(); 
      virtual ~UnivariateFunctor();
      virtual double operator() (const double x);
      virtual int junk_func() { return 1; }
};


				// Find a minimum of function f
				// over the interval [a,b] with the
				// accuracy tol.
				// Returns an approx. to the min location
double fminbr(const double a, const double b, 
	      UnivariateFunctor& f, const double tol=DBL_EPSILON);

inline int nintf(const double x) { return (fmod(x,1.0) < 0.5) ? (int) x : ((int)x) + 1; } //!< return int closest to a given double

#endif // !defined (MATH_NUM_H)
