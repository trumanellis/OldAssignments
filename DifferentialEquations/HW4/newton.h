///////////////////////////////////////////////////////////////////////////////
// NEWTON'S METHOD FOR 2-DIMENSIONAL FUNCTIONS
//
// This function uses Newton's method to solve the 2-dimensional system
//   f(x) = 0
// We assume that functions exist to find f(x) and f'(x) (called df(x)).
//
// Inputs:
//   x       The initial guess at the solution.
//   f       The name of the function ( interpret as f(input x, output value) )
//   df      The derivative of the function
//   tol     The convergence tolerance (must be > 0).
//   maxIter The maximum number of iterations that can be taken.
// Outputs:
//   x       The solution.
//   iter    The number of iterations taken.
// Return:
//   status  An error status code.
//      SUCCESS=0   Sucessful termination.
//      WONT_STOP=1 Error: Exceeded maximum number of iterations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEWTON_H
#define NEWTON_H

typedef double (vec)[2];
typedef double (mat)[2][2];

enum state { SUCCESS=0, WONT_STOP };

state newton(vec& x, void f(const vec&,vec&), void df(const vec&,mat&),
	     int& iter, double tol = 1e-6, int maxIter = 100);

#endif
