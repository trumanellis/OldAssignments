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

//#define NEWTON_TEST

#include <math.h>
#include "newton.h"
#ifdef NEWTON_TEST
#include <iostream.h>
#endif

inline static void invert(mat& a) {
  double detInv = 1/(a[0][0]*a[1][1] - a[1][0]*a[0][1]);
  double a00 = a[0][0];

  a[0][0] = a[1][1]*detInv;
  a[1][1] = a00*detInv;
  a[1][0] = -a[1][0]*detInv;
  a[0][1] = -a[0][1]*detInv;
}

inline static double normSq(const vec& x) {
  return x[0]*x[0] + x[1]*x[1];
}

state newton(vec& x, void f(const vec&,vec&), void df(const vec&,mat&),
	     int& iter, double tol, int maxIter) {
  vec dx, fx;
  mat dfx;
  tol *= tol; // We work with squares

  for(iter = 1; iter <= maxIter; iter++) {
    df(x,dfx);
    invert(dfx);
    f(x,fx);

    dx[0] = dfx[0][0]*fx[0] + dfx[0][1]*fx[1];
    dx[1] = dfx[1][0]*fx[0] + dfx[1][1]*fx[1];

    x[0] -= dx[0];
    x[1] -= dx[1];

    if(normSq(dx) < tol) return SUCCESS;
  }
  return WONT_STOP;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef NEWTON_TEST

void f(const double (&x)[2], double (&fx)[2]) {
  fx[0] = x[0] - exp(-x[0]);
  fx[1] = x[1]*(x[1]-5.0);
}

void df(const double (&x)[2], double (&dfx)[2][2]) {
  dfx[0][0] = 1.0 + exp(-x[0]);
  dfx[1][1] = 2*x[1]-5.0;
  dfx[1][0] = 0;
  dfx[0][1] = 0;
}

void main () {
  double x[2], fx[2], tol;
  int iter, maxIter;

  // INPUT

  cout << "Enter tolerance: ";
  cin >> tol;
  cout << "Enter max iterations: ";
  cin >> maxIter;
  cout << "Enter an initial guess: ";
  cin >> x[0] >> x[1];

  // NEWTON SOLVER

  int state = newton(x, f, df, iter, tol, maxIter);
  f(x,fx);

  // STATUS AND OUTPUT

  switch(state) {
  case SUCCESS:
    cout << "An approximate root is (" << x[0] << ", " << x[1] <<")\n";
    cout << "f(approximate root) is (" << fx[0] << ", " << fx[1] <<")\n";
    cout << "The number of iterations is " << iter << endl;
    break;
  case WONT_STOP:
  default:
    cout << "ERROR: Failed to converge in "<<maxIter<<" iterations!" << endl;
    break;
  }
}

#endif
