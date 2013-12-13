///////////////////////////////////////////////////////////////////////////////
// Solve a system of 2 equations satisfying
//   y' = f(t,y), a < t <= b
//   y(a) = y_a
// using backward euler
//     y_{n+1} = y_n + h f(t_{n+1},y_{n+1})
//   ==>
//     F(y_{n+1}) = 0  where  F(y) = y - y_n - h f(t_{n+1},y)
//                            dF(y) = I - h Df(t_{n+1},y)
// 
// Inputs:  y,       the initial solution at time a
//          fcn,     the name of the function in the problem
//          dfcn,    the name of the derivative of the function in the problem
//          a and b, the initial and final times
//          nSteps,  the number of euler steps to take
//          tol,     the Newton tolerance
//          maxIter, the max number of iterations allowed for Newton
// Outputs: y,       the final solution at time b
//          b,       the final time
// Return:  state,   the return state
//       SUCCESS=0   Sucessful termination.
//       WONT_STOP=1 Error: Exceeded maximum number of iterations.
//
// Also solve with Euler's method
//   y_{n+1} = y_n + h f(t_n,y_n)
//
// Inputs:  y,       the initial solution at time a
//          fcn,     the name of the function in the problem
//          a and b, the initial and final times
//          nSteps,  the number of euler steps to take
// Output:  y,       the final solution at time b
///////////////////////////////////////////////////////////////////////////////

#ifndef BACKEULER_H
#define BACKEULER_H

#include "newton.h"

state backEuler(vec& y, void (*f)(double, const vec&, vec&),
		void (*df)(double, const vec&, mat&),
		double a, double& b, int nSteps, double tol, int maxIter);

void euler(vec& y, void (*f)(double, const vec&, vec&),
	   double a, double b, int nSteps);

#endif
