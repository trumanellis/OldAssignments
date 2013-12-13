/******************************************************************************
 * Solve a system of n equations satisfying
 *   y' = f(t,y), a < t <= b
 *   y(a) = y_a
 * 
 * Inputs:  y,       the initial solution at time a
 *          n,       number of equations
 *          fcn,     the name of the function in the problem
 *          a and b, the initial and final times
 *          nSteps,  the number of euler steps to take
 * Outputs: y,       the final solution at time b
 *****************************************************************************/

#ifndef EULER_H
#define EULER_H

void euler(double* y, int n, void (*fcn)(double,double*,double*),
	   double a, double b, int nSteps);

#endif
