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

#define TEST

#ifdef TEST
#include <iostream>
#include <math.h>
using namespace std;
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#include "euler.h"

void euler(double* y, int n, void (*fcn)(double,double*,double*),
	   double a, double b, int nSteps) {
  int i,j;
  double h = (b-a)/nSteps;

  static int size=0;
  static double* f=0;
  if(n > size) {
    size = n;
    f = new double[size];
  }

  for(j=0; j<nSteps; j++) {
    fcn(a,y,f);
    for(i=0; i<n; i++) y[i] += h*f[i];
    a += h;
  }
}

#ifdef TEST
static const int n = 2;

// f(t,y) = (-2y[0],-3y[1]), d=2
void testFcn1(double t, double* y, double* f) {
  f[0] = -2*y[0];
  f[1] = -3*y[1];
}
void testY1(double t0, double t1, double* y) {
  double dt = t1-t0;
  y[0] = y[0]*exp(-2*dt);
  y[1] = y[1]*exp(-3*dt);
}

// f(t,y) = (1/2, -t), d=2
void testFcn2(double t, double* y, double* f) {
  f[0] = 0.5;
  f[1] = -t;
}
void testY2(double t0, double t1, double* y) {
  y[0] = y[0] + 0.5*(t1-t0);
  y[1] = y[1] - 0.5*(t1*t1-t0*t0);
}

// f(t,y) = (y[0]^2t,y[0]*y[1]), d=2
void testFcn3(double t, double* y, double* f) {
  f[0] = -y[0]*y[0]*t;
  f[1] = 0;
}
void testY3(double t0, double t1, double* y) {
  y[0] = y[0]/(1 + 0.5*(t1*t1-t0*t0)*y[0]);
}

// f(t,y) = (y[0]^2t,y[0]*y[1]), d=2
void testFcn4(double t, double* y, double* f) {
  f[0] = -y[0];
  f[1] = y[0];
}
void testY4(double t0, double t1, double* y) {
  double e = exp(t0-t1);
  double yy = y[0];
  y[0] = yy*e;
  y[1] = yy*(1-e)+y[1];
}

int main() {
  int i,j,nTimes,nSteps,testCase;
  double a,b;
  double y[n], yy[n];
  void (*fcn)(double,double*,double*);
  void (*trueY)(double,double,double*);

  // Enter data

 retry:
  cout << "Enter the test case number (1, 2, 3 or 4): ";
  cin >> testCase;
  switch(testCase) {
  case 1: {
    fcn = testFcn1;
    trueY = testY1;
    break;
  }
  case 2: {
    fcn = testFcn2;
    trueY = testY2;
    break;
  }
  case 3: {
    fcn = testFcn3;
    trueY = testY3;
    break;
  }
  case 4: {
    fcn = testFcn4;
    trueY = testY4;
    break;
  }
  default:
    goto retry;
  }

  cout << "Enter times a and b: ";
  cin >> a >> b;
  cout << "Enter the " << n << " initial conditions: ";
  for(i=0; i<n; i++) cin >> y[i];

  cout << "Enter number of outputs and number of euler steps between outputs: ";
  cin >> nTimes >> nSteps;

  // Loop over times to print solution

  double h = (b-a)/nTimes;
  double t0 = a, t1 = t0;
  for(i=0; i<n; i++) yy[i] = y[i];
  double error = 0;
  cout << "The solution is:" << endl;

  j=0;
  while(1) {
    trueY(t0,t1,yy);
    for(i=0; i<n; i++) error = max(error,fabs(yy[i]-y[i]));

    cout << "t = " << t1 << ":  ";
    for(i=0; i<n; i++) cout << y[i] << "  ";
    cout << "(  true  ";
    for(i=0; i<n; i++) cout << yy[i] << "  ";
    cout << ")" << endl;

    if(j==nTimes) break;
    j++;

    t0 = t1;
    t1 = t0+h;
    euler(y,n,fcn,t0,t1,nSteps);
  }

  cout << "ERROR: log(|error|) " << log(error)
    << ", log(1/h) " << -log(h/nSteps) << endl;
}
#endif
