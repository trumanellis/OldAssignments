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

//
// #define TEST

#ifdef TEST
#include <iostream>
#include <math.h>
using namespace std;
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif
#include "backEuler.h"

static vec yOld;
static double h;
static double t;
static void (*fcn)(double, const vec&, vec&);
static void (*dfcn)(double, const vec&, mat&);

void euler(vec& y, void (*f)(double,const vec&, vec&),
	   double a, double b, int nSteps) {
  int i,j;
  double h = (b-a)/nSteps;
  vec ff;

  for(j=0; j<nSteps; j++) {
    f(a,y,ff);
    y[0] += h*ff[0];
    y[1] += h*ff[1];
    a += h;
  }
}

void F(const vec& y, vec& Fy) {
  fcn(t,y,Fy);
  Fy[0] = y[0] - yOld[0] - h*Fy[0];
  Fy[1] = y[1] - yOld[1] - h*Fy[1];
  return;
}

void dF(const vec& y, mat& dFy) {
  dfcn(t,y,dFy);
  dFy[0][0] = 1 - h*dFy[0][0];
  dFy[1][0] =   - h*dFy[1][0];
  dFy[0][1] =   - h*dFy[0][1];
  dFy[1][1] = 1 - h*dFy[1][1];
  return;
}

state backEuler(vec& y, void (*f)(double, const vec&, vec&),
		void (*df)(double, const vec&, mat&),
		double a, double& b, int nSteps, double tol, int maxIter) {
  int i,j,iter;

  t = a;
  h = (b-a)/nSteps;
  fcn = f;
  dfcn = df;

  double fy[2];
  double dfy[2][2];

  for(j=0; j<nSteps; j++) {
    yOld[0] = y[0];
    yOld[1] = y[1];
    state st = newton(y,F,dF,iter,tol,maxIter);
    if(st != SUCCESS) {
      y[0] = yOld[0];
      y[1] = yOld[1];
      b = t;
      return st;
    }
    t += h;
  }
  return SUCCESS;
}

#ifdef TEST

// f(t,y) = (-y[0] + y[1],-100y[1])
void testFcn1(double t, const vec& y, vec& f) {
  f[0] = -y[0] + y[1];
  f[1] = -100*y[1];
}
void dTestFcn1(double t, const vec& y, mat& df) {
  df[0][0] = -1;
  df[1][0] = 0;
  df[0][1] = 1;
  df[1][1] = -100;
}
void testY1(double t, const vec& yInit, vec& y) {
  y[0] = (yInit[0] + yInit[1]/99.0)*exp(-t) - (yInit[1]/99.0)*exp(-100*t);
  y[1] = yInit[1]*exp(-100*t);
}

int main() {
  int i,j,nTimes,nSteps,testCase,maxIter;
  double a,b,tol;
  vec y, yy, yInit;
  void (*fcn)(double, const vec&, vec&);
  void (*dfcn)(double, const vec&, mat&);
  void (*trueY)(double, const vec&, vec&);

  // Enter data

 retry:
  cout << "Enter the test case number (1): ";
  cin >> testCase;
  switch(testCase) {
  case 1: {
    fcn = testFcn1;
    dfcn = dTestFcn1;
    trueY = testY1;
    break;
  }
  default:
    goto retry;
  }

  cout << "Enter times a and b: ";
  cin >> a >> b;
  cout << "Enter the 2 initial conditions: ";
  cin >> yInit[0] >> yInit[1];
  y[0] = yInit[0];
  y[1] = yInit[1];

  cout << "Enter number of outputs & number of backward euler steps per output: ";
  cin >> nTimes >> nSteps;

  cout << "Enter Newton tolerance and max number iterations: ";
  cin >> tol >> maxIter;

  // Loop over times to print solution

  double h = (b-a)/nTimes;
  double t0 = a, t1 = t0;
  double error = 0;
  cout << "The solution is:" << endl;

  j=0;
  while(1) {
    trueY(t1,yInit,yy);
    error = max(error,max(fabs(yy[0]-y[0]),fabs(yy[1]-y[1])));

    cout << "t = " << t1 << ":  ";
    cout << y[0] << "  " << y[1];
    cout << " (true  ";
    cout << yy[0] << "  " << yy[1];
    cout << ")" << endl;

    if(j==nTimes) break;
    j++;

    t0 = t1;
    t1 = t0+h;
    state st = backEuler(y,fcn,dfcn,t0,t1,nSteps,tol,maxIter);
    if(st != SUCCESS) {
      cout << "ERROR: Newton failed to converge!" << endl;
      return 1;
    }
  }

  cout << "ERROR: log(|error|) " << log(error)
    << ", log(1/h) " << -log(h/nSteps) << endl;
  return 0;
}
#endif
