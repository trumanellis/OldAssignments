/******************************************************************************
 * Solve a system of 2 equations satisfying
 *   y' = f(t,y), a < t <= b
 *   y(a) = y_a
 * using euler and backward euler to study stability
 *****************************************************************************/

#include <iostream>
#include <iomanip>
#include <math.h>
#define max(a,b) ((a) >= (b) ? (a) : (b))
#include "backEuler.h"
#include <stdio.h>
//#include "euler.h"
using namespace std;

// f(t,y) = (-y[0] + y[1],-100y[1])
void testFcn1(double t, const vec& y, vec& f) {
  f[0] = -y[0] + y[1];
  f[1] = -100*y[1];
}
void dTestFcn1(double t, const vec& y, mat& df) {
  df[0][0] = -1;
  df[0][1] = 1;
  df[1][0] = 0;
  df[1][1] = -100;
}
void testY1(double t, const vec& yInit, vec& y) {
  y[0] = (yInit[0] + yInit[1]/99.0)*exp(-t) - (yInit[1]/99.0)*exp(-100*t);
  y[1] = yInit[1]*exp(-100*t);
}

// f(t,y) = (197*y[0]-594*y[1],99*y[0]-298*y[1])
void testFcn2(double t, const vec& y, vec& f) {
  f[0] = 197*y[0] - 594*y[1];
  f[1] =  99*y[0] - 298*y[1];
}
void dTestFcn2(double t, const vec& y, mat& df) {
  df[0][0] =  197;
  df[0][1] = -594;
  df[1][0] =   99;
  df[1][1] = -298;
}
void testY2(double t, const vec& yInit, vec& y) {
  double E1 = exp(-t), E2 = exp(-100*t);
  y[0] = 3*(yInit[0] - 2*yInit[1])*E1 + 2*(-yInit[0] + 3*yInit[1])*E2;
  y[1] =   (yInit[0] - 2*yInit[1])*E1 +   (-yInit[0] + 3*yInit[1])*E2;
}

// f(t,y) = (-y[0]*y[0],-y[0]*y[1])
void testFcn3(double t, const vec& y, vec& f) {
  f[0] = -y[0]*y[0];
  f[1] = -y[0]*y[1];
}
void dTestFcn3(double t, const vec& y, mat& df) {
  df[0][0] = -2*y[0];
  df[0][1] = 0;
  df[1][0] = -y[1];
  df[1][1] = -y[0];
}
void testY3(double t, const vec& yInit, vec& y) {
  y[0] = yInit[0]/(yInit[0]*t + 1);
  y[1] = yInit[1]/(yInit[0]*t + 1);
}

int main() {
  int i,j,nTimes,nSteps,testCase,maxIter;
  double a,b,tol;
  vec y_beuler, y_euler, y_true, yInit;
  void (*fcn)(double, const vec&, vec&);
  void (*dfcn)(double, const vec&, mat&);
  void (*trueY)(double, const vec&, vec&);

  // Enter data

 retry:
  cout << "Enter the test case number (1, 2 or 3): ";
  cin >> testCase;
  switch(testCase) {
  case 1: {
    fcn = testFcn1;
    dfcn = dTestFcn1;
    trueY = testY1;
    break;
  }
  case 2: {
    fcn = testFcn2;
    dfcn = dTestFcn2;
    trueY = testY2;
    break;
  }
  case 3: {
    fcn = testFcn3;
    dfcn = dTestFcn3;
    trueY = testY3;
    break;
  }
  default:
    goto retry;
  }

  cout << "Enter times a and b: ";
  cin >> a >> b;
  cout << "Enter the 2 initial conditions: ";
  cin >> yInit[0] >> yInit[1];
  y_beuler[0] = yInit[0];
  y_beuler[1] = yInit[1];
  y_euler[0] = yInit[0];
  y_euler[1] = yInit[1];

  cout << "Enter number of outputs & number of backward euler steps per output: ";
  cin >> nTimes >> nSteps;

  //cout << "Enter Newton tolerance and max number iterations: ";
  //cin >> tol >> maxIter;
  tol = 1e-6; maxIter = 20;

  // Loop over times to print solution

  double h = (b-a)/nTimes;
  double t0 = a, t1 = t0;
  double error_euler = 0;
  double error_beuler = 0;

  cout.setf(ios::scientific);
  cout.setf(ios::right);
  cout.precision(4);
  cout << "The solution is:" << endl;

  j=0;
  while(1) {
    trueY(t1,yInit,y_true);
    //    printf("true= %g\t\t%g\n", y_true[0], y_true[1]);
    error_euler = max(error_euler,max(fabs(y_true[0]-y_euler[0]),
				      fabs(y_true[1]-y_euler[1])));
    error_beuler = max(error_beuler,max(fabs(y_true[0]-y_beuler[0]),
					fabs(y_true[1]-y_beuler[1])));

    cout << "t:" <<setw(12)<< t1;
    cout << "  E:" <<setw(12)<< y_euler[0]  << " " <<setw(12)<< y_euler[1];
    cout << "  B:" <<setw(12)<< y_beuler[0] << " " <<setw(12)<< y_beuler[1];
    cout << "  T:" <<setw(12)<< y_true[0]   << " " <<setw(12)<< y_true[1];
    cout << endl;

    if(j==nTimes) break;
    j++;

    t0 = t1;
    t1 = t0+h;

    euler(y_euler,fcn,t0,t1,nSteps);

    state st = backEuler(y_beuler,fcn,dfcn,t0,t1,nSteps,tol,maxIter);
    if(st != SUCCESS) {
      cout << "ERROR: Newton failed to converge!" << endl;
      return 1;
    }
  }

  cout << "Final max log error:";
  cout << "  E: " << log(error_euler) << "  B: " << log(error_beuler)
    << "  log(1/h): " << -log(h/nSteps) << endl;
  cout << "Final max error:    ";
  cout << "  E: " << error_euler << "  B: " << error_beuler
    << "  log(1/h): " << h/nSteps << endl;
}
