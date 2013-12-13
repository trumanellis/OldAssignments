#ifndef __VERLETINTEGRATOR_HPP
#define __VERLETINTEGRATOR_HPP

#include "builtin.hpp"
#include "random.hpp"
#include "math.hpp"
#include "time.hpp"
#include "VerletIntegrator.hpp"

using namespace __shedskin__;
namespace __VerletIntegrator__ {

extern str *const_0, *const_1, *const_2, *const_3, *const_4;

using __math__::exp;
using __math__::sqrt;
using __random__::random;



extern __ss_int N, __4, __5, n;
extern double E, U0, V2, Vx2, Vy2, dt, m, vx0, vy0, x0, y0;
extern list<double> *XYave, *Xave, *Yave, *dUn, *dUnp1, *vx, *vy, *x, *xy, *y;
extern str *__name__;
extern file *f;

double sign(double x);
double gauss1(double x, double y);
list<double> *dgauss1(double x, double y);
double gauss2(double x, double y);
list<double> *dgauss2(double x, double y);
double gauss3(double x, double y);
list<double> *dgauss3(double x, double y);
double gauss4(double x, double y);
list<double> *dgauss4(double x, double y);
double muller(tuple2<double, double> *__0);
list<double> *dmuller(tuple2<double, double> *__2);

} // module namespace
#endif
