#include "VerletIntegrator.hpp"

namespace __VerletIntegrator__ {

str *const_0, *const_1, *const_2, *const_3, *const_4;

__ss_int N, __4, __5, n;
double E, U0, V2, Vx2, Vy2, dt, m, vx0, vy0, x0, y0;
list<double> *XYave, *Xave, *Yave, *dUn, *dUnp1, *vx, *vy, *x, *xy, *y;
str *__name__;
file *f;


double sign(double x) {
    
    if ((x>0.0)) {
        return 1.0;
    }
    else {
        return (-1.0);
    }
    return 0;
}

double gauss1(double x, double y) {
    double g;

    g = ((-200)*exp((((-(x-1))*(x-1))-((10.0*y)*y))));
    return g;
}

list<double> *dgauss1(double x, double y) {
    double g;
    list<double> *dg;

    dg = (new list<double>(2,0.0,0.0));
    g = gauss1(x, y);
    dg->__setitem__(0, (((-2)*(x-1))*g));
    dg->__setitem__(1, (((-20)*y)*g));
    return dg;
}

double gauss2(double x, double y) {
    double g;

    g = ((-100)*exp((((-x)*x)-((10*(y-0.5))*(y-0.5)))));
    return g;
}

list<double> *dgauss2(double x, double y) {
    double g;
    list<double> *dg;

    dg = (new list<double>(2,0.0,0.0));
    g = gauss2(x, y);
    dg->__setitem__(0, (((-2)*x)*g));
    dg->__setitem__(1, (((-20)*(y-0.5))*g));
    return dg;
}

double gauss3(double x, double y) {
    double g;

    g = ((-170)*exp((((((-6.5)*(x+0.5))*(x+0.5))+((11*(x+0.5))*(y-1.5)))-((6.5*(y-1.5))*(y-1.5)))));
    return g;
}

list<double> *dgauss3(double x, double y) {
    double g;
    list<double> *dg;

    dg = (new list<double>(2,0.0,0.0));
    g = gauss3(x, y);
    dg->__setitem__(0, ((((-13)*(x+0.5))+(11*(y-1.5)))*g));
    dg->__setitem__(1, ((((-13)*(y-1.5))+(11*(x+0.5)))*g));
    return dg;
}

double gauss4(double x, double y) {
    double g;

    g = (15*exp(((((0.7*(x+1))*(x+1))+((0.6*(x+1))*(y-1)))+((0.7*(y-1))*(y-1)))));
    return g;
}

list<double> *dgauss4(double x, double y) {
    double g;
    list<double> *dg;

    dg = (new list<double>(2,0.0,0.0));
    g = gauss4(x, y);
    dg->__setitem__(0, (((1.4*(x+1))+(0.6*(y-1)))*g));
    dg->__setitem__(1, (((1.4*(y-1))+(0.6*(x+1)))*g));
    return dg;
}

double muller(tuple2<double, double> *__0) {
    tuple2<double, double> *__1;
    double g1, g2, g3, g4, x, y;

    __1 = __0;
    x = __1->__getfirst__();
    y = __1->__getsecond__();
    g1 = gauss1(x, y);
    g2 = gauss2(x, y);
    g3 = gauss3(x, y);
    g4 = gauss4(x, y);
    return (((g1+g2)+g3)+g4);
}

list<double> *dmuller(tuple2<double, double> *__2) {
    tuple2<double, double> *__3;
    double x, y;
    list<double> *dg1, *dg2, *dg3, *dg4;

    __3 = __2;
    x = __3->__getfirst__();
    y = __3->__getsecond__();
    dg1 = dgauss1(x, y);
    dg2 = dgauss2(x, y);
    dg3 = dgauss3(x, y);
    dg4 = dgauss4(x, y);
    return (new list<double>(2,(((dg1->__getfast__(0)+dg2->__getfast__(0))+dg3->__getfast__(0))+dg4->__getfast__(0)),(((dg1->__getfast__(1)+dg2->__getfast__(1))+dg3->__getfast__(1))+dg4->__getfast__(1))));
}

void __init() {
    const_0 = new str("__main__");
    const_1 = new str("ParticlePath.txt");
    const_2 = __char_cache[119];;
    const_3 = __char_cache[32];;
    const_4 = __char_cache[10];;

    __name__ = new str("__main__");

    m = 1.0;
    E = 10.0;
    x0 = (-0.558223628359);
    y0 = 1.44172584684;
    if (__eq(__name__, const_0)) {
        U0 = muller((new tuple2<double, double>(2,x0,y0)));
        V2 = ((2.0/m)*(E-U0));
        Vx2 = (random()*V2);
        Vy2 = (V2-Vx2);
        vx0 = (sign((0.5-random()))*sqrt(Vx2));
        vy0 = (sign((0.5-random()))*sqrt(Vy2));
        N = 10000;
        dt = 0.004;
        x = (new list<double>(1,x0));
        y = (new list<double>(1,y0));
        xy = (new list<double>(1,(x0*y0)));
        vx = (new list<double>(1,vx0));
        vy = (new list<double>(1,vy0));
        Xave = (new list<double>(1,x0));
        Yave = (new list<double>(1,y0));
        XYave = (new list<double>(1,(x0*y0)));
        f = open(const_1, const_2);
        f->write(__add_strs(14, __str(x->__getfast__(0)), const_3, __str(y->__getfast__(0)), const_3, __str(vx->__getfast__(0)), const_3, __str(vy->__getfast__(0)), const_3, __str(Xave->__getfast__(0)), const_3, __str(Yave->__getfast__(0)), const_3, __str(XYave->__getfast__(0)), const_4));

        FAST_FOR(n,1,N,1,4,5)
            dUn = dmuller((new tuple2<double, double>(2,x->__getfast__((n-1)),y->__getfast__((n-1)))));
            x->append(((x->__getfast__((n-1))+(vx->__getfast__((n-1))*dt))-((((0.5*dt)*dt)/m)*dUn->__getfast__(0))));
            y->append(((y->__getfast__((n-1))+(vy->__getfast__((n-1))*dt))-((((0.5*dt)*dt)/m)*dUn->__getfast__(1))));
            dUnp1 = dmuller((new tuple2<double, double>(2,x->__getfast__(n),y->__getfast__(n))));
            vx->append((vx->__getfast__((n-1))-(((0.5*dt)/m)*(dUn->__getfast__(0)+dUnp1->__getfast__(0)))));
            vy->append((vy->__getfast__((n-1))-(((0.5*dt)/m)*(dUn->__getfast__(1)+dUnp1->__getfast__(1)))));
            Xave->append((__sum(x)/n));
            Yave->append((__sum(y)/n));
            xy->append((x->__getfast__(n)*y->__getfast__(n)));
            XYave->append((__sum(xy)/n));
            f->write(__add_strs(14, __str(x->__getfast__(n)), const_3, __str(y->__getfast__(n)), const_3, __str(vx->__getfast__(n)), const_3, __str(vy->__getfast__(n)), const_3, __str(Xave->__getfast__(n)), const_3, __str(Yave->__getfast__(n)), const_3, __str(XYave->__getfast__(n)), const_4));
        END_FOR

        f->close();
    }
}

} // module namespace

int main(int, char **) {
    __shedskin__::__init();
    __math__::__init();
    __time__::__init();
    __random__::__init();
    __shedskin__::__start(__VerletIntegrator__::__init);
}
