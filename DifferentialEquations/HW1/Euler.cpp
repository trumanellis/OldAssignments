// Euler's Method
// Written by: Truman Ellis
// Numerical Treatment of Differential Equations
// Spring 2011

#include "Python.h"
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <math.h>

#define EULERTEST

using namespace std;

class DifferentialEquation
{
    public:
        DifferentialEquation(){}
        virtual void eval (double t, vector<double> y, vector<double> &f) = 0;
        virtual void exact(double t0, double t1, vector<double> &y) = 0;
};


class DE0 : public DifferentialEquation
{
    // f(t,y) = -1/y, d=1
    public:
        void eval (double t, vector<double> y, vector<double> &f)
        {
            f[0] = -1/y[0];
        }

        void exact(double t0, double t1, vector<double> &y)
        {
            y[0] = sqrt(1-2*t1);
        }
};


class DE1 : public DifferentialEquation
{
    // f(t,y) = (-2y[0],-3y[1]), d=2
    public:
        void eval (double t, vector<double> y, vector<double> &f)
        {
            f[0] = -2*y[0];
            f[1] = -3*y[1];
        }

        void exact(double t0, double t1, vector<double> &y)
        {
            double dt = t1-t0;
            y[0] = y[0]*exp(-2*dt);
            y[1] = y[1]*exp(-3*dt);
        }
};


class DE2 : public DifferentialEquation
{
    // f(t,y) = (1/2, -t), d=2
    public:
        void eval (double t, vector<double> y, vector<double> &f)
        {
            f[0] = 0.5;
            f[1] = -t;
        }

        void exact(double t0, double t1, vector<double> &y)
        {
            y[0] = y[0] + 0.5*(t1-t0);
            y[1] = y[1] - 0.5*(t1*t1-t0*t0);
        }
};


class DE3 : public DifferentialEquation
{
    // f(t,y) = (y[0]^2t,y[0]*y[1]), d=1
    public:
        void eval (double t, vector<double> y, vector<double> &f)
        {
            f[0] = -y[0]*y[0]*t;
        }

        void exact(double t0, double t1, vector<double> &y)
        {
            y[0] = y[0]/(1 + 0.5*(t1*t1-t0*t0)*y[0]);
        }
};


class DE4 : public DifferentialEquation
{
    // f(t,y) = (y[0]^2t,y[0]*y[1]), d=2
    public:
        void eval (double t, vector<double> y, vector<double> &f)
        {
            f[0] = -y[0];
            f[1] = y[0];
        }

        void exact(double t0, double t1, vector<double> &y)
        {
            double e = exp(t0-t1);
            double yy = y[0];
            y[0] = yy*e;
            y[1] = yy*(1-e)+y[1];
        }
};


class EulerSolution
{
    private:
        shared_ptr<DifferentialEquation> de;
        double tStart;
        double tStop;
        int Nsteps;
        double h;
        vector<double> y0;
        int dim;
        vector< vector<double> > solution;
        vector<double> time;
        vector< vector<double> > EXsolution;
        double error;

    public:
        EulerSolution(shared_ptr<DifferentialEquation> de_, double tStart_,
                    double tStop_, int Nsteps_, vector<double> y0_): de(de_),
                    tStart(tStart_), tStop(tStop_), Nsteps(Nsteps_), y0(y0_)
        {
            h = (tStop-tStart)/double(Nsteps-1);
            dim = y0.size();
            solution.push_back(y0);
            time.push_back(tStart);
        }

        double solve()
        {
            vector<double> f(2);
            for (int n=0; n < Nsteps; ++n){
                double t = tStart + n*h;
                time.push_back(t+h);
                de->eval(t, solution[n], f);
                for (int d = 0; d < dim; ++d){
                    f[d] = solution[n][d] + f[d]*h;
                }
                solution.push_back(f);
            }
            error = exact();
            return error;
        }

        double exact()
        {
            double t0 = tStart;
            double t1 = tStart;
            vector<double> EXy = y0;
            error = 0;
            for (int n=0; n < Nsteps; ++n){
                de->exact(t0, t1, EXy);
                EXsolution.push_back(EXy);
                for (int d=0; d< dim; ++d){
                    error = max(error, fabs(EXy[d]-solution[n][d]));
                }
                t0 = t1;
                t1 += h;
            }
            return error;
        }

        void write()
        {
            ofstream output("output.txt");
            output << "#time        Euler: y[0]  Exact: y[0]";
            cout << "#time        Euler: y[0]  Exact: y[0]";
            for (int d = 1; d < dim; ++d)
            {
                cout << "  Euler: y["<<d<<"]  Exact: y["<<d<<"]";
                output << "  Euler: y["<<d<<"]  Exact: y["<<d<<"]";
            }
            cout << endl;
            output << endl;
            cout.setf(ios::fixed,ios::floatfield);
            output.setf(ios::fixed,ios::floatfield);
            cout.precision(10);
            output.precision(10);
            for(int i=0; i < Nsteps; i++){
                output << time[i] << " ";
                cout << time[i] << " ";
                for(int d=0; d < dim; ++d){
                    output << solution[i][d] << " " << EXsolution[i][d] << " ";
                    cout << solution[i][d] << " " << EXsolution[i][d] << " ";
                }
                output << endl;
                cout << endl;
            }
            output.close();
        }
};

#ifdef EULERTEST

int main()
{
    int probNum;
    int probDim;
    int Nsteps;
    double tStart;
    double tStop;
    vector<double> y0;
    vector<double> f;
    shared_ptr<DifferentialEquation> de;

    cout << "Which test problem (0, 1, 2, 3, 4): ";
    cin >> probNum;
    switch(probNum)
    {
        case 0:
            de = shared_ptr<DifferentialEquation> (new DE0);
            break;
        case 1:
            de = shared_ptr<DifferentialEquation> (new DE1);
            break;
        case 2:
            de = shared_ptr<DifferentialEquation> (new DE2);
            break;
        case 3:
            de = shared_ptr<DifferentialEquation> (new DE3);
            break;
        case 4:
            de = shared_ptr<DifferentialEquation> (new DE4);
            break;
        default:
            return 1;
    }
    cout << "Please enter problem dimension: ";
    cin >> probDim;
    cout << "Please enter start time: ";
    cin >> tStart;
    cout << "Please enter stop time: ";
    cin >> tStop;
    cout << "Please enter number of steps: ";
    cin >> Nsteps;
    cout << "Please enter "<<probDim<<" initial conditions: ";
    for (int d = 0; d < probDim; ++d){
        double ic;
        cin >> ic;
        y0.push_back(ic);
    }

    EulerSolution euler(de, tStart, tStop, Nsteps, y0);
    double error = euler.solve();
    euler.write();
    cout << "error: " << error << endl;

    vector<double> errorArray;
    vector<double> refArray;

    ofstream output("errors.txt");
    output.setf(ios::fixed,ios::floatfield);
    output.precision(10);
    output << Nsteps << " " << error << endl;
    refArray.push_back(Nsteps);
    errorArray.push_back(error);
    for (int ref = 1; ref <= 5; ++ref){
        Nsteps = 2*Nsteps;
        EulerSolution eulerRef(de, tStart, tStop, Nsteps, y0);
        refArray.push_back(Nsteps);
        error = eulerRef.solve();
        errorArray.push_back(error);
        output << Nsteps << " " << error << endl;
    }
    output.close();

    Py_Initialize();
    FILE *fp = fopen ("EulerPlot.py", "r");
    PyRun_SimpleFile(fp, "EulerPlot.py");
    Py_Exit(0);


    return 0;
}

#endif