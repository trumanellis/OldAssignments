# One-dimensional Finite Element Solver
# Solves -u(x)'' + cu(x)  = f(x,y)
# Written by: Truman Ellis
# Numerical Treatment of Differential Equations
# Spring 2011

from pylab import *

close('all')

# Define Initial Refinement
nx = 19

# Define problem domain
xmin = 0.
xmax = 1.

# Define boundary conditions
uL = -3.
uR = 2.

# Define scalar c
c = 0.

# Define exact solution for convergence tests
def exact(x):
    return sin(4*pi*x**3) + 5*x - 3

def dexact(x):
    return 4*pi*cos(4*pi*x**3)*3*x**2 + 5

def d2exact(x):
    return -16*pi**2*sin(4*pi*x**3)*9*x**4 + 4*pi*cos(4*pi*x**3)*6*x

# Define forcing function to allow for convergence tests
def f(x):
    return -d2exact(x) + c*exact(x)

# Set number of refinement steps for convergence test
nref = 4
error = zeros(nref)
derror = zeros(nref)
ndofs = zeros(nref)
hs = zeros(nref)
for r in range(0,nref):
    if r > 0:
        nx = 2*(nx+1)-1
    ndofs[r] = nx
    h = (xmax - xmin)/(nx + 1)
    hs[r] = h

    X = linspace(xmin, xmax, nx+2)
    x = X[1:-1]
    xh = (X[1:] + X[0:-1])/2.
    A = zeros((nx,nx))
    b = zeros(nx)
    u = zeros(nx)
    du = zeros(nx+1)
    u_ex = zeros(nx)
    du_ex = zeros(nx+1)

    for i in range(0,nx):
            A[i,i] = 2./h + c*2.*h/3.
            b[i] += h*f(x[i])
            if (i < nx-1):
                A[i,i+1] = -1./h + c*h/6.
            else:
                b[i] += -uR*(-1./h + c*h/6.)
            if (i > 0):
                A[i,i-1] = -1./h + c*h/6.
            else:
                b[i] += -uL*(-1./h + c*h/6.)

            # Calculate exact solution
            u_ex[i] = exact(x[i])

    # Solve for finite difference solution
    u = linalg.solve(A,b)
    du[1:-1] = (u[1:] - u[0:-1])/h
    du[0] = (u[0]-uL)/h
    du[-1] = (uR - u[-1])/h
    du_ex = dexact(xh)
    error[r] = sqrt(((u-u_ex)**2).sum()/nx)
    derror[r] = sqrt(((du-du_ex)**2).sum()/(nx+1))

    # Plotting
    if r == 2:
        U = zeros(nx+2)
        E = zeros(nx+2)
        U[1:-1] = u
        U[0] = uL
        U[-1] = uR
        E[1:-1] = u_ex
        E[0] = uL
        E[-1] = uR

        figure(1)
        plot(X, E, 'k')
        plot(X, U, 'o')
        legend(('Exact','FEM Solution'),loc='best')
        show()

        figure(2)
        plot(xh, du_ex, 'k')
        plot(xh, du, 'o')


figure()
loglog(hs,error,'-o')
(m,b) = polyfit(log(hs),log(error),1)
xlabel('h')
ylabel('$L_2$ Error Norm')
grid()
legend(('rate=%.2f' % (m,),), loc='best')

figure()
loglog(hs,derror,'-o')
(dm,b) = polyfit(log(hs),log(derror),1)
xlabel('h')
ylabel('Quasi-Optimal Error Norm')
grid()
legend(('rate=%.2f' % (dm,),), loc='best')
