# Statistical Mechanics Homework 5
# Spring 2011
# Written by Truman Ellis

from pylab import *
close('all')

def gauss1(x, y):
    dg = zeros(2)
    d2g = zeros((2,2))
    g = -200*exp(-(x-1)*(x-1)-10*y*y)
    dg[0] = -2*(x-1)*g
    dg[1] = -20*y*g
    d2g[0,0] = (-2+4*(x-1)*(x-1))*g
    d2g[1,1] = (-20 + 400*y*y)*g
    d2g[0,1] = 40*(x-1)*y*g
    d2g[1,0] = d2g[0,1]
    return (g, dg, d2g)

def gauss2(x, y):
    dg = zeros(2)
    d2g = zeros((2,2))
    g = -100*exp(-x*x-10*(y-.5)*(y-.5))
    dg[0] = -2*x*g
    dg[1] = -20*(y-.5)*g
    d2g[0,0] = (-2+4*x*x)*g
    d2g[1,1] = (-20 + 400*(y-.5)*(y-.5))*g
    d2g[0,1] = 40*x*(y-.5)*g
    d2g[1,0] = d2g[0,1]
    return (g, dg, d2g)

def gauss3(x, y):
    dg = zeros(2)
    d2g = zeros((2,2))
    g = -170*exp(-6.5*(x+.5)*(x+.5)+11*(x+.5)*(y-1.5)-6.5*(y-1.5)*(y-1.5))
    dg[0] = (-13*(x+.5)+11*(y-1.5))*g
    dg[1] = (-13*(y-1.5)+11*(x+.5))*g
    d2g[0,0] = ((-13*(x+.5)+11*(y-1.5))**2-13)*g
    d2g[1,1] = ((-13*(y-1.5)+11*(x+.5))**2-13)*g
    d2g[0,1] = (11+(-13*(x+.5)+11*(y-1.5))*(-13*(y-1.5)+11*(x+.5)))*g
    d2g[1,0] = d2g[0,1]
    return (g, dg, d2g)

def gauss4(x, y):
    dg = zeros(2)
    d2g = zeros((2,2))
    g = 15*exp(.7*(x+1)*(x+1)+.6*(x+1)*(y-1)+.7*(y-1)*(y-1))
    dg[0] = (1.4*(x+1)+0.6*(y-1))*g
    dg[1] = (1.4*(y-1)+0.6*(x+1))*g
    d2g[0,0] = ((1.4*(x+1)+.6*(y-1))**2+1.4)*g
    d2g[1,1] = ((1.4*(y-1)+.6*(x+1))**2+1.4)*g
    d2g[0,1] = (.6+(1.4*(x+1)+.6*(y-1))*(1.4*(y-1)+.6*(x+1)))*g
    d2g[1,0] = d2g[0,1]
    return (g, dg, d2g)

def muller((x, y)):
    (g1, dg1, d2g1) = gauss1(x, y)
    (g2, dg2, d2g2) = gauss2(x, y)
    (g3, dg3, d2g3) = gauss3(x, y)
    (g4, dg4, d2g4) = gauss4(x, y)
    return g1 + g2 + g3 + g4

def dmuller((x, y)):
    (g1, dg1, d2g1) = gauss1(x, y)
    (g2, dg2, d2g2) = gauss2(x, y)
    (g3, dg3, d2g3) = gauss3(x, y)
    (g4, dg4, d2g4) = gauss4(x, y)
    return dg1 + dg2 + dg3 + dg4

def VerletIntegrator(Xi, Vi, m, dt):
    dUi = dmuller(Xi)
    Xip1 = Xi + Vi*dt - 0.5*dt*dt/m*dUi
    dUip1 = dmuller(Xip1)
    Vip1 = Vi - 0.5*dt/m*(dUi + dUip1)

    return (Xip1, Vip1)

# Initialize particle at lowest energy point
m = 1.
E = -39
X0 = (-0.55822362835932715, 1.4417258468422827)
U0 = muller(X0)
# Calculate required velocity to produce E0
V2 = 2/m*(E-U0)
# Randomly set Vx and Vy as fractions of V2
Vx2 = rand()*V2
Vy2 = V2 - Vx2
# Also randomize signs of velocities
V0 = (sign(.5-rand())*sqrt(Vx2), sign(.5-rand())*sqrt(Vy2))

N = 200000
dt = 4e-3
X = zeros((N,2))
V = zeros((N,2))
Xave = zeros(N)
Yave = zeros(N)
XYave = zeros(N)
X[0] = X0
V[0] = V0
Xave[0] = X0[0]
Yave[0] = X0[1]
XYave[0] = X0[0]*X0[1]
for n in range(1,N):
    (X[n], V[n]) = VerletIntegrator(X[n-1], V[n-1], m, dt) 
    Xave[n] = average(X[0:n+1,0])
    Yave[n] = average(X[0:n+1,1])
    XYave[n] = average(X[0:n+1,0]*X[0:n+1,1])

xmin = min(X[:,0])-.01
xmax = max(X[:,0])+.01
ymin = min(X[:,1])-.01
ymax = max(X[:,1])+.01
# Compute background pseudocolor plot of Muller potential
nx = 100
ny = 101
mu = zeros( (nx,ny) )
x = linspace(xmin, xmax, 100)
y = linspace(ymin, ymax, 101)
for i in range(0,nx):
    for j in range(0,ny):
        mu[i,j] = muller((x[i],y[j]))

figure(1)
pcolor(x,y,mu)
hold(True)
plot(X[:,0],X[:,1],'k', linewidth=1)
axis([xmin,xmax,ymin,ymax])

figure(2)
hold(True)
plot(Xave)
plot(Yave)
plot(XYave)
legend(('Average x = %.4f'%(Xave[-1],),'Average y = %.4f'%(Yave[-1],),'Average \
    xy = %.4f'%(XYave[-1],)),loc='best')
show()

#Calculate energy error
Ef = muller(X[-1]) + .5*m*sum(V[-1]**2)
EnergyError = abs((E - Ef)/E)
print 'Energy error: %.6f'%(EnergyError,)
