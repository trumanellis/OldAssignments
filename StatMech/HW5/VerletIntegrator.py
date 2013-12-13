# Statistical Mechanics Homework 5
# Spring 2011
# Written by Truman Ellis

from math import exp, sqrt
from random import random

m = 1.
E = -39.

# Initialize particle at lowest energy point
x0 = -0.55822362835932715
y0 = 1.4417258468422827

def sign(x):
    if x > 0.0:
        return 1.0
    else:
        return -1.0

def gauss1(x, y):
    g = -200*exp(-(x-1)*(x-1)-10.*y*y)
    return g
def dgauss1(x, y):
    dg = [0.,0.]
    g = gauss1(x,y)
    dg[0] = -2*(x-1)*g
    dg[1] = -20*y*g
    return dg

def gauss2(x, y):
    g = -100*exp(-x*x-10*(y-.5)*(y-.5))
    return g
def dgauss2(x, y):
    dg = [0.,0.]
    g = gauss2(x,y)
    dg[0] = -2*x*g
    dg[1] = -20*(y-.5)*g
    return dg

def gauss3(x, y):
    g = -170*exp(-6.5*(x+.5)*(x+.5)+11*(x+.5)*(y-1.5)-6.5*(y-1.5)*(y-1.5))
    return g
def dgauss3(x, y):
    dg = [0.,0.]
    g = gauss3(x,y)
    dg[0] = (-13*(x+.5)+11*(y-1.5))*g
    dg[1] = (-13*(y-1.5)+11*(x+.5))*g
    return dg

def gauss4(x, y):
    g = 15*exp(.7*(x+1)*(x+1)+.6*(x+1)*(y-1)+.7*(y-1)*(y-1))
    return g
def dgauss4(x, y):
    dg = [0.,0.]
    g = gauss4(x,y)
    dg[0] = (1.4*(x+1)+0.6*(y-1))*g
    dg[1] = (1.4*(y-1)+0.6*(x+1))*g
    return dg

def muller((x, y)):
    g1 = gauss1(x, y)
    g2 = gauss2(x, y)
    g3 = gauss3(x, y)
    g4 = gauss4(x, y)
    return g1 + g2 + g3 + g4
def dmuller((x, y)):
    dg1 = dgauss1(x, y)
    dg2 = dgauss2(x, y)
    dg3 = dgauss3(x, y)
    dg4 = dgauss4(x, y)
    return [dg1[0] + dg2[0] + dg3[0] + dg4[0], \
            dg1[1] + dg2[1] + dg3[1] + dg4[1]]

if __name__ == "__main__":
    U0 = muller((x0,y0))
    # Calculate required velocity to produce E0
    V2 = 2./m*(E-U0)
    # Randomly set Vx and Vy as fractions of V2
    Vx2 = random()*V2
    Vy2 = V2 - Vx2
    # Also randomize signs of velocities
    vx0 = sign(.5-random())*sqrt(Vx2)
    vy0 = sign(.5-random())*sqrt(Vy2)

    N = 20000
    dt = 4e-3
    x = [x0]
    y = [y0]
    xy = [x0*y0]
    vx = [vx0]
    vy = [vy0]
    Xave = [x0]
    Yave = [y0]
    XYave = [x0*y0]
    f = open('ParticlePath.txt', 'w')
    f.write(str(x[0])+' '+str(y[0])+' '+str(vx[0])+' '+str(vy[0])+' '\
        +str(Xave[0])+' '+str(Yave[0])+' '+str(XYave[0])+'\n')
    for n in range(1,N):
        dUn = dmuller((x[n-1],y[n-1]))
        x.append( x[n-1] + vx[n-1]*dt - 0.5*dt*dt/m*dUn[0] )
        y.append( y[n-1] + vy[n-1]*dt - 0.5*dt*dt/m*dUn[1] )
        dUnp1 = dmuller((x[n],y[n]))
        vx.append( vx[n-1] - 0.5*dt/m*(dUn[0] + dUnp1[0]) )
        vy.append( vy[n-1] - 0.5*dt/m*(dUn[1] + dUnp1[1]) )
        Xave.append( sum(x)/n )
        Yave.append( sum(y)/n )
        xy.append( x[n]*y[n] )
        XYave.append( sum(xy)/n )

        f.write(str(x[n])+' '+str(y[n])+' '+str(vx[n])+' '+str(vy[n])+' '\
            +str(Xave[n])+' '+str(Yave[n])+' '+str(XYave[n])+'\n')

    f.close()

