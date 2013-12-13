# -*- coding: utf-8 -*-
# Newton's method with backward Euler
# Written by Truman Ellis

from pylab import *

def DE(t, y):
    return 3*t**2*y
    #return -t*y**2

def dDE(t, y):
    return 3*t**2
    #return -2*t*y

def Exact(t, y0):
    return y0*exp(t**3)
    #return 1./(1./y0 + 0.5*t**2)

Nsteps = 21
y0 = 1;
y=zeros(Nsteps)
yr=zeros(Nsteps)
y_exact=zeros(Nsteps)
y_exact[0]=y0
t=zeros(Nsteps)
y[0]=y0
error = zeros(5)
refinement = zeros(5)
refinement[0] = Nsteps

maxIter = 10;
for n in range(0,Nsteps-1):
    h = 1./Nsteps
    t[n+1] = t[n]+h
    wi = y[n] + h*DE(t[n],y[n])
    for i in range(0,maxIter):
        wip1 = wi - (wi - h*DE(t[n+1],wi)-y[n])/(1 - h*dDE(t[n+1],wi))
        #print str(wi)+" "+str(wip1)
        #print (wi - h*DE(t[n+1],wi)-y[n])
        #print (1-h*dDE(t[n+1],wi))
        if (abs(wip1 - wi) < 1e-25):
            print str(i)+" "+str(abs(wip1-wi))
            break
        wi = wip1
    y[n+1] = wip1
    y_exact[n+1] = Exact(t[n+1],y0)
    error[0] = max(error[0],abs(y[n+1]-y_exact[n+1]))

plot(t,y,'o')
plot(t,y_exact)

#for r in range(1,5):
#    Nsteps = 2*Nsteps
#    refinement[r] = Nsteps
#    y=zeros(Nsteps)
#    y_exact=zeros(Nsteps)
#    y_exact[0]=y0
#    t=zeros(Nsteps)
#    y[0]=y0
#    for n in range(0,Nsteps-1):
#        h = 1./Nsteps
#        t[n+1] = t[n]+h
#        if n > 2:
#            #y[n+1] = y[n]+h*(3./2.*DE(t[n],y[n])-1./2.*DE(t[n-1],y[n-1]))
#            #y[n+1] = y[n]+h*(23./12.*DE(t[n],y[n])-4./3.*DE(t[n-1],y[n-1])+5./12.*DE(t[n-2],y[n-2]))
#            y[n+1] = y[n]+h*(55./24.*DE(t[n],y[n])-59./24.*DE(t[n-1],y[n-1])\
#                +37./24.*DE(t[n-2],y[n-2])-3./8.*DE(t[n-3],y[n-3]))
#        else:
#            #h = h*h
#            y[n+1] = Exact(t[n+1],y0)
#        #y[n+1] = y[n]+h*DE(t[n],y[n])
#        y_exact[n+1] = Exact(t[n+1],y0)
#        error[r] = max(error[r],abs(y[n+1]-y_exact[n+1]))
#
#figure()
#loglog(refinement,error)
#(m,b) = polyfit(log(refinement),log(error),1)
#xlabel('Degrees of Freedom')
#ylabel('Error')
#legend(('rate=%.2f' % (-m,),), loc='best')
show()
