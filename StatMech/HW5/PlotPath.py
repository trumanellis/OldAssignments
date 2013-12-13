# Statistical Mechanics Homework 5
# Spring 2011
# Written by Truman Ellis

from pylab import *
from VerletIntegrator import muller, E, m
close('all')

pdata = loadtxt('ParticlePath.txt')
x = pdata[:,0]
y = pdata[:,1]
vx = pdata[:,2]
vy = pdata[:,3]
Xave = pdata[:,4]
Yave = pdata[:,5]
XYave = pdata[:,6]

xmin = min(x)-.01
xmax = max(x)+.01
ymin = min(y)-.01
ymax = max(y)+.01
# Compute background pseudocolor plot of Muller potential
nx = 100
ny = 101
mu = zeros( (nx,ny) )
xp = linspace(xmin, xmax, 100)
yp = linspace(ymin, ymax, 101)
for i in range(0,nx):
    for j in range(0,ny):
        mu[i,j] = muller((xp[i],yp[j]))

figure(1)
pcolor(xp,yp,mu)
hold(True)
plot(x,y,'k', linewidth=1)
axis([xmin,xmax,ymin,ymax])

figure(2)
hold(True)
semilogx(Xave)
semilogx(Yave)
semilogx(XYave)
legend(('Average  x = %.4f'%(Xave[-1],),'Average  y = %.4f'%(Yave[-1],), \
'Average xy = %.4f'%(XYave[-1],)),loc='best')
xlabel('Iteration Number')
show()

#Calculate energy error
Ef = muller((x[-1],y[-1])) + .5*m*sum(vx[-1]**2+vy[-1]**2)
EnergyError = abs((E - Ef)/E)
print 'Energy error: %.6f'%(EnergyError,)
