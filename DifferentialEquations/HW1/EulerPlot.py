# -*- coding: utf-8 -*-
# Euler's Method
# Written by: Truman Ellis
# Numerical Treatment of Differential Equations
# Spring 2011

from pylab import *

data=loadtxt('output.txt')

for d in range(1, (data.shape[1]+1)/2):
    figure(d)
    plot(data[:,0],data[:,2*d],'k-', label='Exact Solution')
    plot(data[:,0],data[:,2*d-1],'o', label="Euler's Method")
    legend(loc='best')
    xlabel('time')
    ylabel('y['+str(d-1)+']')

figure()
convergenceData = loadtxt('errors.txt')
(m,b) = polyfit(log(convergenceData[:,0]),log(convergenceData[:,1]),1)
loglog(convergenceData[:,0],convergenceData[:,1])
xlabel('Degrees of Freedom')
ylabel('Error')
legend(('rate=%.2f' % (-m,),), loc='best')
show()