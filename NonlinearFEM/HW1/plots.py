from pylab import *
close('all')


c = 10.
M = -5.
Q = -5.
EI = 100.

A = -c/(6)-M/(2)-Q
B = -c/(24)-M/(6)-Q/(2)-A

d1 = 1./24*(3*c+12*M+8*Q)
d2 = 1./6*(-c-6*M-3*Q)

x = linspace(0,1,100)
uh = d1*(x-1)**2*(2*x+1)+d2*x*(x-1)**2
uhx = d1 * 6*x*(x-1) + d2 * (1-4*x+3*x**2)
uhxx = d1 * (12*x-6) + d2 * (6*x-4)
ue = c/(24)*x**4+M/(6)*x**3+Q/(2)*x**2+A*x+B
uex = c/(6)*x**3+M/(2)*x**2+Q*x+A
uexx = c/(2)*x**2+M*x+Q

x1 = 0.
x2 = 0.5
h = x2 - x1
xl = linspace(-1, 1,100)
N1l = .25*(xl-1)**2*(xl+2)
N2l = h/8*(xl+1)*(xl-1)**2
N3l = -.25*(xl-2)*(xl+1)**2
N4l = h/8*(xl+1)**2*(xl-1)
xg = linspace(x1, x2,100)
N1g = -(xg-x2)**2*(-h+2*(x1-xg))/h**3
N2g = (xg-x1)*(xg-x2)**2/h**2
N3g = (xg-x1)**2*(h+2*(x2-xg))/h**3
N4g = (xg-x1)**2*(xg-x2)/h**2

figure()
plot(xl, N1l)
plot(xl, N2l)
plot(xl, N3l)
plot(xl, N4l)
legend( (r'$N_1(\xi)$',r'$N_2(\xi)$',r'$N_3(\xi)$',r'$N_4(\xi)$'), loc='best')
title('Local')
xlabel(r'$\xi$', fontsize=20)

figure()
plot(xg, N1g)
plot(xg, N2g)
plot(xg, N3g)
plot(xg, N4g)
legend( ('$N_1(x)$','$N_2(x)$','$N_3(x)$','$N_4(x)$'), loc='best')
title('Global')
xlabel(r'$x$', fontsize=20)

#figure()
#plot(x,ue)
#plot(x,uh)
#ylabel('EI u')
#xlabel('x')
#legend(('Exact','Approximate'), loc='best')
#
#figure()
#plot(x,uex)
#plot(x,uhx)
#ylabel('EI u,x')
#xlabel('x')
#legend(('Exact','Approximate'), loc='best')
#
#figure()
#plot(x,uexx)
#plot(x,uhxx)
#Bp1 = 0.5*(1 - 1/sqrt(3))
#Bp2 = 0.5*(1 + 1/sqrt(3))
#y = [-5.833, -5.833]
##y = [.0218, .3125]
#plot([Bp1, Bp2], y, 'rx', markersize=10, markeredgewidth=3)
#ylabel('EI u,xx')
#xlabel('x')
#legend(('Exact','Approximate','Barlow curvature points'), loc='best')

show()
