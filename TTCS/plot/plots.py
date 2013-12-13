#!/bin/python

from pylab import *

# Set to use TeX
rc('text', usetex=True)

# Exponential Plot
figure(1)
x = linspace(0,5,50)
y=exp(x)

plot(x,y,'-o')
xlabel(r'$x$', fontsize=18)
legend((r'$\exp(x)$',))
text(0.5, 15, r'\Large Matplotlib is {\Huge $\mathbf{e}^x$}\\better than Gnuplot')

# Subplots
fig=figure(2)
x  = linspace(-pi, pi, 200)
y1 = sin(x)
y2 = sin(2*x)
y3 = sin(4*x)
y4 = sin(8*x)

subplot(2,2,1, title='Subplots')
plot(x,y1)
xlabel('x')
ylabel('sin(x)')
axis('tight')
grid(True)
annotate('Something\nimportant', xy=(0, 0), xytext=(-2.5, .25),
            arrowprops=dict(width=1, headwidth=8, ec='k', fc='k'),
            )



subplot(2,2,2)
plot(x,y2)
xlabel('x')
ylabel('sin(2x)')
axis('tight')
grid(True)

subplot(2,2,3)
plot(x,y3)
xlabel('x')
ylabel('sin(4x)')
axis('tight')
grid(True)

subplot(2,2,4)
plot(x,y4)
xlabel('x')
ylabel('sin(8x)')
axis('tight')
grid(True)

# Plotting Data from a File
data = loadtxt('data.dat')
show()
