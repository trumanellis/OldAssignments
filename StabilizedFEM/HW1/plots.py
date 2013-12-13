from pylab import *
from mpl_toolkits.mplot3d import Axes3D

close('all')

ymid = 0.5
x = linspace(0, 1, 201)

alphas = [0.1, 1, 2, 5, 10, 100]
count = 0
for alphae in alphas:
  count = count + 1

  def g(x,y):
    C1 = (1-exp(-2*alphae*(1-y)))/(1-exp(-2*alphae))
    C2 = (exp(2*alphae*y)-1)/(1-exp(-2*alphae))
    if x <= y:
      return C1*(1-exp(-2*alphae*x))
    else:
      return C2*(exp(-2*alphae*x)-exp(-2*alphae))
  g = vectorize(g)
  def b(x):
    return -(1-exp(2*alphae*x))/(1-exp(2*alphae))+x
  b = vectorize(b)

  figure(1)
  plot(x, g(x,ymid), label=r'$\alpha^e='+str(alphae)+'$')

  X, Y = meshgrid(x, x)
  G = g(X, Y)
  fig = figure(count+2)
  ax = fig.gca(projection='3d')
  ax.view_init(20,-120)
  surf = ax.plot_surface(X,Y,G)
  title(r'$\alpha^e='+str(alphae)+'$')

  figure(2)
  plot(x, b(x), label=r'$\alpha^e='+str(alphae)+'$')
figure(1)
legend()

figure(2)
legend(loc='Best')

show()
