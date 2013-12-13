from pylab import *

def NR():
  d = zeros((2,N))
  N1_NR = zeros(N)
  iterations_convergence = [0]

  unconverged = False
  for n in range(1, N):
    if unconverged:
      break
    di = d[:, n]
    Fext = array([F1ext(n), F2ext(n)])
    R0 = Fext - Fint(di)
    R = R0
    for i in range(0, 15):
      K = StiffnessMatrix(di, x)
      Dd = linalg.solve(K, R)
      di += Dd
      R = Fext - Fint(di)

      if (norm(R) < epsilon*norm(R0)):
        iterations_convergence.append(i+1)
        break
      if (i == 14):
        unconverged = True

    if (n < N-1):
      d[:, n+1] = di

    N1_NR[n] = N1(d[:,n])

  figure(1)
  max_load_step = len(iterations_convergence)
  plot(d[0,0:max_load_step], N1_NR[0:max_load_step], '-s', linewidth=2,
      markersize=5, label="Newton-Rhapson")
  figure(2)
  plot(range(0,max_load_step),iterations_convergence, label="Newton-Rhapson")
