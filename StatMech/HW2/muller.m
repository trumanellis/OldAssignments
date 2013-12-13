function [mu,dmu,d2mu] = mueller(x,y);
% compute the mueller potential and it first and second derivatives
[g1,dg1,d2g1] = gauss1(x,y);
[g2,dg2,d2g2] = gauss2(x,y);
[g3,dg3,d2g3] = gauss3(x,y);
[g4,dg4,d2g4] = gauss4(x,y);
mu = g1 + g2 + g3 + g4;
dmu = dg1 + dg2 + dg3 + dg4;
d2mu = d2g1 + d2g2 + d2g3 + d2g4;
