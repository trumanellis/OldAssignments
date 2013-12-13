function [g,dg,d2g] = gauss1(x,y);
% computing the first gassuain in the Mueller potential
g = -200*exp(-(x-1)*(x-1)-10*y*y);
dg(1)=-2*(x-1)*g;
dg(2)=-20*y*g;
d2g(1,1)=(-2+4*(x-1)*(x-1))*g;
d2g(2,2)=(-20+400*y*y)*g;
d2g(1,2)=40*(x-1)*y*g;
d2g(2,1)=d2g(1,2);