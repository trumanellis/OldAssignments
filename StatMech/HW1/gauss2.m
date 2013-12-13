function [g,dg,d2g] = gauss2(x,y);
% computing the second gaussian in the Mueller potential
g = -100*exp(-x*x-10*(y-0.5)*(y-0.5));
dg(1)=-2*x*g;
dg(2)=-20*(y-0.5)*g;
d2g(1,1)=(-2+4*x*x)*g;
d2g(2,2)=(-20+400*(y-0.5)*(y-0.5))*g;
d2g(1,2)=40*x*(y-0.5)*g;
d2g(2,1)=d2g(1,2);