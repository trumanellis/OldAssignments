function [g,dg,d2g] = gauss4(x,y);
% computing the fourth gaussian in the Mueller potential
g = 15*exp(0.7*(x+1)*(x+1)+0.6*(x+1)*(y-1)+0.7*(y-1)*(y-1));
dg(1)=(1.4*(x+1)+0.6*(y-1))*g;
dg(2)=(1.4*(y-1)+0.6*(x+1))*g;
d2g(1,1)=(1.4+(1.4*(x+1)+0.6*(y-1))^2)*g;
d2g(2,2)=(1.4+(1.4*(y-1)+0.6*(x+1))^2)*g;
d2g(1,2)=(0.6+(1.4*(x+1)+0.6*(y-1))*(1.4*(y-1)+0.6*(x+1)))*g;
d2g(2,1)=d2g(1,2);