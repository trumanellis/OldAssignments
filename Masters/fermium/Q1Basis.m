function [w]=Q1Basis(x,y,index)

w=zeros(4,1);

w(1)=(1-x)*(1-y);
w(2)=x*(1 - y);
w(3)=x*y;
w(4)=(1-x)*y;