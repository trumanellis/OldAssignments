function [w]=P1Basis(x,y,index)

w=zeros(3,1);

w(1)=1-x-y;
w(2)=x;
w(3)=y;