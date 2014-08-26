function [w]=Q0r4Basis(x,y,index)

w=zeros(16,1);
d=[0 .25 .5 .75 1];

ix=find(d >= x);
ix=ix(1)-1;
iy=find(d >= y);
iy=iy(1)-1;

if x==0
    ix=1;
end
if y == 0;
    iy=1;
end

w((iy-1)*4+ix)=1;