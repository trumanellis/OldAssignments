function Jmat=JacobianQ1bBasis(nodes,localPt)

X1=nodes(:,1)';
X2=nodes(:,2)';
X3=nodes(:,3)';
X4=nodes(:,4)';
X5=nodes(:,5)';
x=localPt(1);
y=localPt(2);

Jmat=[X2 - X1 + x*(-32*X5*y + 8*X1*y + 8*X2*y + 8*X3*y + 8*X4*y ...
    - 8*X1*y*y - 8*X2*y*y - 8*X3*y*y - 8*X4*y*y + 32*X5*y*y) ...
    + y*(-5*X2 - 5*X4 - 3*X1 - 3*X3 + 16*X5) + y*y*(-16*X5 + 4*X1 ...
    + 4*X2 + 4*X3 + 4*X4);
    X4 - X1 + x*(-5*X2 - 5*X4 - 3*X1 - 3*X3 + 16*X5 - 32*X5*y + 8*X1*y ...
    + 8*X2*y + 8*X3*y + 8*X4*y) + x*x*(-16*X5 + 4*X1 + 4*X2 + 4*X3 ...
    + 4*X4 - 8*X1*y - 8*X2*y - 8*X3*y - 8*X4*y + 32*X5*y)];