function Jmat=JacobianMat(nodes,localPt)

X1=nodes(:,1)';
X2=nodes(:,2)';
X3=nodes(:,3)';
X4=nodes(:,4)';
x=localPt(1);
y=localPt(2);

Jmat=[X2 - X1 + y*(X1 + X3 - X2 - X4);
      X4 - X1 + x*(X1 + X3 - X2 - X4)];