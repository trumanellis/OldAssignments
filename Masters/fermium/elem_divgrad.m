function cF=elem_divgrad(QP,wgts2d,zonevel,GradVBasis,invJ,detJ,ndofpZ)

% Compute corner forces
stiffmat=zeros(ndofpZ);

for n=1:QP
    % Grad p Corner Force
    invJGradV=invJ(:,:,n)*GradVBasis(:,:,n)';
    stiffmat=stiffmat+wgts2d(n)*invJGradV'*invJGradV*detJ(n);
end
cF=-zonevel*stiffmat;