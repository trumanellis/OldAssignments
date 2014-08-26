function [cF]=WeightedCornerForce(nodes,QP,wgts2d,GradVBasis,PBasis,invJ,detJ,...
    ndofpZ,npdof,Pquad,Dquad,zonevel,gamma,len,qlin,qquad,Qfrac)


% Compute corner forces
cF=zeros(2,ndofpZ,npdof);

% Tensor artificial viscosity coefficient
if Qfrac == 1
    stiffmat=zeros(ndofpZ,ndofpZ,npdof);
    SMlocal=zeros(ndofpZ,ndofpZ,QP);
    mu=zeros(1,QP);
end
for n=1:QP
    % Grad p Corner Force
    invJGradV=invJ(:,:,n)*GradVBasis(:,:,n)';
    for j=1:npdof
        cF(:,:,j)=cF(:,:,j)+wgts2d(n)*Pquad(n)*invJGradV*PBasis(j,n)*detJ(n);
    end
    
    % Stiffness matrix for tensor Q force
    if Qfrac == 1
        dV=zonevel*invJGradV';
        DivV=dV(1,1)+dV(2,2);
        if DivV < 0
            % The viscosity coefficient
            SMlocal(:,:,n)=wgts2d(n)*invJGradV'*invJGradV*detJ(n);
            stiffmat(:,:,1)=stiffmat(:,:,1)+SMlocal(:,:,n);
            for j=1:npdof
                mu(n)=Dquad(n)*len*(-qquad*DivV*len+qlin*sqrt(gamma*Pquad(n)/Dquad(n)));
            end
        end
    end
end
if Qfrac == 1
    smoothForce=zonevel*stiffmat(:,:,1);
    gp=.25*det(JacobianQ2Basis(nodes,[.5 .5]))*sqrt(gamma*mean(Pquad)/mean(Dquad))/(len*len);
    fp=sqrt(smoothForce(1,end)^2+smoothForce(2,end)^2);
    phi0=1-exp(-fp/(0.5*gp));
    stiffmat=zeros(ndofpZ,ndofpZ,npdof);
    for j=1:4
        for n=1:QP
            stiffmat(:,:,j)=stiffmat(:,:,j)+phi0*mu(n)*SMlocal(:,:,n)*PBasis(j,n);
        end
        % Add artificial viscosity
        cF(:,:,j)=cF(:,:,j)-zonevel*stiffmat(:,:,j);
    end
end