function [massM,cF,mu]=elemQ1(nodes,QP,wgts2d,refMV,GradVBasis,invJ,detJ,...
    massUpdate,forceUpdate,Pquad,Dquad,zonevel,massZ,gamma,len,dt,qlin,qquad,Qfrac,hgfrac)

% Assemble local mass matrix
massM=zeros(4,4);
if massUpdate
    for n=1:QP
        massM=massM+wgts2d(n)*Dquad(n)*refMV(:,:,n)*detJ(n);
    end
end

mu=0;
if forceUpdate
    % Compute corner forces
    cF=zeros(2,4);

    % Tensor artificial viscosity coefficient
    if Qfrac == 1
        stiffmat=zeros(4,4);
        dV=zonevel*GradQ1Basis(.5,.5)*inv(JacobianQ1Basis(nodes,[.5 .5]))';
        DivV=dV(1,1)+dV(2,2);
        if DivV < 0
            mu=mean(Dquad)*len*(-qquad*DivV*len+qlin*sqrt(gamma*mean(Pquad)/mean(Dquad)));
        else
            mu=0;
        end
    end
    for n=1:QP
        % Grad p Corner Force
        invJGradV=invJ(:,:,n)*GradVBasis(:,:,n)';
        cF=cF+wgts2d(n)*Pquad(n)*invJGradV*detJ(n);

        % Stiffness matrix for tensor Q force
        if Qfrac == 1 %&& DivV < 0
            stiffmat=stiffmat+wgts2d(n)*invJGradV'*invJGradV*detJ(n);
        end
    end
    if Qfrac == 1 
        % Add artificial viscosity
        cF=cF-mu*zonevel*stiffmat;
    end
    if hgfrac > 0
        Vx=zonevel(1,:);
        Vy=zonevel(2,:);

        fx=0.25*(Vx(2)-Vx(3)+Vx(4)-Vx(1));
        fy=0.25*(Vy(2)-Vy(3)+Vy(4)-Vy(1));

        hgForces=[fx -fx fx -fx; fy -fy fy -fy];

        cF=cF+hgfrac/(8*dt)*massZ*hgForces;
    end
else
    cF=[];
end