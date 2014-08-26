function [massM,cF,psi0,mu]=elem_J(nodes,QP,wgts2d,refMQ2,GradVBasis,invJ,detJ,...
    massUpdate,forceUpdate,ndofpZ,Pquad,Dquad,zonevel,massZ,gamma,len,dt,qlin,qquad,Qfrac,hgfrac)
% Higher order elements
psi0=0;
mu=0;

% Assemble local mass matrix
massM=zeros(ndofpZ,ndofpZ);
if massUpdate
    for n=1:QP
        massM=massM+wgts2d(n)*Dquad(n)*refMQ2(:,:,n);
    end
end

if forceUpdate
    % Compute corner forces
    cF=zeros(2,ndofpZ);

    % Tensor artificial viscosity coefficient
    if Qfrac == 1 || hgfrac > 0
        stiffmat=zeros(ndofpZ,ndofpZ);
        SMlocal=zeros(ndofpZ,ndofpZ,QP);
        mu=zeros(QP,1);
        temp3=zeros(QP,1);
        cs=zeros(QP,1);
    end
    for n=1:QP
        % Grad p Corner Force
        invJGradV=invJ(:,:,n)*GradVBasis(:,:,n)';
        cF=cF+wgts2d(n)*Pquad(n)*invJGradV;

        % Stiffness matrix for tensor Q force
        if Qfrac == 1 || hgfrac > 0
            dV=zonevel*invJGradV';
            DivV=dV(1,1)+dV(2,2);
            CurlV=dV(2,1)-dV(1,2);
            
            temp1 = 1.0;
            temp2 = 1.0;
            if DivV > 0
                temp1 = exp(-10*DivV);
                temp2 = exp(-100*DivV);
            end
            if DivV == 0 && CurlV == 0
                temp3(n)=1;
            else
                temp3(n)=abs(DivV)/(abs(DivV)+abs(CurlV));
            end
            cs(n)=sqrt(gamma*Pquad(n)/Dquad(n));
            SMlocal(:,:,n)=wgts2d(n)*invJGradV'*invJGradV;
            stiffmat=stiffmat+SMlocal(:,:,n);
            % The viscosity coefficient
            mu(n)=Dquad(n)*len*(temp2*qquad*abs(DivV)*len+temp1*qlin*cs(n));
        end
    end
    if Qfrac == 1 
        if ndofpZ == 9
            smoothForce=zonevel*stiffmat;
            gp=.25*det(JacobianQ2Basis(nodes,[.5 .5]))*sqrt(gamma*mean(Pquad)/mean(Dquad))/(len*len);
            fp=sqrt(smoothForce(1,9)^2+smoothForce(2,9)^2);
            psi0=1-exp(-fp/(0.005*gp));
        else
            psi0=1;
        end
        stiffmat=zeros(ndofpZ,ndofpZ);
        for n=1:QP
            stiffmat=stiffmat+psi0*mu(n)*SMlocal(:,:,n);
        end
        % Add artificial viscosity
        cF=cF-zonevel*stiffmat;
    end
    if hgfrac > 0
        stiffmat=zeros(ndofpZ,ndofpZ);
        for n=1:QP
            stiffmat=stiffmat+Dquad(n)*len*temp3(n)*cs(n)*SMlocal(:,:,n);
        end
        % Add HG Filter
        cF=cF-hgfrac*zonevel*stiffmat;
    end
else
    cF=[];
end