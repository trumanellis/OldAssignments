% Assemble Mass Matrix and Force Vector
if MassUpdate
    MMQuad=sparse(ndofpD,ndofpD);
end
forceN=zeros(2,ndofpD);
for i=1:NZ
    if StrongMass
        [massM,cornerForce(:,:,i),psi0Z(i),mu(:,i)]=elem_J(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,...
            QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),MassUpdate,true,...
            ndofpZ,QuadPressure(:,i),QuadDensity(:,i),OLDvelocity(:,Quadmap(:,i)),...
            massZ(i),gamma,len(i),dt,qlin,qquad,Qfrac,hgfrac);
    else
        switch func2str(VBasis)
            case 'Q1Basis'
                [massM,cornerForce(:,:,i),mu(:,i)]=elemQ1(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,...
                    QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),MassUpdate,true,...
                    QuadPressure(:,i),QuadDensity(:,i),OLDvelocity(:,Quadmap(:,i)),...
                    massZ(i),gamma,len(i),dt,qlin,qquad,Qfrac,hgfrac);
			otherwise
                [massM,cornerForce(:,:,i),psi0Z(i),mu(:,i)]=elemQ2(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,...
                    QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),MassUpdate,true,...
                    ndofpZ,QuadPressure(:,i),QuadDensity(:,i),OLDvelocity(:,Quadmap(:,i)),...
                    massZ(i),gamma,len(i),dt,qlin,qquad,Qfrac,hgfrac);
        end
    end
    if MassUpdate
        MMQuad(Quadmap(:,i),Quadmap(:,i))=MMQuad(Quadmap(:,i),Quadmap(:,i))+massM;
    end
    forceN(:,Quadmap(:,i))=forceN(:,Quadmap(:,i))+cornerForce(:,:,i);
end
if ~FullMassMatrixSolve && MassUpdate
    massN=full(sum(MMQuad));
end
MMX=MMQuad;
MMY=MMQuad;