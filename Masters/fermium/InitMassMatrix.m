% Assemble Nodal Masses if Mass is not updated every cycle
if ~MassUpdate
    MMQuad=zeros(ndofpD,ndofpD);
    for i=1:NZ
        if StrongMass
            [massM]=elem_J(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),...
                true,false,ndofpZ,[],QuadDensity(:,i),[],[],[],[],[],[]);
        else
            switch func2str(VBasis)
                case 'Q1Basis'
                    [massM]=elemQ1(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),...
                        true,false,[],QuadDensity(:,i),[],[],[],[],[],[],[],[]);
				otherwise
                    [massM]=elemQ2(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,refMV,QuadGradVBasis,QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),...
                        true,false,ndofpZ,[],QuadDensity(:,i),[],[],[],[],[],[],[],[],[]);
            end
                
        end
        MMQuad(Quadmap(:,i),Quadmap(:,i))=MMQuad(Quadmap(:,i),Quadmap(:,i))+massM;
    end
    massN=sum(MMQuad);
end