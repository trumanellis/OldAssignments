% UpdateMesh
% Integrate accelerations to get velocities
NEWvelocity=OLDvelocity+dt*acceleration;

% ----- Apply Source term to velocities -----
if exist('drivingPt')
    for i=1:length(drivingPt)
        NEWvelocity(1,drivingPt(i))=drivingFnx(allnodes(1,drivingPt(i)),allnodes(2,drivingPt(i)),t+dt);
        NEWvelocity(2,drivingPt(i))=drivingFny(allnodes(1,drivingPt(i)),allnodes(2,drivingPt(i)),t+dt);
    end
end

allnodes=allnodes+dt*NEWvelocity;
for N=1:NZ
    for n=1:QP
        QuadJacobian(:,:,n,N)=Jacobian(allnodes(:,Quadmap(:,N)),QuadPts2d(n,:));
        QuadInvJacobian(:,:,n,N)=inv(QuadJacobian(:,:,n,N));
        QuadDetJacobian(n,N)=det(QuadJacobian(:,:,n,N));
        if QuadDetJacobian(n,N) < 0
            fprintf('Warning Negative Jacobian at cell %4.0f\n',N);
            stop=true;
        end
    end
end