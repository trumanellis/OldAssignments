%% Variable Preallocation

% Initialize Quad state variables
acceleration=zeros(2,ndofpD);
OLDvelocity=zeros(2,ndofpD);
NEWvelocity=zeros(2,ndofpD);

% Initialize the Zonal state variables
if HighOrderEnergy
    OLDenergy=zeros(npdof,NZ);
    NEWenergy=zeros(npdof,NZ);
else
    OLDenergy=zeros(1,NZ);
    NEWenergy=zeros(1,NZ);
end
massZ=zeros(1,NZ);

% Preallocate thermodynamic variables
OLDdensity=zeros(npdof,NZ);
NEWdensity=zeros(npdof,NZ);
pressure=zeros(npdof,NZ);

cornerMass=zeros(npdof,NZ);
cornerForce=zeros(2,ndofpZ,NZ);

% Matrix for converting subcess thermal properties to Q1 dofs
localMP=zeros(npdof,npdof,NZ);
localMV=zeros(ndofpZ,ndofpZ,NZ);

% Preallocate Storage for force quadrature points
QP=QuadOrder^2;
QuadPressure=zeros(QP,NZ);
QuadDensity=zeros(QP,NZ);
QuadJacobian=zeros(2,2,QP,NZ);
QuadInvJacobian=zeros(2,2,QP,NZ);
QuadDetJacobian=zeros(QP,NZ);
QuadVBasis=zeros(ndofpZ,QP);
QuadGradVBasis=zeros(ndofpZ,2,QP);
[QuadPts2d QuadWgts2d]= GaussLegendreWeights2d(QuadOrder);
if npdof ~= 2
    refMP=zeros(npdof,npdof,QP);
    QuadPBasis=zeros(npdof,QP);
else
    refMP=zeros(npdof,npdof,QP,2);
    QuadPBasis=zeros(npdof,QP,2);
end
refMV=zeros(ndofpZ,ndofpZ,QP);
for n=1:QP
    QuadVBasis(:,n)=feval(VBasis,QuadPts2d(n,1),QuadPts2d(n,2));
    QuadGradVBasis(:,:,n)=feval(GradVBasis,QuadPts2d(n,1),QuadPts2d(n,2));
    refMV(:,:,n)=QuadVBasis(:,n)*QuadVBasis(:,n)';
    if npdof ~= 2
        QuadPBasis(:,n)=PBasis(QuadPts2d(n,1),QuadPts2d(n,2));
        refMP(:,:,n)=QuadPBasis(:,n)*QuadPBasis(:,n)';
    else
        QuadPBasis(:,n,1)=PBasis(QuadPts2d(n,1),QuadPts2d(n,2),1);
        QuadPBasis(:,n,2)=PBasis(QuadPts2d(n,1),QuadPts2d(n,2),2);
        refMP(:,:,n,1)=QuadPBasis(:,n,1)*QuadPBasis(:,n,1)';
        refMP(:,:,n,2)=QuadPBasis(:,n,2)*QuadPBasis(:,n,2)';
    end
end
for N=1:NZ
    for n=1:QP
        QuadJacobian(:,:,n,N)=Jacobian(allnodes(:,Quadmap(:,N)),QuadPts2d(n,:));
        QuadInvJacobian(:,:,n,N)=inv(QuadJacobian(:,:,n,N));
        QuadDetJacobian(n,N)=det(QuadJacobian(:,:,n,N));
    end
end

len=zeros(1,NZ);
cen=zeros(2,NZ);
psi0Z=zeros(1,NZ);
mu=zeros(QP,NZ);
% Construct mu Map
muMap=reshape(1:QP*NZ',QP,NZ)';
switch QuadOrder
    case 2
        muBasis=@Q1dBasis;
    case 3
        muBasis=@Q2dBasis;
    case 4
        muBasis=@Q3dBasis;
end

% Initialize time history arrays for integrated quantities
timedata=zeros(ceil(maxcycle),1);
internalenergy=zeros(ceil(maxcycle),1);
kineticenergy=zeros(ceil(maxcycle),1);
sumMassN=zeros(ceil(maxcycle),1);
cycle=1;

% Initialize looping conditions
stop=false;
t=tstart;