function [h,pts,rho,pressure]=PlotStrongMassContourf(nodes,rhodetJ,ptmapping,varmapping,npts,ptBasis,varBasis,varargin)

NZ=size(ptmapping,2);
nl=npts^2;
temppts=linspace(0,1,npts);
localpts=zeros(2,nl);
for i=1:npts
    for j=1:npts
        localpts(1,i+(j-1)*npts)=temppts(i);
        localpts(2,i+(j-1)*npts)=temppts(j);
    end
end

pts=zeros(2,nl*NZ);
rho=zeros(1,nl*NZ);
pressure=zeros(1,nl*NZ);

globalpts=zeros(2,nl);
globalrho=zeros(1,nl);
for N=1:NZ
    for i=1:nl
        globalpts(1,i)=LocalToGlobal(nodes(1,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalpts(2,i)=LocalToGlobal(nodes(2,ptmapping(:,N)),localpts(:,i),ptBasis);
        massZ=LocalToGlobal(rhodetJ(varmapping(:,N)),localpts(:,i),varBasis);
        switch func2str(ptBasis)
            case 'Q1Basis'
                globalrho(i)=massZ/det(JacobianQ1Basis(nodes(:,ptmapping(:,N)),localpts(:,i)));
            case 'Q2Basis'
                globalrho(i)=massZ/det(JacobianQ2Basis(nodes(:,ptmapping(:,N)),localpts(:,i)));
        end
    end
    pts(:,(N-1)*nl+(1:nl))=globalpts;
    rho((N-1)*nl+(1:nl))=globalrho;
end
surfmapping=zeros((npts-1)^2*NZ,4);
c=1;
c1=1;
for N=1:NZ
    for n1=1:npts-1
        for n2=1:npts-1
            surfmapping(c1,:)=[c c+1 c+npts+1 c+npts];
            c=c+1;
            c1=c1+1;
        end
        c=c+1;
    end
    c=c+npts;
end
h=quadsurf(surfmapping,pts(1,:),pts(2,:),zeros(1,nl*NZ),rho,varargin{:});