function [h,pts,vals]=PlotFEMContourf2(nodes,projvector,ptmapping,varmapping,npts,ptBasis,varBasis,varargin)

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
vals=zeros(2,nl*NZ);

globalpts=zeros(2,nl);
globalvals=zeros(2,nl);
for N=1:NZ
    for i=1:nl
        globalpts(1,i)=LocalToGlobal(nodes(1,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalpts(2,i)=LocalToGlobal(nodes(2,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalvals(1,i)=LocalToGlobal(projvector(1,varmapping(:,N)),localpts(:,i),varBasis);
        globalvals(2,i)=LocalToGlobal(projvector(2,varmapping(:,N)),localpts(:,i),varBasis);
    end
    pts(:,(N-1)*nl+(1:nl))=globalpts;
    vals(:,(N-1)*nl+(1:nl))=globalvals;
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
vals=sqrt(vals(1,:).^2+vals(2,:).^2);
h=quadsurf(surfmapping,pts(1,:),pts(2,:),zeros(1,nl*NZ),vals,varargin{:});