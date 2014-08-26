function h=PlotFEMVectorField(nodes,projvector,mapping,npts,Basis,varargin)

NZ=size(mapping,2);
nl=npts^2;
temppts=linspace(0,1,npts);
localpts=zeros(2,size(temppts,2)^2);
for i=1:length(temppts)
    for j=1:length(temppts)
        localpts(1,j+(i-1)*length(temppts))=temppts(i);
        localpts(2,j+(i-1)*length(temppts))=temppts(j);
    end
end

pts=zeros(2,nl*NZ);
vals=zeros(2,nl*NZ);

globalpts=zeros(2,nl);
globalvals=zeros(2,nl);
for N=1:NZ
    for i=1:nl
        globalpts(1,i)=LocalToGlobal(nodes(1,mapping(:,N)),localpts(:,i),Basis);
        globalpts(2,i)=LocalToGlobal(nodes(2,mapping(:,N)),localpts(:,i),Basis);
        globalvals(1,i)=LocalToGlobal(projvector(1,mapping(:,N)),localpts(:,i),Basis);
        globalvals(2,i)=LocalToGlobal(projvector(2,mapping(:,N)),localpts(:,i),Basis);
    end
    pts(:,(N-1)*nl+(1:nl))=globalpts;
    vals(:,(N-1)*nl+(1:nl))=globalvals;
end

h=quiver(pts(1,:),pts(2,:),vals(1,:),vals(2,:),varargin{:});