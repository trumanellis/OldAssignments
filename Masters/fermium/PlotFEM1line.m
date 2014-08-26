function [pts,vals]=PlotFEM1line(nodes,projvector,ptmapping,varmapping,npts,ptBasis,varBasis,varargin)

NZ=size(ptmapping,2);
nl=npts;
localpts=[linspace(0,1,npts);.5*ones(1,npts)];

pts=zeros(2,nl*NZ);
vals=zeros(1,nl*NZ);

globalpts=zeros(2,nl);
globalvals=zeros(1,nl);
for N=1:NZ
    for i=1:nl
        globalpts(1,i)=LocalToGlobal(nodes(1,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalpts(2,i)=LocalToGlobal(nodes(2,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalvals(i)=LocalToGlobal(projvector(varmapping(:,N)),localpts(:,i),varBasis);
    end
    pts(:,(N-1)*nl+(1:nl))=globalpts;
    vals((N-1)*nl+(1:nl))=globalvals;
end
% surfmapping=zeros((npts-1)^2*NZ,4);
% c=1;
% c1=1;
% for N=1:NZ
%     for n1=1:npts-1
%         for n2=1:npts-1
%             surfmapping(c1,:)=[c c+1 c+npts+1 c+npts];
%             c=c+1;
%             c1=c1+1;
%         end
%         c=c+1;
%     end
%     c=c+npts;
% end
% h=quadsurf(surfmapping,pts(1,:),pts(2,:),zeros(1,nl*NZ),vals,varargin{:});