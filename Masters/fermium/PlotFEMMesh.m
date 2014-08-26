function PlotFEMMesh(nodes,ptmapping,npts,ptBasis,varargin)

NZ=size(ptmapping,2);
nl=4*npts;
localpts=[linspace(0,1,npts), ones(1,npts), linspace(1,0,npts), zeros(1,npts);
          zeros(1,npts), linspace(0,1,npts), ones(1,npts),linspace(1,0,npts)];

globalpts=zeros(2,npts*4);
hold on;
for N=1:NZ
    for i=1:nl
        globalpts(1,i)=LocalToGlobal(nodes(1,ptmapping(:,N)),localpts(:,i),ptBasis);
        globalpts(2,i)=LocalToGlobal(nodes(2,ptmapping(:,N)),localpts(:,i),ptBasis);
    end
    plot(globalpts(1,:),globalpts(2,:),varargin{:});
end
hold off