function [sumerr,h]=PlotL2Norm(nodes,ptmapping,ptBasis,var,varmapping,varBasis,...
    quadorder,Jacobian,exactFn)
% Compute L2 Norm 

[pts2d wts2d] = GaussLegendreWeights2d(quadorder);
NZ=size(ptmapping,2);
NQ=length(wts2d);
err=zeros(1,NZ);

for N=1:NZ
    for n=1:NQ
        xglob=LocalToGlobal(nodes(1,ptmapping(:,N)),pts2d(n,:),ptBasis);
        yglob=LocalToGlobal(nodes(2,ptmapping(:,N)),pts2d(n,:),ptBasis);
        varglob=[LocalToGlobal(var(1,varmapping(:,N)),pts2d(n,:),varBasis);
                 LocalToGlobal(var(2,varmapping(:,N)),pts2d(n,:),varBasis)];
        detJ=det(Jacobian(nodes(:,ptmapping(:,N)),pts2d(n,:)));
        err(N)=err(N)+wts2d(n)*(feval(exactFn,xglob,yglob)-varglob)'*(feval(exactFn,xglob,yglob)-varglob)*detJ;
    end
end
err=sqrt(err);
h=PlotFEMContourf(nodes,err,ptmapping,1:NZ,2,ptBasis,@Q0Basis,'Linestyle','none');
sumerr=sqrt(sum(err(:).^2));