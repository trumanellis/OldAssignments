function err=L2Norm1(nodes,ptmapping,ptBasis,var,varmapping,varBasis,...
    quadorder,Jacobian,exactFn)
% Compute L2 Norm 

[pts2d wts2d] = GaussLegendreWeights2d(quadorder);
NZ=size(ptmapping,2);
NQ=length(wts2d);
err=0;

for N=1:NZ
    for n=1:NQ
        xglob=LocalToGlobal(nodes(1,ptmapping(:,N)),pts2d(n,:),ptBasis);
        yglob=LocalToGlobal(nodes(2,ptmapping(:,N)),pts2d(n,:),ptBasis);
        varglob=LocalToGlobal(var(varmapping(:,N)),pts2d(n,:),varBasis,N);
        detJ=det(Jacobian(nodes(:,ptmapping(:,N)),pts2d(n,:)));
        err=err+wts2d(n)*(feval(exactFn,xglob,yglob)-varglob)'*(feval(exactFn,xglob,yglob)-varglob)*detJ;
    end
end
err=sqrt(err);