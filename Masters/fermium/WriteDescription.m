% function WriteDescription(figurefile,userno,dtInit,dtMax,FullMassMatrixSolve,MassUpdate,...
%     massQuadOrder,stiffQuadOrder,Qfrac,hgfrac,qquad,qlin,NZx,NZy,jitter)

fid=fopen([SaveLocation,figurefile,'/summary.txt'],'w');
fprintf(fid,['VBasis: ',func2str(VBasis),'\n']);
fprintf(fid,['PBasis: ',func2str(PBasis),'\n']);
fprintf(fid,'NZx: %4.0f, NZy: %4.0f\n',NZx,NZy);
fprintf(fid,'dtInit: %8.6f, dtMax: %8.6f \n',dtInit,dtMax);
fprintf(fid,'FullMassMatrixSolve: %1.0f, MassUpdate: %1.0f \nQuadOrder: %1.0f, \n',...
    FullMassMatrixSolve,MassUpdate,QuadOrder);
fprintf(fid,'Qfrac: %4.2f, hgfrac: %4.2f, qquad: %6.4f, qlin: %6.4f, jitter: %8.6f \n',Qfrac,hgfrac,qquad,qlin,jitter);