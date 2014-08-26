clear all;
close all;
clc;

VBasis=@Q2Basis;

pnodes=[0 1 1.5 .3 .3 1 1 .5 .75;0 0 1.2 1 -.2 .5 1.3 .5 .5];
cnodes=[0 1 1 0 .5 1 .5 0 .5;0 0 1 1 0 .5 1 .5 .5];
vals=[0 0 0 1 0 0 0 0 0];

numpts=9;

figure
[h3,cXpts,cZvals]=PlotFEMContourf(cnodes,vals,(1:9)',(1:9)',numpts,VBasis,VBasis,'Linestyle','none');
view(0,90)
axis equal tight
colorbar

figure
[h3,pXpts,pZvals]=PlotFEMContourf(pnodes,vals,(1:9)',(1:9)',numpts,VBasis,VBasis,'Linestyle','none');
view(0,90)
axis equal tight
colorbar

figure
plot(cXpts(1,:),cXpts(2,:),'ko')
axis equal tight

figure
plot(pXpts(1,:),pXpts(2,:),'ko')
axis equal tight