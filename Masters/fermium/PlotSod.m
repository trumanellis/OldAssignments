% PlotSod
clear all
close all
sodExact=load('sod_solution.dat');
% figure(1)
% hold on
% plotvarN=NEWvelocity(1,:);
% [X2pts,Nvals]=PlotFEM1line(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
% plot(sodExact(:,2),sodExact(:,5),'k');
% plot(X2pts(1,:),Nvals);
% 
% figure(2)
% hold on
% presplot=NEWdensity(:)';
% [X1pts,Zvals]=PlotFEM1line(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
% plot(sodExact(:,2),sodExact(:,3),'k');
% plot(X1pts(1,:),Zvals);

% save(['Figures/Sod/',func2str(VBasis),func2str(PBasis),2], 'X1pts', 'Zvals', 'X2pts', 'Nvals')

load('Figures/Sod/Q1BasisQ0Basis');
Q1Q0pts1=X1pts;
Q1Q0Z=Zvals;
Q1Q0pts2=X2pts;
Q1Q0N=Nvals;
load('Figures/Sod/Q1BasisQ0Basis2');
Q1Q0pts1r=X1pts;
Q1Q0Zr=Zvals;
Q1Q0pts2r=X2pts;
Q1Q0Nr=Nvals;
load('Figures/Sod/Q2BasisQ1dBasis');
Q2Q1pts1=X1pts;
Q2Q1Z=Zvals;
Q2Q1pts2=X2pts;
Q2Q1N=Nvals;
load('Figures/Sod/Q2BasisQ1dBasis2');
Q2Q1pts1r=X1pts;
Q2Q1Zr=Zvals;
Q2Q1pts2r=X2pts;
Q2Q1Nr=Nvals;

figure(1)
hold on
plot(sodExact(:,2),sodExact(:,5),'k');
% plot(Q1Q0pts2(1,:),Q1Q0N,Q1Q0pts2r(1,:),Q1Q0Nr,Q2Q1pts2(1,:),Q2Q1N,Q2Q1pts2r(1,:),Q2Q1Nr)
plot(Q1Q0pts2r(1,:),Q1Q0Nr,Q2Q1pts2(1,:),Q2Q1N,Q2Q1pts2r(1,:),Q2Q1Nr)
% legend('Exact','Q_1-Q_0','Q_2-P_1','Q_2-Q_1','Q_2-Q_2')

figure(2)
hold on
plot(sodExact(:,2),sodExact(:,3),'k');
% plot(Q1Q0pts1(1,:),Q1Q0Z,Q1Q0pts1r(1,:),Q1Q0Zr,Q2Q1pts1(1,:),Q2Q1Z,Q2Q1pts1r(1,:),Q2Q1Zr)
plot(Q1Q0pts1r(1,:),Q1Q0Zr,Q2Q1pts1(1,:),Q2Q1Z,Q2Q1pts1r(1,:),Q2Q1Zr)
% legend('Exact','Q_1-Q_0','Q_2-P_1','Q_2-Q_1','Q_2-Q_2')