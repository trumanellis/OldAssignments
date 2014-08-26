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
% 
% save(['Figures/Sod/',func2str(VBasis),func2str(PBasis)], 'X1pts', 'Zvals', 'X2pts', 'Nvals')

load('Figures/Sod/Q1BasisQ0Basis');
Q1Q0pts1=X1pts;
Q1Q0Z=Zvals;
Q1Q0pts2=X2pts;
Q1Q0N=Nvals;
% load('Figures/Sod/Q2BasisP1qBasis');
% Q2P1pts1=X1pts;
% Q2P1Z=Zvals;
% Q2P1pts2=X2pts;
% Q2P1N=Nvals;
load('Figures/Sod/Q2BasisQ1dBasis');
Q2Q1pts1=X1pts;
Q2Q1Z=Zvals;
Q2Q1pts2=X2pts;
Q2Q1N=Nvals;
% load('Figures/Sod/Q2BasisQ2dBasis');
% Q2Q2pts1=X1pts;
% Q2Q2Z=Zvals;
% Q2Q2pts2=X2pts;
% Q2Q2N=Nvals;

figure(1)
hold on
plot(sodExact(:,2),sodExact(:,5),'k');
% plot(Q1Q0pts2(1,:),Q1Q0N,Q2P1pts2(1,:),Q2P1N,Q2Q1pts2(1,:),Q2Q1N,Q2Q2pts2(1,:),Q2Q2N)
plot(Q1Q0pts2(1,:),Q1Q0N,Q2Q1pts2(1,:),Q2Q1N)
% legend('Exact','Q_1-Q_0','Q_2-P_1','Q_2-Q_1','Q_2-Q_2')

figure(2)
hold on
plot(sodExact(:,2),sodExact(:,3),'k');
plot(Q1Q0pts1(1,:),Q1Q0Z,Q2Q1pts1(1,:),Q2Q1Z)
% plot(Q1Q0pts1(1,:),Q1Q0Z,Q2P1pts1(1,:),Q2P1Z,Q2Q1pts1(1,:),Q2Q1Z,Q2Q2pts1(1,:),Q2Q2Z)
% legend('Exact','Q_1-Q_0','Q_2-P_1','Q_2-Q_1','Q_2-Q_2')