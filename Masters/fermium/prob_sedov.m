% Sedov Problem
% Author: Truman Ellis
% Email: ellis35@llnl.gov

clear all; close all; clc;
format short g
format compact
scrsz = get(0,'ScreenSize');
% h0=sfigure(1,'Position',[0 0 scrsz(3) 7/8*scrsz(4)]);
% h0=sfigure(1,'Position',[1 100 600 500]);
% h0=sfigure(1,'Position',[1 scrsz(4)/8 scrsz(3) scrsz(4)*7/9]);
% set(h0,'renderer','painters')
SetSaveLocation
recover=false;%'FiguresAcousticWave5_2009_07_10_15_02/save0.5000.mat'
sedovExact=load('sedov1.dat');
if recover
    load(['/p/lscratchd/ellis35/FigureFiles/',recover]);
    cycle=cycle+1;
    tic
else
    tic

    SaveFigures=false;
    nplots=10;
    nsaves=1;
    % Set initial time increment
    dtInit=1e-4;
    dtMax=5e-2;

    % Set Kinematic and Thermodynamic Basis
    VBasis=@Q1Basis;
    PBasis=@Q1dBasis;
    HighOrderEnergy=false;
    StrongMass=false;
    
    FullMassMatrixSolve=false;
    MassUpdate=false;
    % Set Quadrature Order
    QuadOrder=3;

    % Use this to turn artificial viscosity on/off
    Qfrac=1;

    % Use this to turn anti-hourglass forces on/off (Only works with Q1Q0)
    hgfrac=0;

    % Define the number of zones
    NZx=10;
    NZy=10;

    % Define the mesh size
    isFullMesh=false;
    xmin=0;
    xmax=1.2;

    ymin=xmin;
    ymax=xmax;

    % Define the start and stop time
    tstart = 0;
    tstop = 1;
% 	tplot=[dtInit,linspace(tstart,tstop,nplots+1)];
    tplot=tstart+tstop/nplots:tstop/nplots:tstop;
    tsave=tstart+tstop/nsaves:tstop/nsaves:tstop;
    pcounter=1;
% 	nplots = nplots+1;
    if SaveFigures    
        scounter=1;
    else
        scounter=nsaves;
    end

    %Define a random "jittering" factor -- this determines the magnitude of
    %mesh distortion
    jitter=0.0;
    % Define Mesh rotations
    RotateMesh=0;
    
    % Define axis and value for lineout before transformation
    lineoutAxis=[2;0];
    vref=5; pref=5;

    % Specify the maximum number of time steps to use
    maxcycle=1e5;

    % These numbers are the coefficients for the artificial viscosity
    qquad=1;
    qlin=1;

    if SaveFigures
        figurefile=['Sedov',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
        if ~exist([SaveLocation,figurefile])
            unix(['mkdir ',SaveLocation,figurefile]);
        else
            figurefile=[figurefile,datestr(clock, '_SS')];
            unix(['mkdir ',SaveLocation,figurefile]);
        end
        WriteDescription
    end

    %% Mesh Assembly
    MeshAssembly
    
%     thetas=linspace(0,pi/2,5);
%     r=xmax/NZx;
%     allnodes(1,Quadmap([2 6 3 7 4],1))=r*cos(thetas);
%     allnodes(2,Quadmap([2 6 3 7 4],1))=r*sin(thetas);
%     allnodes(1,Quadmap([5 9 8],1))=.5*r*cos([0 pi/4 pi/2]);
%     allnodes(2,Quadmap([5 9 8],1))=.5*r*sin([0 pi/4 pi/2]);

    %% Variable Preallocation
    Preallocate

    %% Initialize Properties
    % Set the EOS gamma law constant for an ideal gas
    gamma=1.4;

    rhoInit=1;
    eInit=1e-6;
    if func2str(VBasis) == 'Q1Basis'
        eBigInit=(xmax-xmin)*(ymax-ymin)/(4*xmax*ymax)/det(JacobianQ1Basis(allnodes(:,topo(:,1)),[.5 .5]));
    else
        eBigInit=(xmax-xmin)*(ymax-ymin)/(4*xmax*ymax)/det(JacobianQ2Basis(allnodes(:,Quadmap(:,1)),[.5 .5]));
    end
    pInit=(gamma-1)*rhoInit*eInit;
    pBigInit=(gamma-1)*rhoInit*eBigInit;

    if isFullMesh
        if mod(NZx,2)
            if strcmp(func2str(VBasis),'Q2Basis')
                origin=intersect(find(allnodes(1,Quadmap(9,:)) == 0),find(allnodes(2,Quadmap(9,:)) == 0));
            else
                origin=intersect(find((allnodes(1,Quadmap(1,:))+allnodes(1,Quadmap(3,:))) == 0),find((allnodes(2,Quadmap(1,:))+allnodes(2,Quadmap(3,:))) == 0));
            end
        else
            origin=[intersect(find(allnodes(1,Quadmap(1,:)) == 0),find(allnodes(2,Quadmap(1,:)) == 0)),
                    intersect(find(allnodes(1,Quadmap(2,:)) == 0),find(allnodes(2,Quadmap(2,:)) == 0)),
                    intersect(find(allnodes(1,Quadmap(3,:)) == 0),find(allnodes(2,Quadmap(3,:)) == 0)),
                    intersect(find(allnodes(1,Quadmap(4,:)) == 0),find(allnodes(2,Quadmap(4,:)) == 0))];
            eBigInit=.25*eBigInit;
            pBigInit=(gamma-1)*rhoInit*eBigInit;
        end
    else
        origin=intersect(find(allnodes(1,Quadmap(1,:)) == 0),find(allnodes(2,Quadmap(1,:)) == 0));
    end
    
    OLDenergy(:)=eInit;
    OLDenergy(:,origin)=eBigInit;
    OLDdensity(:)=rhoInit;
    QuadDensity(:)=rhoInit;
    pressure(:)=pInit;
    pressure(1:npdof,origin)=pBigInit;
    QuadPressure(:)=pInit;
    QuadPressure(:,origin)=pBigInit;
    
    for N=1:NZ
        if ~StrongMass
            for n=1:QP
                localMP(:,:,N)=localMP(:,:,N)+QuadWgts2d(n)*refMP(:,:,n)*QuadDetJacobian(n,N);
            end
            cornerMass(:,N)=localMP(:,:,N)*OLDdensity(:,N);
            massZ(N)=sum(cornerMass(:,N),1);
        else
            for n=1:npdof
                detJ=det(Jacobian(allnodes(:,Quadmap(:,N)),pdofpts(:,n)));
                OLDdensity(n,N)=OLDdensity(n,N)*detJ;
                pressure(n,N)=pressure(n,N)*detJ;
            end
            NEWdensity=OLDdensity;
            QuadDensity(:,N)=(OLDdensity(:,N)'*QuadPBasis)';
            QuadPressure(:,N)=(pressure(:,N)'*QuadPBasis)';
            for n=1:QP
                localMP(:,:,N)=localMP(:,:,N)+QuadWgts2d(n)*refMP(:,:,n);
            end
            cornerMass(:,N)=localMP(:,:,N)*OLDdensity(:,N);
            massZ(N)=sum(cornerMass(:,N),1);
        end
    end
    
    % Assemble Nodal Masses if Mass is not updated every cycle
    InitMassMatrix
end
    
%% === Generation Visualization Phase ===

% Step 2) Now we loop over time steps and solve the hydro equations at each
% Lagrangian time step

% === Lagrangian Phase ===

%% === Time Stepping Loop ===
for cycle=cycle:maxcycle
    if stop==true
        break;
    end
    
    % Compute stable time increment
    StableTimeStep
    
%% === Acceleration Phase ===
    
    % Assemble Mass Matrix and Force Vector
    UpdateMassAndForce

    if FullMassMatrixSolve
        %% Mass Matrix Assembly
        
        if ~isFullMesh
            % Impose Boundary Conditions
            % Impose symmetric boundary at bottom wall
            for i=1:sum(nbnodes(1))
                index=bdof(i);
                MMY(index,:)=0;
                MMY(:,index)=0;
                MMY(index,index)=1;
                forceN(2,index)=0;
            end
            % Impose symmetric boundary at left wall
            for i=1+sum(nbnodes(1:3)):sum(nbnodes)
                index=bdof(i);
                MMX(index,:)=0;
                MMX(:,index)=0;
                MMX(index,index)=1;
                forceN(1,index)=0;
            end
        end

        % Solve Linear System
        acceleration=[(MMX\forceN(1,:)')'; (MMY\forceN(2,:)')'];
        
    else
        
        % Solve for acceleration
        acceleration(1,:)=forceN(1,:)./massN;
        acceleration(2,:)=forceN(2,:)./massN;
        
        if ~isFullMesh
            % Impose Boundary Conditions
            % Impose symmetric boundary at bottom wall
            for i=1:sum(nbnodes(1))
                index=bdof(i);
                acceleration(2,index)=0;
            end
            % Impose symmetric boundary at left wall
            for i=1+sum(nbnodes(1:3)):sum(nbnodes)
                index=bdof(i);
                acceleration(1,index)=0;
            end
        end
    end

%% Update Mesh
   UpdateMesh
   
%% === Work and EOS Phase ===   
   WorkandEOS
   
%% === VISUALIZATION AND POST-PROCESSING PHASE ===
   
    %% =========== PLOT RESULTS ===========
   if cycle > 0 && t == tplot(pcounter) || stop
	   if cycle > 2
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  dt: %8.6f,  Total Energy: %14.12f,  Compute Time: %9.3f \n',cycle,t,timedata(cycle-1)-timedata(cycle-2),internalenergy(cycle)+kineticenergy(cycle),toc);
	   end
	   % paper = 0, presentation = 1
	   plotpres=1;
	   
	   if plotpres
		   h=sfigure(1,'Position',[0 0 .7*scrsz(3) 7/8*scrsz(4)]);
	   else
		   h1=sfigure(1,'Position',[0 0 scrsz(3) 7/8*scrsz(4)]);
		   h2=sfigure(2,'Position',[0 0 scrsz(3) 7/8*scrsz(4)]);
	   end

       plotvarN=NEWvelocity;
               
	   if plotpres
		   subplot('position',[.05,.55,.5,.4]);
	   else
		   sfigure(1);
		   subplot('position',[.05,.1,.5,.8]);
% 		   subplot('position',[.025,.1,.45,.8]);
	   end
       [h1,X2pts,Nvals]=PlotFEMContourf2(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
       PlotFEMMesh(allnodes,Quadmap,10,VBasis,'k');
       axis tight
       daspect([1 1 1e-1]);
       title('Velocity Magnitude','FontWeight','demi','FontSize',12)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
       lineoutX2=sqrt(X2pts(1,:).^2+X2pts(2,:).^2);
       lineoutN=Nvals;
       
	   if plotpres
		   subplot('position',[.55,.55,.4,.4]);
	   else
		   sfigure(1);
		   subplot('position',[.6,.1,.35,.8]);
	   end
       plot(sedovExact(:,1),sedovExact(:,5),lineoutX2,lineoutN,'r.')
       legend('Exact','FEM','Location','NorthEast')
       title('Velocity Magnitude','FontWeight','demi','FontSize',12)
	   grid on
               
	   if plotpres
		   subplot('position',[.05,.05,.5,.4]);
	   else
		   sfigure(2);
		   subplot('position',[.05,.1,.5,.8]);
% 		   subplot('position',[.525,.1,.45,.8]);
	   end
       presplot=NEWdensity(:)';
       if ~StrongMass
           [h3,X1pts,Zvals]=PlotFEMContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
       else
           [h3,X1pts,Zvals]=PlotStrongMassContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
       end
       PlotFEMMesh(allnodes,Quadmap,10,VBasis,'w');
       axis tight
       daspect([1 1 1e-1]);
       title('Density','FontWeight','demi','FontSize',12)
       view(az,el);
       colorbar
               
       lineoutX=sqrt(X1pts(1,:).^2+X1pts(2,:).^2);
       lineoutZ=Zvals;
       
	   if plotpres
		   subplot('position',[.55,.05,.4,.4]);
	   else
		   sfigure(2);
		   subplot('position',[.6,.1,.35,.8]);
	   end
       plot(sedovExact(:,1),sedovExact(:,3),lineoutX,lineoutZ,'r.')
       legend('Exact','FEM','Location','NorthEast')
       title('Density','FontWeight','demi','FontSize',12)
	   grid on

% 	   if plotpres
% 		   subplot_title(['Sedov: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
% 	   end
       
       pcounter=pcounter+1;
       if SaveFigures
		   set(1,'PaperPositionMode','auto');
           % Save each plot frame for conversion to animated gif
           saveas(1,[SaveLocation,figurefile,'/plot',num2str(t,'%06.4f'),'.png']);
       else
           drawnow
       end
   end
   if SaveFigures && (t == tsave(scounter) || stop)
       save([SaveLocation,figurefile,'/save',num2str(t,'%06.4f'),'.mat']);
       scounter=scounter+1;
   end
end
if ~plotpres
	if isFullMesh
		FullMeshPlot='/FullMesh';
	else
		FullMeshPlot='';
	end
	Vstr=func2str(VBasis);
	Dstr=func2str(PBasis);
	set(1,'PaperPositionMode','auto'); 
	set(2,'PaperPositionMode','auto'); 
	if hgfrac
		if StrongMass
			saveas(1,['../Thesis/Figures',FullMeshPlot,'/SedovV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
			saveas(2,['../Thesis/Figures',FullMeshPlot,'/SedovD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
		else
			saveas(1,['../Thesis/Figures',FullMeshPlot,'/SedovV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
			saveas(2,['../Thesis/Figures',FullMeshPlot,'/SedovD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
		end
	else
		if StrongMass
			saveas(1,['../Thesis/Figures',FullMeshPlot,'/SedovV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
			saveas(2,['../Thesis/Figures',FullMeshPlot,'/SedovD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
		else
			saveas(1,['../Thesis/Figures',FullMeshPlot,'/SedovV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
			saveas(2,['../Thesis/Figures',FullMeshPlot,'/SedovD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
		end
	end
end
runtime=toc
if SaveFigures
	fprintf(fid,'runtime: %8.6f \n',runtime);
	fclose(fid);
end