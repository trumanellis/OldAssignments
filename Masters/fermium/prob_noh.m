% Noh Problem
% Author: Truman Ellis
% Email: ellis35@llnl.gov

clear all; close all; clc;
format short g
format compact
scrsz = get(0,'ScreenSize');
h0=sfigure(1,'Position',[1 scrsz(4)/8 scrsz(3) scrsz(4)*7/9]);
set(h0,'renderer','painters')
SetSaveLocation
% setNumberOfComputationalThreads(16);
recover=false;%'FiguresNoh7_2009_07_31_10_09/save0.1680.mat'
nohExact=inline('16*(r < t/3)+(1+t./r).*(r > t/3)','r','t');
nohExactV=inline('1*(r > t/3)','r','t');
if recover
    load([SaveLocation,recover]);
    cycle=cycle+1;
    tic
else
    tic

    SaveFigures=false;
    nplots=10;
    nsaves=1;
    % Set initial time increment
    dtInit=1e-4;
    dtMax=1e-3;

    % Set Kinematic and Thermodynamic Basis
    VBasis=@Q2Basis;
    PBasis=@Q2dBasis;
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
    NZx=20;
    NZy=20;

    % Define the mesh size
    xmin=0;
    xmax=1;

    ymin=xmin;
    ymax=xmax;

    % Define the start and stop time
    tstart = 0;
    tstop = .6;
    tplot=tstart+tstop/nplots:tstop/nplots:tstop;
    tsave=tstart+tstop/nsaves:tstop/nsaves:tstop;
    pcounter=1;
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
    
    vref=7; pref=7;

    % Specify the maximum number of time steps to use
    maxcycle=1e5;

    % These numbers are the coefficients for the artificial viscosity
    qquad=1;
    qlin=1;

    if SaveFigures
        figurefile=['Noh',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
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

    %% Variable Preallocation
    Preallocate

    %% Initialize Properties
    % Set the EOS gamma law constant for an ideal gas
    gamma=5/3;

    rhoInit=1;
    eInit=1e-6;
    pInit=(gamma-1)*rhoInit*eInit;

    origin=intersect(find(allnodes(1,:) == 0),find(allnodes(2,:) == 0));
    
    OLDenergy(:)=eInit;
    OLDdensity(:)=rhoInit;
    QuadDensity(:)=rhoInit;
    pressure(:)=pInit;
    QuadPressure(:)=pInit;
    OLDvelocity(1,:)=-allnodes(1,:)./sqrt(allnodes(1,:).^2+allnodes(2,:).^2);
    OLDvelocity(2,:)=-allnodes(2,:)./sqrt(allnodes(1,:).^2+allnodes(2,:).^2);
    OLDvelocity(:,origin)=0;
    
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

        % Solve Linear System
        acceleration=[(MMX\forceN(1,:)')'; (MMY\forceN(2,:)')'];
        
    else
        
        % Solve for acceleration
        acceleration(1,:)=forceN(1,:)./massN;
        acceleration(2,:)=forceN(2,:)./massN;
        
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

%% Update Mesh
   UpdateMesh
   
%% === Work and EOS Phase ===   
   WorkandEOS
   
%% === VISUALIZATION AND POST-PROCESSING PHASE ===
   
    %% =========== PLOT RESULTS ===========
   if t == tplot(pcounter) || stop
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  dt: %8.6f,  Total Energy: %14.12f,  Compute Time: %9.3f \n',cycle,t,timedata(cycle-1)-timedata(cycle-2),internalenergy(cycle)+kineticenergy(cycle),toc);
       
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
	   end
       [h1,X2pts,Nvals]=PlotFEMContourf2(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
       PlotFEMMesh(allnodes,Quadmap,10,VBasis,'k');
       axis tight
       daspect([1 1 1e-1]);
	   axis([.3*xmin .3 .3*ymin .3 0 1 0 1])
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
       plot(linspace(1e-3,.5,200),nohExactV(linspace(1e-3,.5,200),t),lineoutX2,lineoutN,'r.')
	   set(gca,'Xlim',[0 .5],'Ylim',[0,1.2])
       legend('Exact','FEM','Location','NorthEast')
       title('Velocity Magnitude','FontWeight','demi','FontSize',12)
	   grid on
               
	   if plotpres
		   subplot('position',[.05,.05,.5,.4]);
	   else
		   sfigure(2);
		   subplot('position',[.05,.1,.5,.8]);
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
	   axis([.3*xmin .3 .3*ymin .3 0 20 0 20])
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
       plot(linspace(1e-3,.5,200),nohExact(linspace(1e-3,.5,200),t),lineoutX,lineoutZ,'r.')
	   set(gca,'Xlim',[0 .5])
       legend('Exact','FEM','Location','NorthEast')
       title('Density','FontWeight','demi','FontSize',12)
	   grid on

	   if plotpres
		   subplot_title(['Noh: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
	   end
       
       pcounter=pcounter+1;
       if SaveFigures
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
	Vstr=func2str(VBasis);
	Dstr=func2str(PBasis);
	set(1,'PaperPositionMode','auto'); 
	set(2,'PaperPositionMode','auto'); 
	if hgfrac
		if StrongMass
			saveas(1,['../Thesis/Figures/NohV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
			saveas(2,['../Thesis/Figures/NohD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
		else
			saveas(1,['../Thesis/Figures/NohV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
			saveas(2,['../Thesis/Figures/NohD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
		end
	else
		if StrongMass
			saveas(1,['../Thesis/Figures/NohV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
			saveas(2,['../Thesis/Figures/NohD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
		else
			saveas(1,['../Thesis/Figures/NohV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
			saveas(2,['../Thesis/Figures/NohD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
		end
	end
end
runtime=toc
if SaveFigures
	fprintf(fid,'runtime: %8.6f \n',runtime);
	fclose(fid);
end