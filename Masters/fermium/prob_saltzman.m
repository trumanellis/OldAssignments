% Saltzman Piston Problem
% Author: Truman Ellis
% Email: ellis35@llnl.gov

clear all; close all; clc;
remote = 1;
format short g
format compact
scrsz = get(0,'ScreenSize');
h0=sfigure(1,'Position',[1 scrsz(4)/8 scrsz(3) scrsz(4)*7/9]);
set(h0,'renderer','painters')
SetSaveLocation
recover=false;%'FiguresSaltzman7_2009_09_08_12_12/save0.9758.mat'
if recover
    load([SaveLocation,recover]);
    cycle=cycle+1;
    tic
else
    tic

    SaveFigures=false;
    nplots=80;
    nsaves=80;
    % Set initial time increment
    dtInit=1e-4;
    dtMaxFinal=1e-6;

    % Set Kinematic and Thermodynamic Basis
    VBasis=@Q1Basis;
    PBasis=@Q0Basis;
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
    NZx=50;
    NZy=10;

    % Define the mesh size
    xmin=0;
    xmax=1;

    ymin=0;
    ymax=.1*xmax;

    % Define the start and stop time
    tstart = 0;
    tstop = .925;
	tplot=tstop*cos(linspace(pi/2-.1,0,nplots));
	tsave=tstop*cos(linspace(pi/2-.1,0,nsaves));
%      tplot=sort(tstop-logspace(log10(.1),log10(tstop),nplots)+.1);
%      tsave=sort(tstop-logspace(log10(.1),log10(tstop),nsaves)+.1);
%     tplot=tstart+tstop/nplots:tstop/nplots:tstop;
%     tsave=tstart+tstop/nsaves:tstop/nsaves:tstop;
    pcounter=1;
    if SaveFigures    
        scounter=1;
    else
        scounter=nsaves;
    end

    %Define a random "jittering" factor -- this determines the magnitude of
    %mesh distortion
    jitter=0.0;
    % Define Saltzman mesh perturbation
    SaltzmanPerturb=1;
    
    vref=5; pref=5;

    % Specify the maximum number of time steps to use
    maxcycle=1e7;
    dtMax=dtMaxFinal;

    % These numbers are the coefficients for the artificial viscosity
    qquad=2;
    qlin=.25;

    if SaveFigures
        figurefile=['Saltzman',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
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

    OLDenergy(:)=eInit;
    OLDdensity(:)=rhoInit;
    QuadDensity(:)=rhoInit;
    pressure(:)=pInit;
    QuadPressure(:)=pInit;
    
    % Set Driver
    drivingPt=[ones(1,length(find(allnodes(1,:) == 0)));find(allnodes(1,:) == 0)];
    drivingFnx=@(a,b,t) 1;
    drivingFny=@(a,b,t) 0;
    for i=1:size(drivingPt,2)
        OLDvelocity(drivingPt(1,i),drivingPt(2,i))=drivingFnx(allnodes(1,drivingPt(2,i)),allnodes(2,drivingPt(2,i)),tstart);
    end
    
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
    dtMax=(dtMaxFinal-dtInit)*t+dtInit;
    
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
        % Impose zero boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(2,index)=0;
        end
        % Impose zeros boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
        end

        % Solve Linear System
        acceleration=[(MMX\forceN(1,:)')'; (MMY\forceN(2,:)')'];
        
    elseif MassUpdate
        
        % Solve for acceleration
		acceleration(1,:)=forceN(1,:)./massN;
		acceleration(2,:)=forceN(2,:)./massN;
		
		% Impose Boundary Conditions
		% Impose symmetric boundary at bottom wall
		for i=1:sum(nbnodes(1))
			index=bdof(i);
			acceleration(2,index)=0;
		end
		% Impose zero boundary at right wall
		for i=1+nbnodes(1):sum(nbnodes(1:2))
			index=bdof(i);
			acceleration(:,index)=[0;0];
		end
		% Impose symmetric boundary at top wall
		for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
			index=bdof(i);
			acceleration(2,index)=0;
		end
		% Impose zero boundary at left wall
		for i=1+sum(nbnodes(1:3)):sum(nbnodes)
			index=bdof(i);
			acceleration(:,index)=[0;0];
		end
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
		% Impose zero boundary at right wall
		for i=1+nbnodes(1):sum(nbnodes(1:2))
			index=bdof(i);
			acceleration(:,index)=[0;0];
		end
		% Impose symmetric boundary at top wall
		for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
			index=bdof(i);
			acceleration(2,index)=0;
		end
		% Impose zero boundary at left wall
		for i=1+sum(nbnodes(1:3)):sum(nbnodes)
			index=bdof(i);
			acceleration(:,index)=[0;0];
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
		
		if ~remote
			if plotpres
				h=sfigure(1);
			else
				h1=sfigure(1,'Position',[0 0 scrsz(3) 7/8*scrsz(4)]);
				h2=sfigure(2,'Position',[0 0 scrsz(3) 7/8*scrsz(4)]);
			end
			
			plotvarN=NEWvelocity(1,:);
			
			if plotpres
				subplot('position',[.05,.55,.5,.4]);
			else
				sfigure(1);
				subplot('position',[.05,.1,.5,.8]);
			end
			[h1,X2pts,Nvals]=PlotFEMContourf(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
			PlotFEMMesh(allnodes,Quadmap,10,VBasis,'w');
			axis tight
			title('x-Velocity','FontWeight','demi','FontSize',14)
			az=0;
			el=90;
			view(az,el);
			colorbar
			
			lineoutX2=X2pts(1,:);
			lineoutN=Nvals(1,:);
			
			if plotpres
				subplot('position',[.6,.55,.35,.4]);
			else
				sfigure(1);
				subplot('position',[.6,.1,.35,.8]);
			end
			plot(lineoutX2,lineoutN,'r.')
			set(gca,'Xlim',[t 1],'Ylim',[0,1.2])
			title('x-Velocity','FontWeight','demi','FontSize',14)
			
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
			title('Density','FontWeight','demi','FontSize',14)
			view(az,el);
			colorbar
			
			lineoutX=sqrt(X1pts(1,:).^2+X1pts(2,:).^2);
			lineoutZ=Zvals;
			
			if plotpres
				subplot('position',[.6,.05,.35,.4]);
			else
				sfigure(2);
				subplot('position',[.6,.1,.35,.8]);
			end
			plot(lineoutX,lineoutZ,'r.')
			set(gca,'Xlim',[t 1])
			title('Density','FontWeight','demi','FontSize',14)
			
			if plotpres
				subplot_title(['Noh: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
			end
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
	if ~remote
		set(1,'PaperPositionMode','auto');
		set(2,'PaperPositionMode','auto');
		if hgfrac
			if StrongMass
				saveas(1,['../Thesis/Figures/SaltzmanV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
				saveas(2,['../Thesis/Figures/SaltzmanD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM'],'png');
			else
				saveas(1,['../Thesis/Figures/SaltzmanV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
				saveas(2,['../Thesis/Figures/SaltzmanD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON'],'png');
			end
		else
			if StrongMass
				saveas(1,['../Thesis/Figures/SaltzmanV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
				saveas(2,['../Thesis/Figures/SaltzmanD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM'],'png');
			else
				saveas(1,['../Thesis/Figures/SaltzmanV_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
				saveas(2,['../Thesis/Figures/SaltzmanD_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF'],'png');
			end
		end
	else
		if hgfrac
			if StrongMass
				save(['./Figures/Saltzman_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON_SM.mat']);
			else
				save(['./Figures/Saltzman_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgON.mat']);
			end
		else
			if StrongMass
				save(['./Figures/Saltzman_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF_SM.mat']);
			else
				save(['./Figures/Saltzman_',Vstr(1:2),Dstr(1:2),'_',num2str(NZx),'x',num2str(NZy),'_hgOFF.mat']);
			end
		end
	end
end
toc
