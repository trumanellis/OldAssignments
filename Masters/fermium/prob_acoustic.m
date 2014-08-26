% Acoustic Wave Problem
% Author: Truman Ellis
% Email: ellis35@llnl.gov

clear all; close all; clc;
format long g
format compact
scrsz = get(0,'ScreenSize');
h0=sfigure(1,'Position',[1 scrsz(4)/8 scrsz(3) scrsz(4)*7/9]);
set(h0,'renderer','painters')
SetSaveLocation
% setNumberOfComputationalThreads(16);
recover=false;%'FiguresAcousticWave7_2009_08_04_08_36/save0.8000.mat'
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
    dtInit=1e-3;
    dtMax=1e-2;

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
    Qfrac=0;

    % Use this to turn anti-hourglass forces on/off (Only works with Q1 elements)
    hgfrac=1;

    % Define the number of zones
    NZx=32;
    NZy=32;

    % Define the mesh size
    xmin=-1.2;
    xmax=1.2;

    ymin=xmin;
    ymax=xmax;

    % Define the start and stop time
    tstart = 0;
    tstop = 1;
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
    
    % Define axis and value for lineout before transformation
    lineoutAxis=[2;0];
    vref=5; pref=5;

    % Specify the maximum number of time steps to use
    maxcycle=1e5;

    % These numbers are the coefficients for the artificial viscosity
    qquad=1;
    qlin=1;

    if SaveFigures
        figurefile=['Acoustic',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
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
    vsound=1;

    rhoInit=1;
    pInit=vsound^2*rhoInit/gamma;
    eInit=pInit/((gamma-1)*rhoInit);

    OLDenergy(:)=eInit;
    OLDdensity(:)=rhoInit;
    QuadDensity(:)=rhoInit;
    pressure(:)=pInit;
    QuadPressure(:)=pInit;
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
    
    % Set Driver
    drivingPt=[intersect(find(allnodes(1,:) == 0),find(allnodes(2,:) == 0))];
    amp=sqrt(1e-6/massN(drivingPt));
%     amp=.02*60/NZx;
% 	amp=.01;
    omega=2*2*pi;
    if exist('RotateMesh') && RotateMesh > 0
        rot=RotMat*[1;0];
    else
        rot=[1;0];
    end
    drivingFnx=@(a,b,t) amp*rot(1)*sin(omega*t);
    drivingFny=@(a,b,t) amp*rot(2)*sin(omega*t);
    OLDvelocity(1,drivingPt(1))=drivingFnx([],[],tstart);
    OLDvelocity(2,drivingPt(1))=drivingFny([],[],tstart);
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
        % Impose zero boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
        end
        % Impose symmetric boundary at left wall
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
        % Impose zero boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            acceleration(:,index)=[0;0];
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
       
       h=sfigure(1);
       
       if exist('RotateMesh') && abs(RotateMesh) > 0
           plotvarN=invRotMat*NEWvelocity;
           plotvarX=invRotMat*allnodes;
       else
           plotvarN=NEWvelocity;
           plotvarX=allnodes;
       end

%        subplot(2,2,1)
	   subplot(1,2,1)
       [h1,X2pts,Nvals]=PlotFEMContourf(allnodes,plotvarN(1,:),Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
%        PlotFEMMesh(allnodes,Quadmap,4,VBasis,'k');
       axis tight
       axis([xmin xmax ymin ymax -amp amp -1e-3 1e-3])
       daspect([1 1 1e-2]);
       title('x-Velocity','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
	   
	   subplot(1,2,2)
	   [h1,X2pts,Nvals]=PlotFEMContourf(allnodes,plotvarN(2,:),Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
%        PlotFEMMesh(allnodes,Quadmap,4,VBasis,'k');
       axis tight
       axis([xmin xmax ymin ymax -amp amp -1e-3 1e-3])
       daspect([1 1 1e-2]);
       title('y-Velocity','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
% %        subplot(2,2,3)
% 	   subplot(1,2,2)
%        presplot=pressure(:)';
% %        presplot=mean(pressure,1);
% %        PressureMap=1:NZ';
%        if ~StrongMass
%            [h3,X1pts,Zvals]=PlotFEMContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
% %            [h3,X1pts,Zvals]=PlotFEMCenterContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
%        else
%            [h3,X1pts,Zvals]=PlotStrongMassContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
%        end
% %        PlotFEMMesh(allnodes,Quadmap,4,VBasis,'k');
%        axis([xmin xmax ymin ymax .5985 .6015 .5985 .6015])
%        axis tight
%        daspect([1 1 1e-2]);
%        title('Pressure','FontWeight','demi','FontSize',14)
%        view(az,el);
%        colorbar

%        subplot(2,2,4)
%        enrgplot=NEWenergy(:)';
%        if HighOrderEnergy
%            EMap=PressureMap';
%            EBasis=PBasis;
%        else
%            EMap=1:NZ;
%            EBasis=@Q0Basis;
%        end
%        [h4]=PlotFEMContourf(allnodes,enrgplot,Quadmap,EMap,pref,VBasis,EBasis,'Linestyle','none');
%        axis tight
%        daspect([1 1 1e-2]);
%        title('Energy','FontWeight','demi','FontSize',14)
%        view(az,el);
%        colorbar
       
       if exist('RotateMesh') && abs(RotateMesh) > 0
           unrotX1=invRotMat*X1pts;
           unrotX2=invRotMat*X2pts;
           lineoutX=unrotX1(1,loz);
           lineoutZ=Zvals(loz);
           lineoutX2=unrotX2(1,lon);
           lineoutN=Nvals(lon);
       else
           lineoutX=X1pts(1,loz);
           lineoutX2=X2pts(1,lon);
           lineoutZ=Zvals(loz);
           lineoutN=Nvals(lon);
       end
%        lineoutX=lineoutX(1:pref*ZSTRIDE);
%        lineoutZ=lineoutZ(1:pref*ZSTRIDE);
       lineoutX2=lineoutX2(1:vref*ZSTRIDE);
       lineoutN=lineoutN(1:vref*ZSTRIDE);
       
%        subplot(2,2,2)
%        plot(lineoutX2,lineoutN,'g-')
%        title('x-Velocity','FontWeight','demi','FontSize',14)
%        axis([xmin xmax -.02 .02])
%        
%        subplot(2,2,4)
%        plot(lineoutX,lineoutZ,'b-')
%        title('Pressure','FontWeight','demi','FontSize',14)
%        axis([xmin xmax .58 .62])
       
%        subplot(2,2,[2 4])
% %        [AX,H1,H2]=plotyy(timedata(1:cycle),internalenergy(1:cycle)-internalenergy(1),timedata(1:cycle),kineticenergy(1:cycle));
% %        set(get(AX(1),'Ylabel'),'String','Additional Internal Energy','FontWeight','demi','FontSize',10)
% %        set(get(AX(2),'Ylabel'),'String','Kinetic Energy','FontWeight','demi','FontSize',10)       
% 
% %        plot(timedata(1:cycle),internalenergy(1:cycle)+kineticenergy(1:cycle),'k',...
% %            timedata(1:cycle),internalenergy(1:cycle),'r',...
% %            timedata(1:cycle),kineticenergy(1:cycle),'b')
% %        legend('Total Energy','Internal Energy','Kinetic Energy','Location','NorthWest')
% %        axis([0,timedata(cycle),0 1.4*(internalenergy(cycle)+kineticenergy(cycle))])
% 
%        [AX,H1,H2]=plotyy(lineoutX,lineoutZ,lineoutX2,lineoutN);
%        set(AX(1),'YLim',[.59 .61])
%        set(AX(2),'YLim',[-.015 .015])
%        set(AX(1),'YTick',linspace(.59,.61,5))
%        set(AX(1),'YTickLabel',num2str(linspace(.59, .61,5)'))
%        set(AX(2),'YTick',linspace(-.015,.015,5))
%        set(AX(2),'YTickLabel',num2str(linspace(-.015, .015,5)'))
%        set(AX(1),'YGrid','on')
%        set(get(AX(1),'Ylabel'),'String','Pressure','FontWeight','demi','FontSize',10)
%        set(get(AX(2),'Ylabel'),'String','x-Velocity','FontWeight','demi','FontSize',10)
%        set(AX(1),'XLim',[xmin xmax])
%        set(AX(2),'XLim',[xmin xmax])
%        set(H1,'LineStyle','-.');%,'Marker','.')
       
%        subplot_title(['Acoustic Wave: t = ',num2str(t,'%4.3f')],[.1 .85 .4 .05]);
%        subplot_title(['Acoustic Wave: t = ',num2str(t,'%4.3f')],[.3 .42 .4 .05]);
	   
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
toc