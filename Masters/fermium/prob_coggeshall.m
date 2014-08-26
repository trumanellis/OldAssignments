% Acoustic Wave Problem
% Author: Truman Ellis
% Email: ellis35@llnl.gov

clear all; close all; clc;
format long g
format compact
scrsz = get(0,'ScreenSize');
h0=sfigure(1,'Position',[1 scrsz(4)/8 scrsz(3) scrsz(4)*7/9]);
set(h0,'renderer','painters')
CoggeshallTransform=true;
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
    VBasis=@Q2Basis;
    PBasis=@Q1dBasis;
    HighOrderEnergy=false;
    StrongMass=false;
    
    FullMassMatrixSolve=false;
    MassUpdate=false;
    % Set Quadrature Order
    QuadOrder=3;

    % Use this to turn artificial viscosity on/off
    Qfrac=1;

    % Use this to turn anti-hourglass forces on/off (Only works with Q1 elements)
    hgfrac=0;

    % Define the number of zones
    NZx=10;
    NZy=4;

    % Define the mesh size
    xmin=.1;
    xmax=1;

    ymin=0;
    ymax=pi/2;

    % Define the start and stop time
    tstart = 0;
    tstop = 0.7;
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
    qquad=2;
    qlin=0.25;

    if SaveFigures
        figurefile=['Coggeshall',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
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
    eInit=3/25;
    pInit=(gamma-1)*rhoInit*eInit;

    origin=find(sqrt(allnodes(1,:).^2+allnodes(2,:).^2) == xmin);
    
    for N=1:NZ
        for n=1:npdof
            xglob=LocalToGlobal(allnodes(1,Quadmap(:,N)),pdofpts(:,n),VBasis);
            yglob=LocalToGlobal(allnodes(2,Quadmap(:,N)),pdofpts(:,n),VBasis);
            OLDdensity(n,N)=rhoInit*sqrt(xglob^2+yglob^2);
        end
        xglob=LocalToGlobal(allnodes(1,Quadmap(:,N)),[.5;.5],VBasis);
        yglob=LocalToGlobal(allnodes(2,Quadmap(:,N)),[.5;.5],VBasis);
        OLDenergy(1,N)=eInit*(xglob^2+yglob^2);
        pressure(:,N)=(gamma-1)*OLDdensity(:,N).*OLDenergy(:,N);
        % Evaluate density and pressure at quadrature points
        QuadDensity(:,N)=(OLDdensity(:,N)'*QuadPBasis)';
        QuadPressure(:,N)=(pressure(:,N)'*QuadPBasis)';
    end

    OLDvelocity(1,:)=-3/5*sqrt(allnodes(1,:).^2+allnodes(2,:).^2).*allnodes(1,:)./sqrt(allnodes(1,:).^2+allnodes(2,:).^2);
    OLDvelocity(2,:)=-3/5*sqrt(allnodes(1,:).^2+allnodes(2,:).^2).*allnodes(2,:)./sqrt(allnodes(1,:).^2+allnodes(2,:).^2);
    
    % Impose Boundary Conditions
    % Impose symmetric boundary at bottom wall
    for i=1:sum(nbnodes(1))
        index=bdof(i);
        OLDvelocity(2,index)=0;
    end
%     % Impose zero boundary at right wall
%     for i=1+nbnodes(1):sum(nbnodes(1:2))
%         index=bdof(i);
%         OLDvelocity(:,index)=[0;0];
%     end
    % Impose symmetric boundary at top wall
    for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
        index=bdof(i);
        OLDvelocity(1,index)=0;
    end
    % Impose zero boundary at left wall
    for i=1+sum(nbnodes(1:3)):sum(nbnodes)
        index=bdof(i);
        OLDvelocity(:,index)=[0;0];
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
    
    % Set Driver
%     drivingPt=[intersect(find(allnodes(1,:) == 0),find(allnodes(2,:) == 0))];
%     amp=sqrt(1e-6/massN(drivingPt));
%     omega=2*2*pi;
%     if exist('RotateMesh') && RotateMesh > 0
%         rot=RotMat*[1;0];
%     else
%         rot=[1;0];
%     end
%     drivingFnx=@(a,b,t) amp*rot(1)*sin(omega*t);
%     drivingFny=@(a,b,t) amp*rot(2)*sin(omega*t);
%     OLDvelocity(1,drivingPt(1))=drivingFnx([],[],tstart);
%     OLDvelocity(2,drivingPt(1))=drivingFny([],[],tstart);
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
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            forceN(1,index)=0;
        end
        % Impose zero boundary at left wall
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
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            acceleration(1,index)=0;
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
       
       plotvarN=sqrt(NEWvelocity(1,:).^2+NEWvelocity(2,:).^2);
               
       subplot(2,2,1)
       [h1,X2pts,Nvals]=PlotFEMContourf(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
       PlotFEMMesh(allnodes,Quadmap,10,VBasis,'k');
       axis tight
%        axis([0 .15 0 .15])
       daspect([1 1 1e-1]);
       title('Velocity Magnitude','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
               
       subplot(2,2,3)
       presplot=NEWdensity(:)';
       if ~StrongMass
           [h3,X1pts,Zvals]=PlotFEMContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
       else
           [h3,X1pts,Zvals]=PlotStrongMassContourf(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
       end
       PlotFEMMesh(allnodes,Quadmap,10,VBasis,'k');
       axis tight
%        axis([0 .15 0 .15])
       daspect([1 1 1e-1]);
       title('Density','FontWeight','demi','FontSize',14)
       view(az,el);
       colorbar
               
       lineoutX=sqrt(X1pts(1,:).^2+X1pts(2,:).^2);
       lineoutZ=Zvals;
       
       subplot(2,2,2)
       plot(lineoutX,lineoutZ,'r.')
       title('Density','FontWeight','demi','FontSize',14)

%        [AX,H1,H2]=plotyy(lineoutX,lineoutZ,lineoutX2,lineoutN);
%        set(AX(1),'YLim',[0 5])
%        set(AX(2),'YLim',[0 4])
%        set(AX(1),'YTick',linspace(0,5,5))
%        set(AX(1),'YTickLabel',num2str(linspace(0, 5,5)'))
%        set(AX(2),'YTick',linspace(0,4,5))
%        set(AX(2),'YTickLabel',num2str(linspace(0, 4,5)'))
%        set(AX(1),'YGrid','on')
%        set(get(AX(1),'Ylabel'),'String','Density','FontWeight','demi','FontSize',10)
%        set(get(AX(2),'Ylabel'),'String','Velocity Magnitude','FontWeight','demi','FontSize',10)
%        set(AX(1),'XLim',[0 1])
%        set(AX(2),'XLim',[0 1])
%        set(H1,'LineStyle','--','Marker','.')
       if userno == 7 || userno == 12 || userno == 15 || userno == 10
           subplot(2,2,4)
%            enrgplot=NEWenergy(:)';
%            if HighOrderEnergy
%                EMap=PressureMap';
%                EBasis=PBasis;
%            else
%                EMap=1:NZ;
%                EBasis=@Q0Basis;
%            end
%            [h4]=PlotFEMContourf(allnodes,enrgplot,Quadmap,EMap,pref,VBasis,EBasis,'Linestyle','none');
%            axis tight
%            daspect([1 1 1e-2]);
%            title('Energy','FontWeight','demi','FontSize',14)
%            view(az,el);
%            colorbar
           [h3,X1pts,Dvals]=PlotFEMContourf(allnodes,phi0Z,Quadmap,1:NZ,pref,VBasis,@Q0Basis,'Linestyle','none');
           PlotFEMMesh(allnodes,Quadmap,10,VBasis,'w');
           view(az,el);
           set(h3,'LineStyle','none')
           axis equal tight
%            axis([xmin xmax ymin ymax 0 1 0 1])
           title('phi0','FontWeight','demi','FontSize',14)
           colorbar
       else
           subplot(2,2,4)
           plot(timedata(1:cycle),internalenergy(1:cycle)+kineticenergy(1:cycle),'k',...
               timedata(1:cycle),internalenergy(1:cycle),'r',...
               timedata(1:cycle),kineticenergy(1:cycle),'b')
           legend('Total Energy','Internal Energy','Kinetic Energy','Location','NorthWest')
           axis([0,timedata(cycle),0 1.4*(internalenergy(cycle)+kineticenergy(cycle))])
       end
       
       subplot_title(['Sedov: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
       
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