,% Sod Shock Tube Problem
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
recover=false;%'FiguresSaltzman7_2009_08_03_16_11/save0.1199.mat';
sodExact=load('sod_solution.dat');
if recover
    load([SaveLocation,recover]);
    cycle=cycle+1;
    tic
else
    tic

    SaveFigures=false;
    nplots=10;
    nsaves=10;
    % Set initial time increment
    dtInit=5e-4;
    dtMaxFinal=1e-4;

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
    NZx=40;
    NZy=1;

    % Define the mesh size
    xmin=0;
    xmax=1;

    ymin=0;
    ymax=.1*xmax;

    % Define the start and stop time
    tstart = 0;
    tstop = .2;
%     tplot=sort(tstop-logspace(log10(.1*tstop/nplots),log10(tstop),nplots)+.1*tstop/nplots);
%     tsave=sort(tstop-logspace(log10(.1*tstop/nsaves),log10(tstop),nsaves)+.1*tstop/nsaves);
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
    % Define Saltzman mesh perturbation
    SaltzmanPerturb=0;
    
    vref=5; pref=5;

    % Specify the maximum number of time steps to use
    maxcycle=1e7;
    dtMax=dtMaxFinal;

    % These numbers are the coefficients for the artificial viscosity
    qquad=1.0;
    qlin=1.0;

    if SaveFigures
        figurefile=['Sod',func2str(VBasis),'-',func2str(PBasis),'_',datestr(clock, 'yyyy_mm_dd_HH_MM')];
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

    rhoInitLHS=1;
    pInitLHS=1;
    eInitLHS=pInitLHS/(gamma-1)/rhoInitLHS;
    rhoInitRHS=.125;
    pInitRHS=.1;
    eInitRHS=pInitRHS/(gamma-1)/rhoInitRHS;
    
    rhs=find(allnodes(1,Quadmap(1,:)) >= .5);

    OLDenergy(:)=eInitLHS;
    OLDdensity(:)=rhoInitLHS;
    QuadDensity(:)=rhoInitLHS;
    pressure(:)=pInitLHS;
    QuadPressure(:)=pInitLHS;
    
    OLDenergy(rhs)=eInitRHS;
    OLDdensity(:,rhs)=rhoInitRHS;
    QuadDensity(:,rhs)=rhoInitRHS;
    pressure(:,rhs)=pInitRHS;
    QuadPressure(:,rhs)=pInitRHS;
    
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
        % Impose symmetric boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            acceleration(1,index)=0;
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            acceleration(2,index)=0;
        end
        % Impose zero boundary at left wall
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

       plotvarN=NEWvelocity(1,:);
       [X2pts,Nvals]=PlotFEM1line(allnodes,plotvarN,Quadmap,Quadmap,vref,VBasis,VBasis,'Linestyle','none');
       presplot=NEWdensity(:)';
       [X1pts,Zvals]=PlotFEM1line(allnodes,presplot,Quadmap,PressureMap',pref,VBasis,PBasis,'Linestyle','none');
       
       lineoutX=X1pts(1,:);
       lineoutZ=Zvals;
       lineoutX2=X2pts(1,:);
       lineoutN=Nvals(1,:);
%        lineoutX3=X3pts(1,:);
%        lineoutpsi=psivals(1,:);

       sfigure(1);
%        plot(lineoutX2,lineoutN,lineoutX,lineoutZ,lineoutX3,lineoutpsi);
       plot(lineoutX2,lineoutN,lineoutX,lineoutZ);
%        legend('Velocity','Density','psi_0')
       
           
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