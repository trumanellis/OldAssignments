% Written by Truman Ellis
clear all; close all; clc;

% Number of beta's to plot
nBeta = 50;
beta = logspace(-2, 1, nBeta);
% Number of Integration Points
nx = 50;
ny = 50;
% Integration Domain
xMin = -5;
xMax = 3;
yMin = -3;
yMax = 5;
% Calculate grid size
hx = (xMax - xMin)/nx;
hy = (yMax - yMin)/ny;
dA = hx*hy;
% Construct mesh
x = linspace(-5, 3, nx);
y = linspace(-3, 5, ny);

% Define offset to make integration easier
Umin = -145;

% Start quadrature sum
dU = zeros(1,nBeta);
% Loop over different beta's
for nb = 1:nBeta    
    % Integrate over the domain
    for i = 1:nx
        for j=1:ny
            % Sum differential volume
           dU(nb) = dU(nb) + exp(-beta(nb)*(muller(x(i), y(j))-Umin))*dA; 
        end
    end 
end
% Calculate free energy: F = -1/beta*log(Z)
F = Umin + -1./beta.*log(dU);

% Plot results
plot(beta,F,'o','MarkerFaceColor','blue','MarkerSize',5)
xlabel('\beta','FontSize',14)
ylabel('Free Energy','FontSize',14,'FontWeight','demi')
