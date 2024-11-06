%% Study of scalings of K
% Initial condition given by a periodic wave ρ₀(x) = sin(2π x₁) cos(2π x₂)

%% Load data
load('Singular_Kernels_Subs_40_epsB.mat', 'Newtonians')
%Newtonians

% Retrieve box and number of collocation points
aBox     = Newtonians.box;          % This is the domain
[N1, N2] = deal(aBox.N1, aBox.N2);
epsilon  = Newtonians.eps;          % This is ε

% For the selected box and value of  we select a convolution matrix based on a subdivision factor:
Newtonians.Level.Available          % Check available matrices
% Select convolution matrix
Conv = Newtonians.Level.n4;         % Name format: 'n' + factor
% Retrieve additional data
Conv_D = Newtonians.NG;

% Define additional structures in the box for defining the differential equation:
% Points, differentiation matrices, integration vector, and indices
% giving masks for the boundary
[Pts,Diff,Int,Ind] = aBox.ComputeAll();    
grad  = Diff.grad;       div = Diff.div;
bound = Ind.bound;    normal = Ind.normal;
Lap_Space = Diff.Lap;

% compute the interpolation matrix from the box points to a uniform
% grid [only for plotting]
Interp = aBox.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);

% Potential and derivatives
F  = @(s) (1+s) .* log(1+s) + (1-s) .* log(1-s);
dF = @(s) log(1+s) - log(1-s);
ddF = @(s) 2 ./ ( 1 - s.^2 );


%% Setup
colores = [255, 200, 87; 186, 45, 11; 86, 22, 67]/255;

% Time interval
tMax = 1.0; % Usually T = 5
n_t  = 100;   % often 100
% Set up a time line where we are currently at 
ge.yMin  = 0.0;    ge.yMax = tMax;    ge.N = n_t;
TimeLine = SpectralLine(ge);
outTimes = TimeLine.Pts.y;
TimeLine.ComputeDifferentiationMatrix;

% CASE A
phi_ic = sin( 2.0 * Pts.y1_kv * pi ) .* cos( 2.0 * Pts.y2_kv * pi );


%% Solve equation iteratively
Ens  = struct('A',[], 'B',[], 'C',[], 'D',[]);
etas = [-1e+2, -1.5e+2, -5e+1, 1.0];
Fields = fieldnames(Ens);

% Iterate over each field
for i = 1:numel(Fields)
    eta = etas(i);
    % Solve equation
    Phi_to = NL_CH_Integrator_DAE( phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, eta, TimeLine);
    % Energy
    Energy = (Int * F(Phi_to)') - 0.5 * eta * Int * ( (Conv * Phi_to' + Conv_D .* Phi_to') .* Phi_to');
    Ens.(Fields{i}) = Energy;

    % Equilibria
    C_Val_Approx = dF(Phi_to(end,:)') - eta * (Conv * Phi_to(end,:)' + Conv_D .* Phi_to(end,:)');
    aBox.plot(C_Val_Approx);
    max(C_Val_Approx) - (max(C_Val_Approx) + min(C_Val_Approx))/2
end

%% Plot struct
%% Plotting: up to time t=0.25 for graphical purposes


% Colour selection
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:7),7);

%
plot(0:0.2/100:0.25, interp1(outTimes,Ens.D,0:0.2/100:0.25,'makima'), 'LineWidth', 1.5, 'Color', 'black', 'LineStyle', ':', 'DisplayName','$+1$')
hold on
plot(0:0.2/100:0.25, interp1(outTimes,Ens.C,0:0.2/100:0.25,'makima'), 'LineWidth', 1.5, 'Color', colores(2,:), 'LineStyle', '--', 'DisplayName','$-50$')
plot(0:0.2/100:0.25, interp1(outTimes,Ens.A,0:0.2/100:0.25,'makima'), 'LineWidth', 1.5, 'Color', colores(3,:), 'LineStyle', '-.', 'DisplayName','$-100$')
plot(0:0.2/100:0.25, interp1(outTimes,Ens.B,0:0.2/100:0.25,'makima'), 'LineWidth', 1.5, 'Color', colores(1,:), 'DisplayName','$-150$')


%colororder(G_colours)

xlabel('$t$','Interpreter','latex'); 
ylabel( '$\mathcal{E}_\eta(t)$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
fontsize(16, "points")
%set(gca, 'XScale', 'log')
lgd = legend('show', 'Interpreter', 'latex', 'Location','southeast');
title(lgd, 'Scaling $\eta$', 'Interpreter', 'latex'); % Title for the legend
%ylim([-2.5 0.5])
%xlim([0 0.2]);

set(gca, 'FontName', 'CMR10')


exportgraphics(figure(1), 'NLCH_Energy_PW.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
hold off


%% Additional plots (optional): These plots can be added to the for loop above to obtain the gallery for this example
plots_in_box(Phi_to, aBox, outTimes,'A_eta=-1.gif', '$\rho(x;t)$')
% Integrate
Mass_Time = Int * Phi_to';
plot(outTimes, Mass_Time, 'LineWidth', 1.5, 'Color', colores(1,:), 'LineStyle', '-', 'Marker','x')


%%
