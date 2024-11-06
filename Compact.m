%% Study of scalings of K
% Initial condition given by a function of compact support

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

% CASE C
Y1 = 0.5 * sin( 3 * pi * Pts.y1_kv) .* cos( 3 * pi * Pts.y2_kv) + 0.25;
Y2 = max( abs(Pts.y1_kv - 0.5), abs(Pts.y2_kv-0.5));         % l_\infty norm
Y1(Y2 > 0.36) = 0.0;


h_a = 1e-1;
Molly  = @(y1,y2) fillmissing( (y1.^2 + y2.^2  < 1.0 ) .* exp( -1.0./( 1 - (y1.^2 + y2.^2 ) ) ) , 'Constant', 0);
Molly_a = @(y1,y2) h_a.^(-2) * Molly(y1/h_a, y2/h_a);
H_a = aBox.ComputeConvolutionMatrix(Molly_a,true);
phi_ic = 3 * (H_a * Y1);
aBox.plot(phi_ic);

eta = 1.0;


%% Solve equation iteratively
Ens  = struct('A',[], 'B',[], 'C',[], 'D',[]); %struct('A',[], 'B',[], 'C',[], 'D',[], 'E',[], 'F',[]);
etas = [500, 100, -50, -100];
Fields = fieldnames(Ens);

% Iterate over each field
for i = 1:numel(Fields)
    eta = etas(i);
    % Solve equation
    %Phi_to = NL_CH_Integrator_DAE( phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, eta, TimeLine);
    % If process gets stuck (-150), then try:
    Phi_to = NL_CH_Integrator_Simple(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, TimeLine, eta, 0.0 );
    % Energy
    Energy = (Int * F(Phi_to)') - 0.5 * eta * Int * ( (Conv * Phi_to' + Conv_D .* Phi_to') .* Phi_to');
    Ens.(Fields{i}) = Energy;

    % Plot solution as an animation
    Out_Name = strcat(['C_eta=',int2str(eta), '.gif']);
    plots_in_box(Phi_to, aBox, outTimes, Out_Name, '$\rho(x;t)$')

    % Equilibria
    C_Val_Approx = dF(Phi_to(end,:)') - eta * (Conv * Phi_to(end,:)' + Conv_D .* Phi_to(end,:)');
    %aBox.plot(C_Val_Approx);
    max(C_Val_Approx) - (max(C_Val_Approx) + min(C_Val_Approx))/2

    % Modify plot_panels for this:
    Out_Name = strcat(['C',int2str(i), '_Solution_eta=',int2str(eta), '.pdf']);
    plot_panels
end

%% Plot struct
%% Plotting: up to time t=0.05 for graphical purposes


% Colour selection
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:7),7);

% Plot for each field in struct
Lstyles = {'-','--','-.','-','--','-.','-'};
for i = 1:numel(Fields)
    eta = etas(i);
    plot(0:0.05/100:0.05, interp1(outTimes,Ens.(Fields{i}),0:0.05/100:0.05,'makima'), 'LineWidth', 1.5, ...
                        'LineStyle', Lstyles{i}, 'DisplayName', strcat(['$',int2str(eta), '$']) )
    hold on
end

colororder(G_colours)

xlabel('$t$','Interpreter','latex'); 
ylabel( '$\mathcal{E}_\eta(t)$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
fontsize(16, "points")
%set(gca, 'XScale', 'log')
lgd = legend('show', 'Interpreter', 'latex', 'Location','southeast');
title(lgd, 'Scaling $\eta$', 'Interpreter', 'latex'); % Title for the legend


set(gca, 'FontName', 'CMR10')


exportgraphics(figure(1), 'NLCH_Energy_C.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
hold off


%% Additional plots (optional): These plots can be added to the for loop above to obtain the gallery for this example
%plots_in_box(Phi_to, aBox, outTimes,'A_eta=-1.gif', '$\rho(x;t)$')
% Integrate
Mass_Time = Int * Phi_to';
plot(outTimes, Mass_Time, 'LineWidth', 1.5, 'Color', colores(1,:), 'LineStyle', '-', 'Marker','x')


%%
