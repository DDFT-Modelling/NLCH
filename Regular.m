%% Study the approximation through regular potentials
% Initial condition is a wave of compact support

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


% Regular potential
dddF = @(s) 4*s ./ (s.^2 -1 ).^2;

h_b = 1e-3;
f_ra = @(s) dF( 1-h_b) +  ddF(1-h_b) * (s - 1 + h_b) + 0.5 * dddF( 1-h_b) * (s - 1 + h_b ).^2;
f_rb = @(s) dF(-1+h_b) + ddF(-1+h_b) * (s + 1 - h_b) + 0.5 * dddF(-1+h_b) * (s + 1 - h_b).^2;

dF_r  = @(s) (s <= -1+h_b ) .* f_rb(s) + (s >= 1-h_b ) .* f_ra(s)  +  dF( ( ( -1 + h_b < s ) & ( 1 - h_b > s ) ) .* s);
ddF_r = @(s) (s <= -1+h_b ) .* ( ddF(-1+h_b) + dddF(-1+h_b) * (s + 1 - h_b) ) + ...
             (s >=  1-h_b ) .* ( ddF(1-h_b) + dddF( 1-h_b) * (s - 1 + h_b ) )  +  ddF( ( ( -1 + h_b < s ) & ( 1 - h_b > s ) ) .* s);


%% Solve equation 

% Use the same initial condition as in case C and set eta to 500 (set all tolerances to 1e-11)
eta = 500;
% Simulate things up to Tmax = 1/3 * 2e-3
Phi_t = NL_CH_Integrator_Simple(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, TimeLine, eta, 0.0 );
% Use final evaluation as initial condition
phi_ic_r = Phi_t(end,:)';
% Also solve for Tmax = 1/3 * 1e-2
Phi_t = NL_CH_Integrator_Simple(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, TimeLine, eta, 0.0 );
phi_T = Phi_t(end,:)';

% Plot solution as an animation
Out_Name = strcat(['D_eta=',int2str(eta), '.gif']);
plots_in_box(Phi_t, aBox, 10*outTimes, Out_Name, '$\rho(x;t)$')
% Modify plot_panels for this:
Out_Name = strcat(['D',int2str(i), '_Solution_eta=',int2str(eta), '.pdf']);
Phi_to = Phi_t;
plot_panels


% Run regularised problem up to time 1/375
Phi_t = NL_CH_Integrator_Simple(phi_ic_r,  aBox, Diff, Ind,  Conv, Conv_D, dF_r, ddF_r, TimeLine, eta, 0.0 );


% Error:
(Int * (Phi_t(end,:)' - phi_T).^2)^0.5

