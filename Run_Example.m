%% Example of the NLCH with singular potentials

%%  Load box and convolution structures
% Load a SK: 10 × 10 box with 𝜀 = 10^-5
load('Singular_Kernels_Subs_20_epsB.mat', 'Newtonians')
%Newtonians

% Retrieve box and number of collocation points
aBox     = Newtonians.box;          % This is the domain
[N1, N2] = deal(aBox.N1, aBox.N2);
epsilon  = Newtonians.eps;          % This is ε

% For the selected box and value of  we select a convolution matrix based on a subdivision factor:
Newtonians.Level.Available          % Check available matrices
% Select convolution matrix
Conv = Newtonians.Level.n8;         % Name format: 'n' + factor
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
Interp = aBox.ComputeInterpolationMatrix(...
                             (-1:0.02:1)',(-1:0.02:1)',true,true);

% Potential and derivatives
F  = @(s) (1+s) .* log(1+s) + (1-s) .* log(1-s);
dF = @(s) log(1+s) - log(1-s);
ddF = @(s) 2 ./ ( 1 - s.^2 );

%% Create a time interval
tMax = 1.0; % Usually T = 5
n_t  = 100;   % often 100
% Set up a time line where we are currently at 
ge.yMin  = 0.0;    ge.yMax = tMax;    ge.N = n_t;
TimeLine = SpectralLine(ge);
outTimes = TimeLine.Pts.y;

%% Define initial condition
Y1 = 0.5 * sin( 3 * pi * Pts.y1_kv) .* cos( 3 * pi * Pts.y2_kv) + 0.25;
Y2 = max( abs(Pts.y1_kv - 0.5), abs(Pts.y2_kv-0.5));         % l_\infty norm
Y1(Y2 > 0.36) = 0.0;

h_a = 1e-1;
Molly  = @(y1,y2) fillmissing( (y1.^2 + y2.^2  < 1.0 ) .* exp( -1.0./( 1 - (y1.^2 + y2.^2 ) ) ) , 'Constant', 0);
Molly_a = @(y1,y2) h_a.^(-2) * Molly(y1/h_a, y2/h_a);
H_a = aBox.ComputeConvolutionMatrix(Molly_a,true);
phi_ic = 3 * (H_a * Y1);
aBox.plot(phi_ic);

%% Modify kernel? (optional)
ConvL = Conv + (0.025) * H_a;


%% Solve equation
% Select a kernel scaling
eta = 300;
% Now solve
Phi_t = NL_CH_Integrator_Simple(phi_ic,  aBox, Diff, Ind,  ConvL, Conv_D, dF, ddF, TimeLine, eta, 0.0 );
%Phi_t = NL_CH_Integrator_DAE( phi_ic,  aBox, Diff, Ind,  ConvL, Conv_D, eta, TimeLine);                % Use when DAE fails

% Compute energy
Energy = (Int * F(Phi_t)') - 0.5 * eta * Int * ( (ConvL * Phi_t' + Conv_D .* Phi_t') .* Phi_t');
% Compute approximation of C_* and \delta
C_Val_Approx = dF(Phi_t(end,:)') - eta * (ConvL * Phi_t(end,:)' + Conv_D .* Phi_t(end,:)');
%aBox.plot(C_Val_Approx);
Int * (C_Val_Approx)                % Approx C_*
1-max( abs(Phi_t(end,:)),[],'all')  % Approx \delta

% Visualise
aBox.plot(Phi_t(end,:)');
% Animation
plots_in_box(Phi_t, aBox, outTimes,'MWE.gif', '$\rho(x;t)$')