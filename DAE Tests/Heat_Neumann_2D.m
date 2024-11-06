function [Phi_t, u_t, Mass_Error, outTimes] = Heat_Neumann_2D(n_x)

%% Set up a box for space
%n_x = 30;    % Number of collocation points in space
geom.N = [n_x;n_x];
geom.y1Min = -0.25; geom.y1Max = 0.5; geom.y2Min = -0.5; geom.y2Max = 0.25;
aBox = Box(geom);
aBox.ComputeAll;

Interp = aBox.ComputeInterpolationMatrix( (-1:0.02:1)',(-1:0.02:1)',true,true);
% aBox.plot(zeros(30*30,1))

%% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% Neumann system
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%% Normal tests
y = ones(n_x^2,1);
% y(1:n_x^2/2) = 0;

% Normal vs standard differentiation
normal_grad = aBox.Ind.normal * aBox.Diff.grad;
% cond(full(normal_grad(:,aBox.Ind.bound))) ~ 2.92
% normal_grad(:,~aBox.Ind.bound) * y(~aBox.Ind.bound) + normal_grad(:,aBox.Ind.bound) * y(aBox.Ind.bound) ~ 1e-12
% -normal_grad(:,aBox.Ind.bound) \ ( normal_grad(:,~aBox.Ind.bound) * y(~aBox.Ind.bound) )
norm(-normal_grad(:,aBox.Ind.bound) \ ( normal_grad(:,~aBox.Ind.bound) * y(~aBox.Ind.bound) ) - 1)    % ~1e-15
% How well bahaved is this system?
% NGb = normal_grad(:,aBox.Ind.bound);
% all(eig(full(NGb)) > 0)   is true
% cond(full(NGb)) ~ 2.9216
% The system is not diagonally dominant by rows but by columns:
% all( full( (2 * abs(diag(NGb))) >= sum(abs(NGb),1)' ) )   is true

% Let us consider a discontinuous function
y(1:n_x^2/2) = 0;

% We want to equate the normal derivative of a solution to the differential
% operator to our data, so let's compute the normal derivate of y
dy_dn = aBox.Ind.normal * aBox.Diff.grad * y;
% So this will be our RHS but let's assume we have only observed y at the
% interior:
y_reconstruct = normal_grad(:,aBox.Ind.bound) \ ( dy_dn - normal_grad(:,~aBox.Ind.bound) * y(~aBox.Ind.bound) );
% Since we have the exact solution at hand, we can compute the exact error
norm(y_reconstruct - y(aBox.Ind.bound))   % ~1e-16

% Let us try once more for a more complicated but smooth solution
y = exp(-aBox.Pts.y1_kv.^2 - aBox.Pts.y2_kv.^2);
dy_dn = aBox.Ind.normal * aBox.Diff.grad * y;
y_reconstruct = normal_grad(:,aBox.Ind.bound) \ ( dy_dn - normal_grad(:,~aBox.Ind.bound) * y(~aBox.Ind.bound) );
norm(y_reconstruct - y(aBox.Ind.bound))   % ~1e-15

%% Boundary matrix
D_bound_inv = normal_grad(:,aBox.Ind.bound)^-1;   % Also: (aBox.Ind.normal * aBox.Diff.grad(:,aBox.Ind.bound))^-1
%% To save time, we can also obtain the forward operator for the inner part
D_inner_to_bound = D_bound_inv * normal_grad(:,~aBox.Ind.bound);


%% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% Differential system
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

%% Time collocation points
tMax = 3;    n_t  = 100;
% Set up a time line where we are currently at 
ge.yMin  = 0;    ge.yMax = tMax;    ge.N = n_t;
TimeLine = SpectralLine(ge);
outTimes = TimeLine.Pts.y;

%% Additional spatial handlers
Lap = aBox.Diff.Lap;        % Laplacian
bound = aBox.Ind.bound;     % Indices of bound terms
normal = aBox.Ind.normal;   % Normal operator
y_1 = aBox.Pts.y1_kv;       % Horizontal mesh
y_2 = aBox.Pts.y2_kv;       % Vertical mesh

%% Initial condition
%u_0 = exp(-y_1 - y_2) / 3.0;                       % Experiment 1
c = 50.0;                                           % Experiment 2
static = exp(-c*y_1.^2 - c*y_2.^2);                 % Experiment 2
u_0 = ((2*exp(1)/(1+exp(2)))^2) * static;           % Experiment 2
u_0_inner = u_0(~bound);
%% Precompute f without time dependance
%dy = [ -exp(-y_1 - y_2); -exp(-y_1 - y_2) ]/3;     % Experiment 1: exp(-t)
dy = -2 * c * [y_1; y_2] .* [static;static];        % Experiment 2
dy_dn = normal * dy;        % Lacks exp(-t)

%% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% PDE Solver
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%% Solve the PDE, the RHS of which is given below

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Stats','on');
              
tic
[~,Phi_t] = ode45(@rhs, outTimes, u_0_inner, opts);
toc

% It might be steep but it worked:
% 137279 successful steps
% 9152 failed attempts
% 878587 function evaluations
% Elapsed time is 26.892985 seconds.

%% Now process the solution for it to also contain the domain
u_t = zeros(n_t,n_x*n_x);
%f_Times = dy_dn * exp(-outTimes)';                   % Experiment 1
f_Times = dy_dn * (1.0 - tanh(outTimes-1.0).^2)';     % Experiment 2

u_bound = (D_bound_inv * f_Times - D_inner_to_bound * Phi_t')';
u_t(:,bound) = u_bound;
u_t(:,~bound) = Phi_t;
% Time 0 might have a different boundary [i.e., it's discontinuous]
u_t(1,bound) = u_0(bound);

% At last define an error
%plot(outTimes, aBox.Int * u_t')
Mass = aBox.Int * u_t';
%Exact_Mass = (aBox.Int * u_t(1,:)') * exp(-outTimes)';                 % Experiment 1
Exact_Mass = (aBox.Int * static) * (1.0 - tanh(outTimes-1.0).^2)';      % Experiment 2
Mass_Error = 1 - Mass./Exact_Mass;

%plot(outTimes, Mass_Error);


%% Define ODE in time
function dudt = rhs(t,u_inner)
    % Compute RHS of Neumann condition
    %f_bound_t = dy_dn * exp(-t);                   % Experiment 1
    f_bound_t = dy_dn * (1.0 - tanh(t-1.0).^2);     % Experiment 2
    % Compute u at bound
    u_bound = D_bound_inv * f_bound_t - (D_inner_to_bound * u_inner);

    % Consolidate full domain information
    u = zeros(n_x*n_x, 1);
    u(bound) = u_bound;
    u(~bound) = u_inner;

    % Apply differential operator
    %Lu = Lap * u - 3 * u;      % Experiment 1

    g = -4 * c^2 * (y_1.^2 + y_2.^2) - 2 * tanh(t-1) + 4*c;
    Lu = Lap * u + g .* u;      % Experiment 2

    % Project to interior of domain
    dudt = Lu(~bound);
end

end