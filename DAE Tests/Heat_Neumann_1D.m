function [Phi_t, u_t, Mass_Error, outTimes] = Heat_Neumann_1D(n_x)

%% Set up a line for space
%n_x = 100;    % Number of collocation points in space
geom.yMin = -3;    geom.yMax = 4;    geom.N = n_x;
aLine = SpectralLine(geom);
aLine.ComputeAll;
Interp = aLine.ComputeInterpolationMatrix((-1:0.02:1)',true);   % Ben shall know


%% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% Neumann system
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%% Normal tests
y = ones(n_x,1);
y(1:n_x/2) = 0;

% Normal vs standard differentiation
%(aLine.Ind.normal * aLine.Diff.Dy - [-1;1] .* aLine.Diff.Dy(aLine.Ind.bound,:)

% aLine.Ind.normal * aLine.Diff.Dy * y 
% [-1;1] .* aLine.Diff.Dy(aLine.Ind.bound,:) * y
all(aLine.Ind.normal * aLine.Diff.Dy * y  - [-1;1] .* aLine.Diff.Dy(aLine.Ind.bound,:) * y == 0);

%% Boundary matrix
% reshape(aLine.Diff.Dy(aLine.Ind.left | aLine.Ind.right,aLine.Ind.left),
% 1, [])    % First row
% aLine.Diff.Dy(aLine.Ind.right, aLine.Ind.left | aLine.Ind.right)  %
% Second row
D_bound = [aLine.Diff.Dy(aLine.Ind.left, aLine.Ind.left | aLine.Ind.right);  ...
            aLine.Diff.Dy(aLine.Ind.right, aLine.Ind.left | aLine.Ind.right)];
% Alternatively
% all(full(aLine.Diff.Dy(aLine.Ind.bound, aLine.Ind.bound) - D_bound == 0),'all')
D_bound = aLine.Diff.Dy(aLine.Ind.bound, aLine.Ind.bound);
D_bound_inv = D_bound^-1;

%% Inner matrix at bound
D_bound_inner = [aLine.Diff.Dy(aLine.Ind.left, ~aLine.Ind.bound); aLine.Diff.Dy(aLine.Ind.right, ~aLine.Ind.bound)];
% Alternatively
% D_bound_inner - aLine.Diff.Dy(aLine.Ind.left | aLine.Ind.right, ~aLine.Ind.bound)
D_bound_inner = aLine.Diff.Dy(aLine.Ind.left | aLine.Ind.right, ~aLine.Ind.bound);

%% To save time, we can also obtain the forward operator for the inner part
D_inner_to_bound = D_bound_inv * D_bound_inner;


%% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% Differential system
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

%% Time collocation points
tMax = 8;    n_t  = 100;
% Set up a time line where we are currently at 
ge.yMin  = 0;    ge.yMax = tMax;    ge.N = n_t;
TimeLine    = SpectralLine(ge);
outTimes = TimeLine.Pts.y;

% Integration weights
Int_Time    = TimeLine.ComputeIntegrationVector;
Interp_Time = TimeLine.ComputeInterpolationMatrix((-1:0.02:1)',true);

%% Additional spatial handlers
Lap = aLine.Diff.DDy;       % Second derivative
bound = aLine.Ind.bound;    % Indices of bound terms
y = aLine.Pts.y;            % Mesh

%% Initial condition
%u_0 = (y.^4).^(1/3);       % Experiment 1
u_0 = y.^2;                 % Experiment 2
c_0 = (aLine.Int * u_0);
u_0 = u_0 / c_0;            % Experiment 2
%u_0 = u_0 / c_0 + 1;       % Experiment 2b
u_0_inner = u_0(~bound);

%% Solve the PDE, the RHS of which is given below

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Stats','on');         % Reduced computation times 10 fold
              
tic
[~,Phi_t] = ode45(@rhs, outTimes, u_0_inner, opts);
toc

% It might be steep but it worked:
% 268639 successful steps
% 17914 failed attempts
% 1.71932e+06 function evaluations
% Elapsed time is 13.320438 seconds.

%% Now process the solution for it to also contain the domain
u_t = zeros(n_t,n_x);
%f_Times = [exp(-outTimes), exp(-outTimes)];        % Experiment 1
%f_Times = (4/(3*c_0)) * [ -(3^(1/3)) * exp(-outTimes), (4^(1/3)) * exp(-outTimes)];        % Experiment 1b
f_Times = (2.0/c_0) * [ -3.0 * exp(-outTimes), 4.0 * exp(-outTimes)];                       % Experiment 2
u_bound = (D_bound_inv * f_Times' - D_inner_to_bound * Phi_t')'; % This is considering the cancelation of sign due to normal * condition on f 
u_t(:,bound) = u_bound;
u_t(:,~bound) = Phi_t;
% Time 0 might have a different boundary [i.e., it's discontinuous]
u_t(1,bound) = u_0(bound);

% Some plots
%plot(y,u_t')
%surf(y,outTimes,u_t,'Edgecolor','none')
contourf(y,outTimes,u_t)
colormap('bone')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
%colorbar('Ticks',[min(u_t,[],'all'), -0.1, 0.3, 0.6, 0.8],...                      % Experiment 1
%colorbar('Ticks',[min(u_t,[],'all'), 0.1, 0.2, 0.27,0.9*max(u_t,[],'all')],...     % Experiment 1b
colorbar('Ticks',[min(u_t,[],'all'), 0.1, 0.25, 0.4, 0.945*max(u_t,[],'all')],...
          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'}, ...
          'TickLabelInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

fontsize(14,"points")
exportgraphics(figure(1), 'HN_1D_Solution_Contour.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)


%Mass_Error = 1-(aLine.Int * u_t')';                                        % Experiment 1
%c_1_s = 1 + (4/(3*c_0)) * (3^(1/3) + 4^(1/3) ) * (1-exp(-outTimes));       % Experiment 1b
%c_1_s = 1 + (8/(3*c_0)) * (3^(5/3) + 4^(5/3) ) * (1-exp(-outTimes));       % Experiment 1c?
c_1_s = exp(-outTimes);                                                    % Experiment 2
%c_1_s = exp(-outTimes) + 7;
%Mass_Error = c_1_s - (aLine.Int * u_t')';                                  % Experiment 1a
%Mass_Error = 1.0 - (aLine.Int * u_t')'./c_1_s;                             % Experiment 1 & 2b
Mass_Error = (c_1_s - (aLine.Int * u_t')') ./ max(c_1_s, 5e-3) ;            % Experiment 2
% plot((aLine.Int * u_t')')
% hold on
% plot(c_1_s)


%% Define ODE in time
function dudt = rhs(t,u_inner)
    % Compute Neumann vector
    f_bound_t = f_bound(t) .* [-1;1];
    % Compute u at bound
    u_bound = D_bound_inv * f_bound_t - (D_inner_to_bound * u_inner);

    % Consolidate full domain information
    u = zeros(n_x,1);
    u(bound) = u_bound;
    u(~bound) = u_inner;

    % Apply differential operator
    Lu = Lap * u;
    % Add source (only needed in experiment 2)
    %g = ( (40/9) - y.^2 ) .* (y.^2).^(1/3) * exp(-t) / c_0;
    g = -(y.^2 + 2.0 ) * exp(-t) / c_0;
    Lu = Lu + g;

    % Project to interior of domain
    dudt = Lu(~bound);
end

%% Boundary condition
function dydn = f_bound(t)
    % For mass conservation we avoid evaluating at x
    % dydn = [-exp(-t);exp(-t)];                                               % Experiment 1
    % dydn = (4/(3*c_0)) * [ (3^(1/3)) * exp(-t); (4^(1/3)) * exp(-t)];        % Experiment 1b
    dydn = (2.0/c_0) * [ 3.0 * exp(-t); 4.0 * exp(-t)];                        % Experiment 2
end

end