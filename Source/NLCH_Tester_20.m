%% Testing the NL CH code for different boxes and values of ε


N_box = 20;
load(strcat(['Singular_Kernels_Subs_',int2str(N_box),'_epsB.mat']), 'Newtonians')

%% Retrieve box and number of collocation points
aBox     = Newtonians.box;          % This is the domain
[N1, N2] = deal(aBox.N1, aBox.N2);
epsilon  = Newtonians.eps;          % This is ε

% Define additional structures from box
[Pts,Diff,Int,Ind] = aBox.ComputeAll();    
grad  = Diff.grad;       div = Diff.div;
bound = Ind.bound;    normal = Ind.normal;
Lap_Space = Diff.Lap;

% compute the interpolation matrix from the box points to a uniform
% grid [only for plotting]
Interp = aBox.ComputeInterpolationMatrix( (-1:0.02:1)',(-1:0.02:1)',true,true);

%% Exact solution
phi_ic = sin( 2.0 * Pts.y1_kv * pi ) .* cos( 2.0 * Pts.y2_kv * pi );

% Time collocation points
tMax = 10.0; % Usually T = 10
n_t  = 70;   % often 100
% Set up a time line where we are currently at 
ge.yMin  = 0;    ge.yMax = tMax;    ge.N = n_t;
TimeLine = SpectralLine(ge);
outTimes = TimeLine.Pts.y;
TimeLine.ComputeDifferentiationMatrix;
D_Time = TimeLine.Diff.Dy;
phi = exp(-outTimes) * phi_ic';

% W function
F  = @(s) (1+s) .* log(1+s) + (1-s) .* log(1-s);
dF = @(s) log(1+s) - log(1-s);
ddF = @(s) -2 ./ ( 1 - s.^2 );

% Retrieve integral of kernel at excluded region
Conv_D = Newtonians.NG;

%% Iterate on available factors
Fields = fieldnames(Newtonians.Level);
for i = 3:numel(Fields)
    % Perform operation on the field
    fprintf('Processing field %s\n', Fields{i});

    % Select convolution matrix
    Conv = Newtonians.Level.(Fields{i});
    
    % For example, calculate the mean
    meanValue = mean(Conv ,'all');
    fprintf('Mean of %s: %.2e\n', Fields{i}, meanValue);

    g_tilde = -phi;
    mu_tilde = dF(phi);    % Notice that time zero is not well defined
    mu_tilde = mu_tilde - (Conv * phi' + Conv_D .* phi')';
    g_tilde = g_tilde - (Lap_Space * mu_tilde')' ;

    % Solve PDE
    %Phi_t = NL_CH_Integrator(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, 0,0, TimeLine );
    Phi_t = NL_CH_Integrator_Source(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, mu_tilde, g_tilde, TimeLine );


    % Plot energy
    Time_Span_ODE = [0: outTimes(end)/100: outTimes(end)];
    
    Mass_Time = Int * Phi_t';
    plot(Time_Span_ODE, Mass_Time, 'LineWidth', 1.5, 'Color', '#e9c46a', 'LineStyle', '-', 'Marker','x')
    ylim([min(Mass_Time)*1.01, max(Mass_Time)*1.2])
    
    xlabel('$t$','Interpreter','latex'); 
    ylabel( strcat(['$\int\limits_\Omega \varphi(x;t) ', ...
        '\mathrm{d}x $']), ...
        'Interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    
    Out_Name = strcat(['NLCH_Diffusion_Test_Mass_N=',int2str(N_box), '_factor_',  Fields{i}, '.pdf']);
    exportgraphics(figure(1), Out_Name,'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
    
    close all


    % Plot solution
    Phi_e = Phi_t * 0.0;

    k = 1;
    for t = Time_Span_ODE
        IP_Time = TimeLine.ComputeInterpolationMatrixPhys(t).InterPol;
        Phi_e(k,:) = IP_Time * phi;
        k = k + 1;
    end
    
    Out_Name = strcat(['NLCH_Diffusion_Test_Sol_N=',int2str(N_box), '_factor_',  Fields{i}, '.gif']);
    plots_in_box(Phi_t, aBox, Time_Span_ODE, Out_Name, '$\varphi(x;t)$')
    close all

    % Plot error
    Out_Name = strcat(['NLCH_Diffusion_Test_Error_N=',int2str(N_box), '_factor_',  Fields{i}, '.gif']);
    plots_in_box(Phi_e - Phi_t, aBox, Time_Span_ODE, Out_Name, '$e(x;t)$')
    close all


    % Plot l_2 error in time
    Err_Int = (Int * ((Phi_t - Phi_e)'.^2)).^0.5;
    
    %figure(1)
    plot(Time_Span_ODE, Err_Int, 'LineWidth', 1.5, 'Color', ...
        '#2a9d8f', 'LineStyle', '-', 'Marker','x')
    ylim([0,max(Err_Int)*1.01])
    
    xlabel('$t$','Interpreter','latex'); 
    ylabel( strcat(['$\left( \int\limits_\Omega [\varphi(x;t) ', ...
        '- \tilde{\varphi}(x;t)]^2 \mathrm{d}x \right)^{1/2} $']), ...
        'Interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    
    Out_Name = strcat(['NLCH_Diffusion_Test_Error_l2_N=',int2str(N_box), '_factor_',  Fields{i}, '.pdf']);
    exportgraphics(figure(1), Out_Name, 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
    close all

end
