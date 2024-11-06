%% Storing 20s          -       Fields = fieldnames(Newtonians.Level); 
% 20 — eps 1e-5 — tol 10^-9
%
% Mass_20_Level.n8 = Mass_Time;
% AErr_20_Level.n8 = Err_Int;
% Sol_20_Level.n8 = Phi_t;

% Mass_20_Level.n2 = Mass_Time;
% AErr_20_Level.n2 = Err_Int;
% Sol_20_Level.n2 = Phi_t;

% Mass_20_Level.n5 = Mass_Time;
% AErr_20_Level.n5 = Err_Int;
% Sol_20_Level.n5 = Phi_t;

% Source_NL_CH_20.Times = Time_Span_ODE;
% Source_NL_CH_20.Mass  = Mass_20_Level;
% Source_NL_CH_20.Error = AErr_20_Level;
% Source_NL_CH_20.Sols  = Sol_20_Level;
% save('Source_NL_CH_20_epsB.mat', 'Source_NL_CH_20');
load('Source_NL_CH_20_epsB.mat', 'Source_NL_CH_20')

%% Order
Mass_20_Level = orderfields(Source_NL_CH_20.Mass);
AErr_20_Level = orderfields(Source_NL_CH_20.Error);
Sol_20_Level  = orderfields(Sol_20_Level);

colores = [255, 200, 87; 186, 45, 11; 86, 22, 67]/255;
%% Conservation of mass

figure(1)

% Iterate on available factors
Fields = fieldnames(Mass_20_Level);
for i = 1:numel(Fields)

    % Select vector
    Mass_Time = Mass_20_Level.(Fields{i});

    % Plot
    plot(Source_NL_CH_20.Times, Mass_Time, 'LineWidth', 1.5, 'Color', colores(i,:), ...
        'LineStyle', '-', 'Marker','x', 'DisplayName', Fields{i}(2:end))
    hold on
end

% '#e9c46a'

xlabel('$t$','Interpreter','latex'); 
ylabel( strcat(['$\int\limits_\Omega \rho(x;t) ', '\mathrm{d}x $']), 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
fontsize(16,"points")
set(gca, 'YScale', 'log')

lgd = legend('show', 'Interpreter', 'latex', 'Location','southeast');
title(lgd, 'Factors', 'Interpreter', 'latex'); % Title for the legend
ylim([1e-15,1e-12])

exportgraphics(figure(1), 'NLCH_Diffusion_Test_Mass_20.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
hold off

%% Error in t

figure(3)

% Iterate on available factors
Fields = fieldnames(Mass_20_Level);
for i = 1:numel(Fields)

    % Select vector
    Err_Int = AErr_20_Level.(Fields{i});

    % Plot
    plot(Source_NL_CH_20.Times, Err_Int, 'LineWidth', 1.5, 'Color', colores(i,:), ...
        'LineStyle', '-', 'Marker','x', 'DisplayName', Fields{i}(2:end))
    hold on
end

xlabel('$t$','Interpreter','latex'); 
ylabel( strcat(['$\big\| \rho(t) - \varphi(t) \big\|_{ L^2(\Omega) } $']), 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
fontsize(16,"points")
set(gca, 'YScale', 'log')


lgd = legend('show', 'Interpreter', 'latex', 'Location','southeast');
title(lgd, 'Factors', 'Interpreter', 'latex'); % Title for the legend

% '#2a9d8f'
ylim([1e-8,5e-5])


exportgraphics(figure(3), 'NLCH_Diffusion_Test_Error_20.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
hold off


%% 




















