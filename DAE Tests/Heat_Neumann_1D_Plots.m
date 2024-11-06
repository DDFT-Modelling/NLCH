%% Plot exact solution
% figure(1)
% u_t = zeros(n_t,n_x);
% c_0 = aLine.Int * (y.^2);
% u_t = exp(-[outTimes; (8.1:0.1:10)']) * (y.^2)' / c_0;
% 
% contourf(y,[outTimes; (8.1:0.1:10)'],u_t)
% colormap bone
% xlabel('$x$', 'Interpreter', 'latex')
% ylabel('$t$', 'Interpreter', 'latex')
% yscale log
% 
% colorbar('Ticks',[min(u_t,[],'all')*1.08, 0.13, 0.26, 0.39, 0.5], 'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'}, ...
%     'TickLabelInterpreter', 'latex');
% 
% fontsize(14,"points")
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'FontName', 'CMR10')
% 
% exportgraphics(figure(1), 'HN_1D_Solution_Contour.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)

%% Iterate errors

for k = [20,40,60,80,100]
    [Phi_t, u_t, Mass_Error, outTimes] = Heat_Neumann_1D(k);
    Errors{k} = Mass_Error;
    Times{k} = outTimes;
end

%% Plot errors

plot(Times{40}, Errors{40}, 'LineWidth', 1.5, 'Color', '#2a9d8f', 'LineStyle', '-', 'Marker','x', 'DisplayName', '40','MarkerSize',9)
hold on
plot(Times{60}, Errors{60}, 'LineWidth', 1.5, 'Color', '#e9c46a', 'LineStyle', '-', 'Marker','o', 'DisplayName', '60','MarkerSize',9)
plot(Times{80}, Errors{80}, 'LineWidth', 1.5, 'Color', '#C84630', 'LineStyle', '-', 'Marker','+', 'DisplayName', '80','MarkerSize',9)
plot(Times{100}, Errors{100}, 'LineWidth', 1.5, 'Color', '#331832', 'LineStyle', '-', 'Marker','*', 'DisplayName', '100','MarkerSize',9)

xscale log

%ylabel( strcat(['$1-\int\limits_\Omega u(x;t)/c_1(t) ', '\mathrm{d}x $']), 'Interpreter','latex');
ylabel( strcat(['$\bigg(c_1(t)-\int\limits_\Omega u(x;t) ', '\mathrm{d}x \bigg)/ \max \big\{c_1(t), 5\times 10^{-3} \big\} $']), 'Interpreter','latex');
%ylabel( strcat(['$c_1(t)-\int\limits_\Omega u(x;t) ', '\mathrm{d}x $']), 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('$t$','Interpreter','latex'); 
fontsize(14,"points")

lgd = legend('show', 'Interpreter', 'latex', 'Location','southwest');
title(lgd, '$N$', 'Interpreter', 'latex'); % Title for the legend

plot([1e-3, Times{40}(2)], [Errors{40}(1), Errors{40}(2)],'LineWidth', 1.5, 'Color', '#2a9d8f', 'HandleVisibility', 'off')
plot([1e-3, Times{60}(2)], [Errors{60}(1), Errors{60}(2)],'LineWidth', 1.5, 'Color', '#e9c46a', 'HandleVisibility', 'off')
plot([1e-3, Times{80}(2)], [Errors{80}(1), Errors{80}(2)],'LineWidth', 1.5, 'Color', '#C84630', 'HandleVisibility', 'off')
plot([1e-3, Times{100}(2)], [Errors{100}(1), Errors{100}(2)],'LineWidth', 1.5, 'Color', '#331832', 'HandleVisibility', 'off')

% Store
exportgraphics(figure(1), 'HN_1D_Solution_Error.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)