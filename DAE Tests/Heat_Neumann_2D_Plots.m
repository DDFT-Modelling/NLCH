for k = [20,30,40]
    [Phi_t, u_t, Mass_Error, outTimes] = Heat_Neumann_2D(k);
    Errors{k} = Mass_Error;
    Times{k} = outTimes;
end

%% Plot errors

plot(Times{20}, abs(Errors{20}), 'LineWidth', 1.5, 'Color', '#2a9d8f', 'LineStyle', '-', 'Marker','x', 'DisplayName', '20','MarkerSize',9)
hold on
plot(Times{30}, abs(Errors{30}), 'LineWidth', 1.5, 'Color', '#e9c46a', 'LineStyle', '-', 'Marker','o', 'DisplayName', '30','MarkerSize',9)
plot(Times{40}, abs(Errors{40}), 'LineWidth', 1.5, 'Color', '#C84630', 'LineStyle', '-', 'Marker','+', 'DisplayName', '40','MarkerSize',9)

xscale log
yscale log

%ylabel( strcat(['$1-\int\limits_\Omega u(x;t)/c_1(t) ', '\mathrm{d}x $']), 'Interpreter','latex');
ylabel( strcat(['$\bigg| 1-\int\limits_\Omega u(x;t)/c_1(t) ', '\mathrm{d}x \bigg|$']), 'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('$t$','Interpreter','latex'); 
fontsize(14,"points")

lgd = legend('show', 'Interpreter', 'latex', 'Location','northwest');
title(lgd, '$N$', 'Interpreter', 'latex'); % Title for the legend

plot([1e-4, Times{20}(2)], abs([Errors{20}(1), Errors{20}(2)]),'LineWidth', 1.5, 'Color', '#2a9d8f', 'HandleVisibility', 'off')
plot([1e-4, Times{30}(2)], abs([Errors{30}(1), Errors{30}(2)]),'LineWidth', 1.5, 'Color', '#e9c46a', 'HandleVisibility', 'off')
plot([1e-4, Times{40}(2)], abs([Errors{40}(1), Errors{40}(2)]),'LineWidth', 1.5, 'Color', '#C84630', 'HandleVisibility', 'off')
%plot([1e-3, Times{100}(2)], [Errors{100}(1), Errors{100}(2)],'LineWidth', 1.5, 'Color', '#331832', 'HandleVisibility', 'off')

% Store
exportgraphics(figure(1), 'HN_2D_Solution_Error.pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)