figure(500)


% Setup
Tile = tiledlayout(1,5,'TileSpacing','compact');

% Axes
for i = 1:5
    h_S(i) = nexttile;
end

% Initial plot
axes(h_S(1))
aBox.plot(phi_ic,{});
xlabel('$x_1$'); ylabel('$x_2$'); caxis([-1,1])
title('$t=0$','Interpreter','latex');    shading interp

set(h_S(1),'fontsize',6,'linewidth',0.1);

h_S(1).CameraPosition = [0.5 0.5 30];

set(h_S(1), 'TickLabelInterpreter', 'latex');

% Interpolate solution at specified times
%TimesEval = [0,0.01,0.05,0.07,1];      % Test A
%TimesEval = [0, 0.01, 0.04, 0.05,1];      % Test B
TimesEval = [0,0.001, 0.0025, 0.005, 1];      % Test C

for i = 2:5
    t = TimesEval(i);
    IP = TimeLine.ComputeInterpolationMatrixPhys(t).InterPol;
    rho_t = (IP * Phi_to)';
    
    % Add panels
    axes(h_S(i))
    aBox.plot(rho_t,{});
    xlabel('$x_1$'); ylabel('');
    
    zlim([-1-4e-1,1+1e-1]); caxis([-1,1]);    shading interp
    title(['$t$ = ' num2str(TimesEval(i))],'Interpreter','latex')

    set(h_S(i),'fontsize',6,'linewidth',0.1,'yticklabel',[]);

    h_S(i).CameraPosition = [0.5 0.5 30];
    set(h_S(i), 'TickLabelInterpreter', 'latex');
end

% Add additional comments
set(h_S(5),'YAxisLocation', 'right')

axes(h_S(5));    ylabel('$\rho(x;t)$');

colormap bone

fontsize(11,"points")

exportgraphics(figure(500), strcat('A_Solution_PS.pdf'), 'BackgroundColor','none','ContentType','vector', 'Resolution', 300)

%colormap default

%cb = colorbar;
%cb.Layout.Tile = 'east';
%cb.Position = [0.92,0.12,0.02,0.8];


