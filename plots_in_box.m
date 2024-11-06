function plots_in_box(Rho_t, aBox,outTimes,filename, title_plot)

arguments
    Rho_t       double;
    aBox        Box;
    outTimes    double;
    filename    char = strcat( 'NLCH_Diffusion_Test.gif' );
    title_plot  char = '$\varphi(x,t)$';
end

    % Plot solution
    h  = figure('Position',[100,100,400,300]);
    axis tight manual % this ensures that getframe() returns a consistent size

    MaxR = max(Rho_t,[],'all');   MinR = min(Rho_t,[],'all');
    % Add a little extra of height with the same magnitude
    MaxR = MaxR + 10^floor( log10(abs(MaxR)) );
    
    hR = subplot(1,1,1);  zlim(hR,[MinR MaxR]);  %caxis(hR,[0 MaxR + abs(MinR) ]);

    set(gca,'Color','none');   set(gcf, 'color', 'white');
    set(gca, 'TickLabelInterpreter', 'latex');      % Change tick labels to LaTeX

    opts = {};  % plotting options - the default are ok for now
    

    for iTime = 1:length(outTimes)

        % extract the solution at this time
        Rt = Rho_t(iTime,:)';

        % plot with built-in class function
        hold off

        % Plot Rho
        axes(hR)
        aBox.plot(Rt,opts);
        zlim(hR,[MinR MaxR]);   %caxis(hR,[0 MaxR + abs(MinR) ]);
        xlabel('$x_1$'); ylabel('$x_2$'); 
        zlabel(title_plot,'Interpreter','latex')
        set(gca, 'TickLabelInterpreter', 'latex');
        colormap('bone')


        % set title with some info
        sgtitle(['$t$ = ' num2str(outTimes(iTime))],'Interpreter','latex','FontSize',20)

        % fix axes
    %         colorbar

        %pause(0.05)

        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File 
        if iTime == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
        end 

    end
    hold off
end