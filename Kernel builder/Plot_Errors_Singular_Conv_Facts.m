function Plot_Errors_Singular_Conv_Facts(Nx, Ny, de, Range_Epsilons, Abs_a, Abs_b, Abs_c, Exact, Factors, store)

    colores = [255, 200, 87; 186, 45, 11; 86, 22, 67]/255;
    %[186, 45, 11; 0, 105, 137; 255, 200, 87]/255; %; 46, 196, 182; 86, 22, 67]/255;
    N = Nx * Ny;

    % Define a colour gradient
    colour_A = [255, 225, 168]/255;
    colour_B = [114, 61, 70]/255;
    GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
    G_colours = GRADIENT_flexible((1:(de/2)),de/2);

    %% First plot: swarms of each error set 
    %%%%%%%%%%%
    figure(2)
    a(1) = plot(nan, nan, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'black', 'DisplayName', 'Median'); hold on
    a(2) = plot(nan, nan, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'black', 'DisplayName', 'Mode');
    a(3) = plot(nan, nan, 'LineWidth', 1, 'LineStyle', '--', 'Color', 'black', 'DisplayName', 'Mean');
    % 10 points
    swarmchart(1.0*ones(N,1), Abs_a{1}, 10, colores(2,:), 'XJitterWidth', 0.3)
    line([1-0.15,1+0.15],[median(Abs_a{1}),median(Abs_a{1})], 'LineWidth', 2, 'LineStyle', '-', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
    line([1-0.15,1+0.15],[mode(Abs_a{1}),mode(Abs_a{1})], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
    line([1-0.15,1+0.15],[mean(Abs_a{1}),mean(Abs_a{1})], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
    yscale log
    %hold on
    for i = 2:de
        loc = (i+1)/2;
        swarmchart(loc*ones(N,1), Abs_a{i}, 10, colores( mod(i,3) + 1,:), 'XJitterWidth', 0.3)
        line([loc-0.15,loc+0.15],[median(Abs_a{i}),median(Abs_a{i})], ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
        line([loc-0.15,loc+0.15],[mode(Abs_a{i}),mode(Abs_a{i})], ...
            'LineWidth', 2, 'LineStyle', ':', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
        line([loc-0.15,loc+0.15],[mean(Abs_a{i}),mean(Abs_a{i})], ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
    end
    xlim([0.75,(de+1)/2 + 0.25]);
    
    MinMax = [1e+10,0];
    for i = 1:de
        [mina, maxa] = bounds(Abs_a{i}(Abs_a{i}>0),'all');
        if mina < MinMax(1)
            MinMax(1) = mina;
        end
        if maxa > MinMax(2)
            MinMax(2) = maxa;
        end
    end
    %ylim([3e-18, 2e-5]);
    ylim(MinMax .* [0.65, 1.5] );
    % Add vertical lines as separators between classes
    xline( 0.5*(2 + 2.5) ,'-.');


    % Additional decoration
    lgd = legend(a, 'Interpreter','latex','Location','southeast');
    title(lgd,'Descriptors','Interpreter','latex')
    % 
    % 
    % 
    xticks([1, 1.5, 2.0, 2.5, 3.0, 3.5 ]) %, 4.0, 4.5])
    % 
    divisions = {};
    for i = 1:de
        divisions{i} = strcat( num2str(Nx), num2str(Factors( mod(i-1,de/2) + 1 ),'/%d'));
    end
    
    %xlabelArray = [ divisions; {'10^{-2}','10^{-2}', '10^{-5}', '10^{-5}'} ];  
    xlabelArray = [ divisions; {' ','10^{-2}',' ',  ' ',' 10^{-5}',' '} ]; 
    xtickLabels = strtrim(sprintf('%s\\newline%s\n', xlabelArray{:}));
    % 
    % %xticklabels({'A','B','C','D','E',  'A','B','C','D','E',  'A','B','C','D','E',  'A','B','C','D','E'})
    % %xticklabels({10, 10, 10, 10, 10,  20, 20, 20, 20, 20,  30, 30, 30, 30, 30,  40, 40, 40, 40, 40})
    xticklabels(xtickLabels)
    set(gca, 'FontName', 'CMR10')
    xlabel('Grid resolution and detail factor','Interpreter','latex')
    ylabel('Absolute errors','Interpreter','latex')
    set(gca, 'YTick', 10.^(-20:1:-1)) %set(gca, 'YTick', logspace(-18,-1, 18))

    if store
        file_name = strcat('Test_Singular_Convolution_Fixed_Swarm[', num2str(Nx), ',', num2str(Factors,'%d,'), num2str(de/2), '].pdf');
        exportgraphics(figure(2), file_name, 'BackgroundColor','none','ContentType','vector')
    end
    hold off

end