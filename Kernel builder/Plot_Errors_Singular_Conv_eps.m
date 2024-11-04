function Plot_Errors_Singular_Conv_eps(Nx, Ny, de, Range_Epsilons, Abs_a, Abs_b, Abs_c, Exact, fact, store)

    colores = [186, 45, 11; 0, 105, 137; 255, 200, 87; 46, 196, 182; 86, 22, 67]/255;
    N = Nx * Ny;

    % Define a colour gradient
    colour_A = [255, 225, 168]/255;
    colour_B = [114, 61, 70]/255;
    GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
    G_colours = GRADIENT_flexible((1:de),de);

    %% First plot: swarms of each error set 
    %%%%%%%%%%%
    figure(2)
    % 10 points
    swarmchart(1.0*ones(N,1), Abs_a{1}, 5, colores(1,:), 'XJitterWidth', 0.4)
    yscale log
    hold on
    for i = 2:de
        loc = (i+1)/2;
        swarmchart(loc*ones(N,1), Abs_a{i}, 5, colores( mod(i,5) + 1,:), 'XJitterWidth', 0.4)
    end
    xlim([0.5,(de+1)/2 + 0.5]);
    
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
    ylim(MinMax .* [0.9, 1.1] );
    
    % exponents = floor( log10( abs(Range_Epsilons(1:1:end)) ) );
    % mantissas = Range_Epsilons(1:1:end) ./ (10.^exponents);
    % xtickLabelz = [];
    % for i = 1:length(mantissas)
    %     xtickLabelz{i} = sprintf('%0.0f × 10^{ %0.0f}', mantissas(i), exponents(i));
    % end
    % xticklabels(xtickLabelz)
    xticklabels([])
    
    set(gca, 'FontName', 'CMR10')
    xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
    ylabel('Absolute errors','Interpreter','latex')
    %set(gca,'xtick',[])
    set(gca, 'YTick', 10.^(-20:1:0)) %set(gca, 'YTick', logspace(-18,-1, 18))

    if store
        file_name = strcat('Test_Singular_Convolution_Eps_AbsE_Swarm[', num2str(Nx), ',', num2str(fact*Nx), '].pdf');
        exportgraphics(figure(2), file_name, 'BackgroundColor','none','ContentType','vector')
    end
    hold off
   
    %% Second plot: error curves per node
    Rel_A = zeros([de,N]);
    for i = 1:de
        Rel_A(i,:) = Abs_a{i};
    end
    [Mab, Iab] = sort(Rel_A(:,1), 'ascend');
    G_colours = G_colours(Iab,:);

    figure(3)
    plot(Range_Epsilons,Rel_A);     hold on
    plot(Range_Epsilons,mean(Rel_A, 2) + std(Rel_A, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','-.')
    plot(Range_Epsilons,mean(Rel_A, 2), 'Color', [39, 39, 39]/255,'LineWidth',2,'LineStyle',':')
    plot(Range_Epsilons,median(Rel_A, 2), 'Color', [195, 47, 39]/255,'LineWidth',2,'LineStyle',':')
    plot(Range_Epsilons,mean(Rel_A, 2) - std(Rel_A, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','--')
    yscale log
    xscale log
    set(gca, 'YTick', 10.^(-20:1:0))
    colororder(G_colours)

    ylim(MinMax .* [0.9, 1.1] );
    xlim([Range_Epsilons(1), Range_Epsilons(end)]);
    set(gca, 'FontName', 'CMR10')
    xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
    ylabel('Absolute errors','Interpreter','latex')

    if store
        file_name = strcat('Test_Singular_Convolution_Eps_AbsE_Cont[', num2str(Nx), ',', num2str(fact*Nx), '].pdf');
        exportgraphics(figure(3), file_name, 'BackgroundColor','none','ContentType','vector')
    end
    hold off

    %% Third plot: Test the whole integral approximation by subtracting the piece that is missing
    Exact_v = zeros([1,de]);
    for i = 1:de
        Exact_v(i) = max(abs(Exact{i}));
    end
    Abs_B = zeros([de,N]);
    for i = 1:de
        Abs_B(i,:) = Abs_b{i};
    end

    figure(4)
    plot(Range_Epsilons,Abs_B);    hold on
    plot(Range_Epsilons,mean(Abs_B, 2) + std(Abs_B, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','-.')
    plot(Range_Epsilons,mean(Abs_B, 2), 'Color', [39, 39, 39]/255,'LineWidth',2,'LineStyle',':')
    plot(Range_Epsilons,median(Abs_B, 2), 'Color', [195, 47, 39]/255,'LineWidth',2,'LineStyle',':')
    plot(Range_Epsilons,mean(Abs_B, 2) - std(Abs_B, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','--')
    plot(Range_Epsilons, Exact_v, 'black','LineWidth',2,'LineStyle','--')
    xlim([Range_Epsilons(1), Range_Epsilons(end)]);
    [Min_b, Max_a] = bounds(Abs_B(Abs_B>0));
    ylim([Min_b, Max_a] .* [0.9, 1.1] );
    xscale log;    yscale log
    set(gca, 'YTick', 10.^(-20:1:0))
    colororder(G_colours)
    
    % Identify largest curve (prints std):
    [Mab, Iab] = max(Abs_B(1:5,:),[],2);   std(Mab)
    plot(Range_Epsilons(:), Abs_B(:,Iab(1) ), 'LineWidth',2);
    set(gca, 'FontName', 'CMR10')
    xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
    ylabel('Absolute errors','Interpreter','latex')
    
    if store
        file_name = strcat('Test_Singular_Convolution_Eps_Domain_Cont[', num2str(Nx), ',', num2str(fact*Nx), '].pdf');
        exportgraphics(figure(4), file_name, 'BackgroundColor','none','ContentType','vector')
    end
    hold off
    
    %% Fourth plot: Plot error against full convolution (adding exact part)

    % Abs_C = zeros([de,N]);
    % for i = 1:de
    %     Abs_C(i,:) = Abs_c{i};
    % end
    % 
    % figure(5)
    % plot(Range_Epsilons,Abs_C);    hold on
    % plot(Range_Epsilons,mean(Abs_C, 2) + std(Abs_C, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','-.')
    % plot(Range_Epsilons,mean(Abs_C, 2), 'Color', [39, 39, 39]/255,'LineWidth',2,'LineStyle',':')
    % plot(Range_Epsilons,median(Abs_C, 2), 'Color', [195, 47, 39]/255,'LineWidth',2,'LineStyle',':')
    % plot(Range_Epsilons,mean(Abs_C, 2) - std(Abs_C, [],2), 'Color', [188, 189, 192]/255,'LineWidth',1,'LineStyle','--')
    % xlim([Range_Epsilons(1), Range_Epsilons(end)]);
    % [Min_b, Max_a] = bounds(Abs_C(Abs_C>0));
    % ylim([Min_b, Max_a] .* [0.9, 1.1] );
    % xscale log;    yscale log
    % set(gca, 'YTick', 10.^(-20:2:-1))
    % colororder(G_colours)
    % 
    % set(gca, 'FontName', 'CMR10')
    % xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
    % ylabel('Absolute errors','Interpreter','latex')
    % 
    % if store
    %     file_name = strcat('Test_Singular_Convolution_Eps_Complete_Cont[', num2str(Nx), ',', num2str(fact*Nx), '].pdf');
    %     exportgraphics(figure(5), file_name, 'BackgroundColor','none','ContentType','vector')
    % end
    % hold off

end