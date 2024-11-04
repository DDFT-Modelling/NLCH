function Plot_Errors_Singular_Conv_Domain(epsilon,errors, Pts, store)
    % Define a colour gradient
    colour_A = [255, 225, 168]/255;
    colour_B = [114, 61, 70]/255;
    GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
    
    [MinA, MaxA]  = bounds( -log( errors(errors > 0) ));
    % Compute number of levels in which we will split the colour map for better
    % understanding of the errors
    de = ceil(MaxA) - floor(MinA) + 1;
    
    % Compute colour map
    G_colours = GRADIENT_flexible((1:de),de);
    % For each positive value of error, pick a colour
    scatter_Cmap = zeros(length(errors), 3);
    mode_index   = ceil(-log(mode(errors)) - floor(MinA));
    scatter_Indx = zeros(length(errors), 1);
    for i = 1:length(errors)
        if errors(i) > 0
            % Create an index in G_colours:
            ind = ceil(-log(errors(i)) - floor(MinA));
            if (mode_index - 0) == ind
                scatter_Cmap(i,:) = [234, 82, 111]/255;
            else
                scatter_Cmap(i,:) = G_colours(ind,:);
            end
            scatter_Indx(i) = ind;
        end
    end

    % Average errors for labels
    Us_Indx  = unique(scatter_Indx);
    MeanErrs = cell(length(Us_Indx),1);
    for i = 1:length(Us_Indx)
        merr =  mean(errors(scatter_Indx == Us_Indx(i)));
        exponent = floor( log10( merr ) );
        mantissa = merr / (10^exponent);
        %StrRep   = sprintf('%0.0f  10^{ %0.0f}', mantissa, exponent);
        %MeanErrs{i} = strcat( ['$\sim ', StrRep, '$'] );
        MeanErrs{i} = strcat( ['$\sim ', sprintf('%0.1f', mantissa), '\times', sprintf('10^{ %0.0f}', exponent), '$'] );
    end
    %MeanErrs
    
    % PLOT
    figure(4)
    s = scatter(Pts.y1_kv, Pts.y2_kv, 'filled');
    s.CData = scatter_Cmap;
    hold on
    
    % Dummy plots for legend
    scatter_Indx = Us_Indx;
    for i = 1:length(scatter_Indx)
        if scatter_Indx(i) == 0
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0, 0, 0], 'DisplayName', '$0$');
        elseif scatter_Indx(i) == mode_index
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [234, 82, 111]/255, 'DisplayName', MeanErrs{i});
        else
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', G_colours(scatter_Indx(i),:), 'DisplayName', MeanErrs{i} );
        end
    end
    
    % Lower
    Line(1) = plot([2*epsilon 1-2*epsilon], [2*epsilon 2*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    Line(2) = plot([2*epsilon 1-2*epsilon], [1-2*epsilon 1-2*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    % Upper
    Line(3) = plot([3*epsilon 1-3*epsilon], [3*epsilon 3*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    Line(4) = plot([3*epsilon 1-3*epsilon], [1-3*epsilon 1-3*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    % Left
    Line(5) = plot([2*epsilon 2*epsilon], [2*epsilon 1-2*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    Line(6) = plot([3*epsilon 3*epsilon], [3*epsilon 1-3*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    % Right
    Line(7) = plot([1-2*epsilon 1-2*epsilon], [2*epsilon 1-2*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    Line(8) = plot([1-3*epsilon 1-3*epsilon], [3*epsilon 1-3*epsilon], 'Color', [1,1,1]/1.5, 'LineWidth', 1);
    
    uistack(Line, 'bottom');
    
    lgd = legend(dummyP,'Location','northeast', 'interpreter', 'latex');
    set(gca, 'FontName', 'CMR10')

    if store
        file_name = strcat('Test_Singular_Convolution_Domain[', num2str(length(errors)), ',', num2str(epsilon), '].pdf');
        %exportgraphics(figure(4), file_name, 'BackgroundColor','none','ContentType','vector')
        exportgraphics(figure(4), file_name, 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
    end

    hold off

end