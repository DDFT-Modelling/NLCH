function Plot_Angles_Domain(angles, Pts, store)
    % Define a colour gradient
    colour_A = [255, 225, 168]/255;
    colour_B = [114, 61, 70]/255;
    GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);

    % Convert angles to degrees for better precision and easier handling
    degrees = rad2deg(angles);      % I am really curious as to why this function exists
    
    [MinA, MaxA]  = bounds( degrees(degrees > 0) );    % We can expect a min angle of 45° and a max angle close to 90°
    % Set number of levels in which we will split the colour map for better
    % understanding of the angles: 45-90
    de = 5;%length(45:5:90);
    Lower_Bounds = [45, 50, 70, 80, 85];      %45:5:85;
    Upper_Bounds = [50, 70, 80, 85, 90];      %50:5:90;
    
    % Compute colour map
    %G_colours = GRADIENT_flexible((1:de),de);
    G_colours = [ colour_A; GRADIENT_flexible(3,de); colour_B; [255, 138, 91]/255; [0,0,0] ];
    % For each positive value of error, pick a colour
    scatter_Cmap = ones(length(angles), 3);
    mode_index   = 5;
    scatter_Indx = zeros(length(angles), 1);
    for i = 1:length(angles)
        if degrees(i) > 0
            % Create an index in G_colours:
            ind = find( degrees(i) >= Lower_Bounds & degrees(i) < Upper_Bounds,1);
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
        merr =  max(degrees(scatter_Indx == Us_Indx(i)));
        exponent = floor( log10( merr ) );
        mantissa = merr / (10^exponent);
        %StrRep   = sprintf('%0.0f  10^{ %0.0f}', mantissa, exponent);
        %MeanErrs{i} = strcat( ['$\sim ', StrRep, '$'] );
        MeanErrs{i} = strcat( ['$\leq ', sprintf('%0.0f', merr),'^{\circ}$'] );
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
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1, 1, 1], 'DisplayName', '---');
        elseif scatter_Indx(i) == mode_index
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [234, 82, 111]/255, 'DisplayName', MeanErrs{i});
        else
            dummyP(i) = scatter(nan, nan, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', G_colours(scatter_Indx(i),:), 'DisplayName', MeanErrs{i} );
        end
    end
    
    lgd = legend(dummyP,'Location','northeast', 'interpreter', 'latex');
    set(gca, 'FontName', 'CMR10')

    if store
        file_name = strcat('Test_Angles[', num2str(Pts.N1), '].pdf');
        %exportgraphics(figure(4), file_name, 'BackgroundColor','none','ContentType','vector')
        exportgraphics(figure(4), file_name, 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300)
    end

    hold off

end