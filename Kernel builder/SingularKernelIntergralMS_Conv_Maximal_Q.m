function [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(epsilon, Nx, Ny, fact)
% Test for convolution weight with maximal number of multishapes around an
% excluded region around each point. 
% We create the following number of MS for each case:
%       Corner:     3 MS
%       Border:     5 MS
%       Inner:      8 MS
    
    arguments
        epsilon double = 1e-3
        Nx      double = 10
        Ny      double = 10
        fact    double = 1
    end

    %----------------------------------------------------------------------
    % Global domain options
    %----------------------------------------------------------------------

    % number of points
    %Nx = 10;
    %Ny = 10;
    
    %----------------------------------------------------------------------
    % Size of your cut-out region
    %----------------------------------------------------------------------
    
    %epsilon = 0.001;

    %----------------------------------------------------------------------
    % Factor of sub domain division
    %----------------------------------------------------------------------

    % factor = 1;
    
    %----------------------------------------------------------------------
    % Set up standard box to do the full computation on
    %----------------------------------------------------------------------
    
    geom.N = [Nx;Ny];
    geom.y1Min = 0;
    geom.y1Max = 1;
    geom.y2Min = 0;
    geom.y2Max = 1;

    % Guarantee that the value of epsilon is always less than minimum
    % length of box
    geom.lengths = [geom.y1Max - geom.y1Min, geom.y2Max - geom.y2Min];
    if any([epsilon, epsilon] >= geom.lengths)
        epsilon = 0.1 * min(geom.lengths);
    end
    
    % Make a copy for plotting
    %geomPlot = geom;
    %geomPlot.N = [40;40];
    
    % Create an abstract box, we will move this along the domain
    aBox = Box(geom);
    %aBox.ComputeAll(geomPlot);   % not needed as we are not computing any
    %convolution or diff
    % aBox.PlotGrid;

    % Add a convolution kernel
    Kernel = @(y1,y2) 0.25 * log(y1.^2 + y2.^2)/pi;    % One operation less to worry about
    
    y1_kv = aBox.Pts.y1_kv;
    y2_kv = aBox.Pts.y2_kv;
    Points = [y1_kv, y2_kv];
    Newtonians = Newtonian(Points, epsilon);
    
    % Convolve a function
    test_1 = @(y1,y2) ones(size(y1));

    %Vals = Kernel(y1_kv, y2_kv);
    %Vals(isinf(Vals) | isnan(Vals)) = -2;    % A small value is required
    %aBox.plot(Vals,{})

    
    % Evaluate function
    fBox = test_1(y1_kv, y2_kv);
    % Store values of convolution outside the box
    t_Full = Newtonians.NI;
    all_sq = Newtonians.NJ;
    smaller = Newtonians.NG;
    
    %----------------------------------------------------------------------
    % Set up epsilon box (centred at zero)
    %----------------------------------------------------------------------
    
    geomEps.N = [Nx;Ny];
    geomEps.y1Min = -epsilon;
    geomEps.y1Max =  epsilon;
    geomEps.y2Min = -epsilon;
    geomEps.y2Max =  epsilon;
    geomEps.shape = 'Box';
    
    %----------------------------------------------------------------------
    % Find the intersection with the epsBox centred at one point
    %----------------------------------------------------------------------

    Corners = [ [aBox.y1Min, aBox.y2Min], % lower left (1)
                [aBox.y1Min, aBox.y2Max], % upper left
                [aBox.y1Max, aBox.y2Min], % lower right
                [aBox.y1Max, aBox.y2Max] ]; % upper right (4)

    errors  = zeros(Nx*Ny,1);     % We are going to store the errors
    experiment = zeros(Nx*Ny,1);  % How well are we approximating the full integral without adding the bits?
    full_thing = zeros(Nx*Ny,1);  % Same but adding the bits
    Conv = zeros(Nx*Ny,Nx*Ny);    % Store convolution matrix as well
    
    %----------------------------------------------------------------------
    % - loop over all of the points in aBox
    % - Compute the intersection with the shifted epsBox, which gives the
    % corners
    % - Construct MS with these corners (depends on which edge/corner
    % you're at)
    % - Use the MS to construct the interpolation
    % - Construct a convolution matrix on aBox
    %----------------------------------------------------------------------
    [nx, ny] = deal( round(fact*Nx), round(fact*Ny) );        % # of divisions in capBox % Nx/2, Ny/2 [!!!!!!]
    %%
    for i = 1:(Nx * Ny)
        % Define point
        [y1, y2] = deal(y1_kv(i), y2_kv(i));
        % centre box at [y1;y2];
        geomEpsPt = ShiftGeometryOrigin(geomEps,[y1;y2]);
        epsBox = Box(geomEpsPt);
        
        % Visualise box centred at (y_1,y_2) across grid
        %hold on
        %epsBox.PlotGrid;
        %pause(0.5)

        % Compute intersection of main box and the small box
        capBox = Intersect_Box_Box(epsBox, aBox);
        %hold on
        %capBox.PlotGrid;
        %pause(0.5)

        % --------------------------------------------------------------
        % --------------------------------------------------------------
        % Set up of MS
        % --------------------------------------------------------------
        % --------------------------------------------------------------
        % Extract general corners
        [left, right] = deal( capBox.y1Min, capBox.y1Max );
        [bottom, top] = deal( capBox.y2Min, capBox.y2Max);
        capLines = [left, right, bottom, top];
        
        %---------------------------------------------------------------
        % Set up MultiShape domain                                            
        %---------------------------------------------------------------
        % Detect which kind of shape we will build
        Checker = ismembertol(capLines, [Corners(1,:), Corners(4,:)]);
        % If any of the logical tests is true, then we are at either a
        % corner or a border
        if any(Checker)
            % Now we have to determine if we are at a corner or border
            Border_Line_Position = find(Checker);
            % By construction, epsilon < geom.lengths so at most two points
            % can be detected simultaneously (corner)
            if sum(Checker) == 1
                % A border is identified: 5 MS
                % Position tells us which border we are dealing with
                % 1. left, 2. right, 3. bottom, 4. top
                switch Border_Line_Position
                    case 1
                        % Lower left
                        q_geom_a.Y = [ 0 0; 0 bottom; right bottom; right 0];
                        q_geom_a.N = [nx; ny];
                        % Lower right
                        q_geom_b.Y = [ right 0; right bottom; 1 bottom; 1 0];
                        q_geom_b.N = [nx; ny];
                        % Centre
                        q_geom_c.Y = [ right bottom; right top; 1 top; 1 bottom];
                        q_geom_c.N = [nx; ny];
                        % Upper right
                        q_geom_e.Y = [ right top; right 1; 1 1; 1 top];
                        q_geom_e.N = [nx; ny];
                        % Upper right
                        q_geom_f.Y = [ 0 top; 0 1; right 1; right top];
                        q_geom_f.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b,q_geom_c,q_geom_e,q_geom_f];
                    case 2
                        % Lower right
                        q_geom_a.Y = [ left 0; left bottom; 1 bottom; 1 0];
                        q_geom_a.N = [nx; ny];
                        % Lower left
                        q_geom_b.Y = [ 0 0; 0 bottom; left bottom; left 0];
                        q_geom_b.N = [nx; ny];
                        % Centre
                        q_geom_c.Y = [ 0 bottom; 0 top; left top; left bottom];
                        q_geom_c.N = [nx; ny];
                        % Upper left
                        q_geom_e.Y = [ 0 top; 0 1; left 1; left top];
                        q_geom_e.N = [nx; ny];
                        % Upper right
                        q_geom_f.Y = [ left top; left 1; 1 1; 1 top];
                        q_geom_f.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b,q_geom_c,q_geom_e,q_geom_f];
                    case 3
                        % Lower left
                        q_geom_a.Y = [ 0 0; left 0; left top; 0 top];
                        q_geom_a.N = [nx; ny];
                        % Upper left
                        q_geom_b.Y = [ 0 top; 0 1; left 1; left top];
                        q_geom_b.N = [nx; ny];
                        % Centre
                        q_geom_c.Y = [ left top; left 1; right 1; right top];
                        q_geom_c.N = [nx; ny];
                        % Upper right
                        q_geom_e.Y = [ right top; right 1; 1 1; 1 top];
                        q_geom_e.N = [nx; ny];
                        % Lower right
                        q_geom_f.Y = [ right 0; right top; 1 top; 1 0];
                        q_geom_f.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b,q_geom_c,q_geom_e,q_geom_f];
                    case 4
                        % Upper left
                        q_geom_a.Y = [ 0 bottom; 0 1; left 1; left bottom];
                        q_geom_a.N = [nx; ny];
                        % Lower left
                        q_geom_b.Y = [ 0 0; 0 bottom; left bottom; left 0];
                        q_geom_b.N = [nx; ny];
                        % Centre
                        q_geom_c.Y = [ left 0; left bottom; right bottom; right 0];
                        q_geom_c.N = [nx; ny];
                        % Lower right
                        q_geom_e.Y = [ right 0; right bottom; 1 bottom; 1 0];
                        q_geom_e.N = [nx; ny];
                        % Upper right
                        q_geom_f.Y = [ right bottom; right 1; 1 1; 1 bottom];
                        q_geom_f.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b,q_geom_c,q_geom_e,q_geom_f];
                end

            else
                % A corner has been identified: Can degenerate to a box split if too large. 
                if sum(Checker) > 2
                    % Degenerate corner identified
                    
                    % Bottom box:
                    if all(Checker == [1 1 1 0]) 
                        % Top
                        q_geom_b.Y = [ 0 top; 0 1; 1 1; 1 top];
                        q_geom_b.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_b];
                    % Top box
                    elseif all(Checker == [1 1 0 1])
                        % Lower
                        q_geom_a.Y = [ 0 0; 0 top; 1 top; 1 0];
                        q_geom_a.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a];
                    % Left box
                    elseif all(Checker == [0 1 1 1])
                        % Right
                        q_geom_a.Y = [ 0 left; 0 1; left 1; left 0];
                        q_geom_a.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a];
                    % Right box
                    elseif all(Checker == [1 0 1 1])
                        % Left
                        q_geom_a.Y = [ right 0; right 1; 1 1; 1 0];
                        q_geom_a.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a];
                    else
                        error('Radius too large! Cannot create empty box.')
                    end
                        
                    
                elseif sum(Checker) == 2
                    % Now we can continue: Corners have 3 MS
                    switch Border_Line_Position * [2,1]' - 4
                    case 1 % [1 3] bottom left
                        % Right
                        q_geom_a.Y = [ right 0; right top; 1 top; 1 0];
                        q_geom_a.N = [nx; ny];
                        % Upper
                        q_geom_b.Y = [ 0 top; 0 1; right 1; right top];
                        q_geom_b.N = [nx; ny];
                        % Diagonal
                        q_geom_c.Y = [ right top; right 1; 1 1; 1 top];
                        q_geom_c.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b, q_geom_c];

                    case 2 % [1 4] top left
                        % Right
                        q_geom_a.Y = [ right bottom; right 1; 1 1; 1 bottom];
                        q_geom_a.N = [nx; ny];
                        % Lower
                        q_geom_b.Y = [ 0 0; left bottom; right bottom; right 0];
                        q_geom_b.N = [nx; ny];
                        % Diagonal
                        q_geom_c.Y = [ right 0; right bottom; 1 bottom; 1 0];
                        q_geom_c.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b, q_geom_c];

                    case 3 % [2 3] bottom right
                        % Left
                        q_geom_a.Y = [ 0 0; 0 top; left top; left bottom];
                        q_geom_a.N = [nx; ny];
                        % Upper
                        q_geom_b.Y = [ left top; left 1; 1 1; right top];
                        q_geom_b.N = [nx; ny];
                        % Diagonal
                        q_geom_c.Y = [ 0 top; 0 1; left 1; left top];
                        q_geom_c.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b, q_geom_c];
                    case 4 % [2 4] top right
                        % Left
                        q_geom_a.Y = [ 0 bottom; 0 1; left 1; left bottom];
                        q_geom_a.N = [nx; ny];
                        % Lower
                        q_geom_b.Y = [ left 0; left bottom; 1 bottom; 1 0];
                        q_geom_b.N = [nx; ny];
                        % Diagonal
                        q_geom_c.Y = [ 0 0; 0 bottom; left bottom; left 0];
                        q_geom_c.N = [nx; ny];
                        % Collect all
                        q_geom = [q_geom_a,q_geom_b, q_geom_c];
                end
                end
            end
            %
        else
            % In this case, we divide the shape into eight areas
            % Lower left
            q_geom_a.Y = [0 0; 0 bottom; left bottom; left 0]; %[ 0 0; 1 0; right bottom; left bottom];
            q_geom_a.N = [nx; ny];
            % Lower centre
            q_geom_b.Y = [ left 0; left bottom; right bottom; right, 0];
            q_geom_b.N = [nx; ny];
            % Lower right
            q_geom_c.Y = [ right 0; right bottom; 1 bottom; 1, 0];
            q_geom_c.N = [nx; ny];

            % Centre left
            q_geom_d.Y = [ 0 bottom; 0 top; left top; left bottom];
            q_geom_d.N = [nx; ny];
            % Centre right
            q_geom_e.Y = [ right bottom; right top; 1 top; 1 bottom];
            q_geom_e.N = [nx; ny];

            % Upper left
            q_geom_f.Y = [0 top; 0 1; left 1; left top];
            q_geom_f.N = [nx; ny];
            % Upper centre
            q_geom_g.Y = [ left top; left 1; right 1; right, top];
            q_geom_g.N = [nx; ny];
            % Upper right
            q_geom_h.Y = [ right top; right 1; 1 1; 1, top];
            q_geom_h.N = [nx; ny];

            % Collect all
            q_geom = [q_geom_a, q_geom_b, q_geom_c, q_geom_d, q_geom_e, q_geom_f, q_geom_g, q_geom_h];
        end
        
        % Now build MS from number of quadrilaterals
        for q_i = 1:length(q_geom)
            shapes(q_i).shape = 'Quadrilateral';
            shapes(q_i).geom = q_geom(q_i);
        end
        %---------------------------------------------------------------
        % Create MS
        %---------------------------------------------------------------
        MS_opts = struct('doDiff',false,'doInd',false,'doInt',true,'doInterp',false);        % Optimised
        %
        MS = MultiShape(shapes,[], MS_opts);
        clear shapes        % Reset this variable

        % Plot Multi shapes
        %%MS.PlotGrid;
        %%pause(0.2)
        %%clf(figure(1))

        
        % Extract Integration weights and Interpolation
        Int = MS.Int;

        %% ---------------------------------------------------------------
        % Test MS
        % ---------------------------------------------------------------
        % Build interpolator from original box to MS
        Interp = aBox.InterpolationMatrix_Pointwise(MS.Pts.y1_kv, MS.Pts.y2_kv);
        % Evaluate kernel outside the point (y1,y2)
        K_diff_x = Kernel(y1 - MS.Pts.y1_kv, y2 - MS.Pts.y2_kv);
        % Convolution weights to act on f
        Conv(i,:) = (Int .* K_diff_x') * Interp;

        %% Time for some numerical experiments

        % K▣ ★ f at point i = ★ at MS ◰
        % Compute ★ at MS
        Conv_MS = Conv(i,:) * fBox;

        % Compare
        errors(i) = abs(t_Full(i) - Conv_MS);
        experiment(i) = abs(all_sq(i) - Conv_MS);
        full_thing(i) = abs(all_sq(i) - Conv_MS - smaller(i));

        %% Note:
        % It seems that relative errors are quite bad, since the exact
        % integrals are always below 1
        %errors(i) = abs(t_Full(i,:) - (Conv_capBox + Conv_MS)) ./ max(abs(t_Full(i,:)), 1e-6);
        
        % plot( errors ./ (abs(t_Full) + 1e-1 ))
        % aBox.plot( errors(:, logical([1 0 1]) ) ./ max( abs(t_Full(:, logical([1 0 1]) ) ), 1e-5) );
        % aBox.plot( errors ./ max( abs(t_Full), 1e-5) );
    end
    
    
    
end