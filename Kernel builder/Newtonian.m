function N = Newtonian(Points, epsilon)
    % Points should be a [n,2] matrix of coordinates
    N.I = @Conv_I;
    %aBox.plot(Newtonians.I(Points))    % Has the right shape and numerical range

    N.J = @Conv_J;
    %aBox.plot(Newtonians.J(Points))    % Has the right shape and numerical range

    N.G = @Conv_G;
    %aBox.plot(Newtonians.G(Points, 1e-4));
    
    % We return three possible outputs based on the number of inputs
    switch nargin
        case 1
            N.JN = Conv_J(Points);
        case 2
            N.NJ = Conv_J(Points);
            N.NG = Conv_G(Points, epsilon);
            N.NI = N.NJ - N.NG;
            % aBox.plot(Newtonians.NI);   % Exact integral [.5,.4 r so cool]
        otherwise
            N;
    end
end

function I_All = Conv_I(Points)
    % Function only suitable for 2D inputs
    if ~ismatrix(Points) || isscalar(Points)
        error('Input array is not 2 dimensional.');
    end

    % Fix dimension in case points are transposed
    [n, m] = size(Points);

    if m ~= 2
        Points = Points.';
        n = m;
    end
    
    I_All = zeros(n,1); % We expect either a point or a cloud of points
    
    % Split data into two cases
    ZP = ~any((Points == 0), 2);    % Identifies points without a zero coordinate
    % The function is 0 at any point with a zero coordinate, by construction we do nothing
    
    % Evaluate at points in mesh without a zero coordinate
    x = Points(ZP,1);
    y = Points(ZP,2);
    I = ( x .* y .* log(x.^2 + y.^2) ) + ( (x.^2 - y.^2) .* atan2(y, x) ) - (3 * x .* y);
    I = I / (4*pi);
    I = I + 0.125 * y.^2;
    
    I_All(ZP) = I;
end

function J_All = Conv_J(Points)
    % Function only suitable for 2D inputs
    if ~ismatrix(Points) || isscalar(Points)
        error('Input array is not 2 dimensional.');
    end
    % Fix dimension in case points are transposed
    m = size(Points,2);
    if m ~= 2
        Points = Points.';
    end

    % Create shifts
    Shift_x = zeros(size(Points));    Shift_y = zeros(size(Points));
    Shift_x(:,1) = -1;                Shift_y(:,2) = -1;
    Shift_xy = 1 - Points;

    
    J_All = Conv_I(Points) + Conv_I(abs(Points + Shift_x)) + Conv_I(abs(Points + Shift_y)) + Conv_I(Shift_xy);
end

function G_All = Conv_G(Points, epsilon)
    % Function only suitable for 2D inputs
    if ~ismatrix(Points) || isscalar(Points)
        error('Input array is not 2 dimensional.');
    end
    % Fix dimension in case points are transposed
    m = size(Points,2);
    if m ~= 2
        Points = Points.';
    end

    % Create shifts
    Shift_x = zeros(size(Points));    Shift_y = zeros(size(Points));
    Shift_x(:,1) = -1;                Shift_y(:,2) = -1;
    Shift_xy = min(1 - Points, epsilon);
    
    
    G_All = Conv_I(min(Points, epsilon)) + ...
            Conv_I(min(abs(Points + Shift_x), epsilon)) + ...
            Conv_I(min(abs(Points + Shift_y), epsilon)) + Conv_I(Shift_xy);
end