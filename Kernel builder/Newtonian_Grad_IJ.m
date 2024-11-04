function N = Newtonian_Grad_IJ(Points)
    % Points should be a [n,2] matrix of coordinates
    N.I = @dI;
    %aBox.plot(Newtonians.I(Points))    % Has the right shape and numerical range

    N.J = @dJ;
    %aBox.plot(Newtonians.J(Points))    % Has the right shape and numerical range
    
    % We return two possible outputs based on the number of inputs
    switch nargin
        case 1
            N.dI = dI(Points);
            N.dJ = dJ(Points);
        otherwise
            N;
    end
end

function dI_All = dI(Points)
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
    
    dI_All = zeros(2*n,1); % We expect either a point or a cloud of points
    
    % Split data into two cases
    ZP = ~all((Points == 0), 2);    % Identifies points without a zero coordinate
    % The function is 0 at any point with a zero coordinate, by construction we do nothing
    
    % Evaluate at points in mesh without a zero coordinate
    x = Points(ZP,1);
    y = Points(ZP,2);
    % Derivative in x
    I = ( 0.5 * y .* log(x.^2 + y.^2) ) - y + ( x .* atan2(y, x) );
    I = I / (2*pi);
    
    dI_All(ZP) = I;

    % Derivative in y
    I = ( 0.5 * x .* log(x.^2 + y.^2) ) - x - ( y .* atan2(y, x) );
    I = (I / (2*pi)) + 0.25 * y ;

    % Store results
    dI_All([ logical(zeros(n,1)) ;ZP]) = I;

end

function dJ_All = dJ(Points)
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

    % Create shifts
    Shift_x = zeros(size(Points));    Shift_y = zeros(size(Points));
    Shift_x(:,1) = -1;                Shift_y(:,2) = -1;
    Shift_xy = 1 - Points;

    Sign_x = ones(2*n,1);        Sign_y = ones(2*n,1);
    Sign_x(1:n) = -1;            Sign_y(n+1:end) = -1;

    
    dJ_All = dI(Points) + Sign_x .* dI(abs(Points + Shift_x)) + Sign_y .* dI(abs(Points + Shift_y)) - dI(Shift_xy);
end