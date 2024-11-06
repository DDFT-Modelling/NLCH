function Phi_t = NL_CH_Integrator_Source_Constant(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF, g_tilde, dnsigma, TimeLine, Tmin )

% Check initial time
if nargin <= 11
  Tmin = 0.0;
end

    % Retrieve differential operators in space
    grad  = Diff.grad;       div = Diff.div;
    bound = Ind.bound;    normal = Ind.normal;
    N = size(Conv,1);  Lap_Space = Diff.Lap;
    y1 = aBox.Pts.y1_kv;
    y2 = aBox.Pts.y2_kv;

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-flux BC) there.
    
    mM        = ones([N,1]);
    mM(bound) = 0;
    
    
    %opts = odeset('RelTol',10^-4,'AbsTol',10^-4,'Stats','on', 'Jacobian', @JacobiFun,'Mass',diag(mM));
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Stats','on', 'Mass',diag(mM),'Vectorized','on');
    opts = odeset('RelTol',10^-7,'AbsTol',10^-7,'Stats','on',   ...
                  'Mass',diag(mM),'MStateDependence','none', 'MassSingular','yes', ...      % Additional options save a bit of time
                  'InitialStep',1e-4, ...       % Allowed to solve the 40 case with more accuracy
                  'Vectorized','on', ...        % Reduced computation times 10 fold
                  'BDF','on' );                 % Had more successful steps when on, I will consider it again later  
    %                  'InitialSlope', dy0, ...      % daeic12 was struggling with this, which makes sense for time 0
    %dy0 = rhs(Tmin,phi_ic);
    %opts.InitialSlope = dy0;
    %opts.Jacobian = @JacobiFun;

    Tmax = TimeLine.yMax;
    h    = (Tmax - Tmin)/99;

    tic
    %[~,Phi_t] = ode15s(@rhs, TimeLine.Pts.y(2:end), phi_ic, opts);
    %[~,Phi_t] = ode15s(@rhs, Tmin:h:Tmax, phi_ic, opts);
    [~,Phi_t] = ode15s(@rhs, TimeLine.Pts.y, phi_ic, opts);
    toc
   

    % Define ODE in time
    function dydt = rhs(t,phi)    
        % Compute the first derivative of mu and apply nonlocal operator
        flux = getFlux(phi);
        % Obtain Laplacian of mu and add the source term
        dydt = (div * flux) + g_tilde;
    
        % no-flux BCs: gives that j.n = mu.n on boundary
        dydt(bound,:) = (normal * flux) - dnsigma;
    end

    function f = getFlux(phi)
        % the flux is the derivative of the difference between F' and the
        % nonlocality
        %
        % Approximate derivative: [n]
        % z = phi;
        % z(z > 1) =  1.0 - 1e-5;
        % z(z <-1) = -1.0 + 1e-5;

        f1 = dF(phi) - (Conv * phi + Conv_D .* phi);
        f = grad * f1;
    end

    function drhs = JacobiFun(t,phi)
        % Spatial part
        Jmu = scalarOperator(ddF(phi) - Conv_D) - Conv;
        drhs = (Lap_Space * Jmu);
        % Boundary
        drhs(bound,:) = normal * grad * Jmu;
    end

    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end

end
% Life lesson «If you have a code that works already: don't touch it.»