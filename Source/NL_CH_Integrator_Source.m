function Phi_t = NL_CH_Integrator_Source(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, dF, ddF,    mu_tilde, g_tilde, TimeLine )

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

    %dy0 = rhs_2(That_Time, phi_ic' );
    dy0 = rhs(1e-4, phi_ic * 0.99 );
    dy0(bound) = 0.0;

    %%%phi_ic(bound) = normal * getFlux(phi_ic');

    %opts = odeset('RelTol',10^-4,'AbsTol',10^-4,'Stats','on', 'Jacobian', @JacobiFun,'Mass',diag(mM));
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Stats','on', 'Mass',diag(mM));
    % opts = odeset('RelTol',10^-5,'AbsTol',10^-5,'Stats','on', ...
    %               'Mass',diag(mM),'MStateDependence','none', 'MassSingular','yes', ...      % Additional options save a bit of time
    %               'InitialStep',1e-4, ...       % Allowed to solve the 40 case with more accuracy
    %               'Vectorized','on', ...        % Reduced computation times 10 fold
    %               'InitialSlope', dy0, ...      % daeic12 was struggling with this, which makes sense for time 0
    %               'BDF','on' );                 % Had more successful steps when on, I will consider it again later  

    tic
    [~,Phi_t] = ode15s(@rhs, 0:TimeLine.yMax/100:TimeLine.yMax, phi_ic, opts);
    toc
    %Phi_t = getFlux(phi_ic');
    %normal * (getFlux(phi_ic') - grad * mu_tilde(3,:)' )  % This is actually satisfied and equal to zero... should I set the boundary to zero then?

    

    % Setup and solve ODE in time
    %opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats, ...
    %        'NonNegative',1:3*N ); 
    %[~,Phi_t] = ode113(@rhs,outTimes,phi_ic,opts);

    % Define ODE in time
    function dydt = rhs(t,phi)
        % Compute interpolator at time t
        IP_Time = TimeLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate mu
        IP_mu = (IP_Time * mu_tilde)';
        % Interpolate g
        IP_g  = (IP_Time * g_tilde)';
    
        % Compute the first derivative of mu and apply nonlocal operator
        flux = getFlux(phi);
        % Obtain Laplacian of mu and add the source term
        %dydt = (div * flux)/(2*pi) + IP_g;
        dydt = (div * flux) + IP_g;
    
        % no-flux BCs: gives that j.n = mu.n on boundary
        %dydt(bound) = normal * (flux - grad * IP_mu/(2*pi));
        dydt(bound) = normal * (flux - grad * IP_mu);
        dydt = dydt(:);
    end

    function dydt = rhs_2(t,phi)
        % Compute interpolator at time t
        %IP_Time = TimeLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate mu
        %IP_mu = (IP_Time * mu_tilde)';
        % Interpolate g
        %IP_g  = (IP_Time * g_tilde)';

        % Evaluate exact solution
        if t < 1e-10
            tt = max(t,1e-10);
            fprintf('%.2e %.2e\n', t,tt);   % daeic12: It is evaluating everything many times at 0, why?
        end
        phi_tilde_2 = exp(-t) * sin( 2.0 * y1 * pi ) .* cos( 2.0 * y2 * pi );
        % Evaluate exact mu
        flux_phi_tilde_2 = getFlux(phi_tilde_2);
        mu_tilde_2 = dF(phi_tilde_2) - (Conv * phi_tilde_2 + Conv_D .* phi_tilde_2);
        % Evaluate exact g
        g_tilde_2 = -phi_tilde_2  - (Lap_Space * mu_tilde_2) / (4 * pi^2);   % This is an instance where computing div * grad vs Lap did make a difference

        %%% Error 
        %err_mu = max( abs(flux_phi_tilde_2 - grad * IP_mu/(2*pi)))
        %err_g  = IP_g - g_tilde_2;
        %[y1(ind),y2(ind)]
        %aBox.plot(flux_phi_tilde_2(1:end/2));
        %drawnow

        % Compute the first derivative of mu and apply nonlocal operator
        flux = getFlux(phi);
        % Obtain Laplacian of mu and add the source term
        dydt = (div * flux)/(2*pi) + g_tilde_2;
    
        % no-flux BCs: gives that j.n = mu.n on boundary
        dydt(bound,:) = normal * (flux - flux_phi_tilde_2);      % Vectorized input needs (bound,:), else can just use (bound)      %% [this is actually very close to 0!!]

        % Initial direction might be tricky when the data is not separated
        TF = isnan(dydt);
        if sum(TF) > 1
            [sum(TF), sum(TF(bound))]
        end

        %dydt = dydt(:);
    end
    
    function f = getFlux(phi)
        % the flux is the derivative of the difference between F' and the
        % nonlocality
        f1 = dF(phi) - (Conv * phi + Conv_D .* phi);
        f = grad * f1; %/(2*pi);
    end

    function drhs = JacobiFun(t,phi)
        % Spatial part
        Jmu = scalarOperator(ddF(phi) - Conv_D) - Conv;
        drhs = (Lap_Space * Jmu)/(4*pi^2);
        % Boundary
        drhs(bound,:) = normal * grad * Jmu/(2*pi);
    end

    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end

end
% Life lesson «If you have a code that works already: don't touch it.»