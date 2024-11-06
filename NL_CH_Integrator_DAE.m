function Rho_t = NL_CH_Integrator_DAE(phi_ic,  aBox, Diff, Ind,  Conv, Conv_D, eta, TimeLine)
    % Solves the Nonlocal Cahn–Hilliard equation with a source 
    % No flux boundary conditions are enforcing using own DAE

    %% Scalar functions
    % F  = @(s) (1+s) .* log(1+s) + (1-s) .* log(1-s);      % Logarithmic potential
    dF = @(s) log(1+s) - log(1-s);                          % Derivative of logarithmic potential
    ddF = @(s) 2 ./ ( 1 - s.^2 );
    Q  = @(s) tanh(0.5 * s);                                % Inverse of dF
    dQ = @(s) 0.5 * sech(0.5 * s).^2;

    %% Retrieve differential operators in space
    grad  = Diff.grad;       div = Diff.div;
    bound = Ind.bound;    normal = Ind.normal;
    N = size(Conv,1);  Lap_Space = Diff.Lap;
    y1 = aBox.Pts.y1_kv;
    y2 = aBox.Pts.y2_kv;
    N_bound = size(aBox.Ind.normal,1);
    % Scale convolutions if needed
    Conv   = eta * Conv;
    Conv_D = eta * Conv_D;

    %% Select submatrices and invert A at boundary indices
    % This is A:
    normal_grad = normal * grad;
    % Multiplying A times the discretisation of the convolution, we obtain B
    B_No_Flux = normal_grad * (Conv + diag( Conv_D )  );    
    % We can see that |B|_2 < 1 (40: 0.0929) via norm(B_No_Flux)
    %norm(B_No_Flux)
    
    % Obtain the inverse of A at the boundary:
    A_bound_inv = normal_grad(:,bound)^-1;
    % To find C just multiply the previous two matrices at the indices of the boundary
    C_bound = A_bound_inv * B_No_Flux(:,bound);
    % We find |C|_2 < 1 (40: 6e-7)
    norm(C_bound)
    
    %% Auxiliary functions for fixed point
    R   = @(u,d)     Q(C_bound * u + d);
    NR  = @(u,d) u - Q(C_bound * u + d);
    dNR = @(u,d) eye(N_bound) - diag( dQ(C_bound * u + d) ) * C_bound;


    %% Solve the PDE, the RHS of which is given below
    % There is no need to define a mass matrix to enforce BC

    
    Tmax = 1; %10.0; %TimeLine.Pts.y(2);
    %Tmax = 8;

    %opts = odeset('RelTol',10^-4,'AbsTol',10^-4,'Stats','on', 'Jacobian', @JacobiFun,'Mass',diag(mM));
    opts = odeset('RelTol',10^-7,'AbsTol',10^-7,'Stats','on','InitialStep',1e-3, 'MaxStep',Tmax * 0.9);
    %opts.Jacobian = @JacobiFun;
    opts.Vectorized = 'on';
    
    %d = ode(ODEFcn=@rhs);          d.InitialTime = 0.0;              d.InitialValue = phi_ic(~bound);
    %d.RelativeTolerance = 1e-7;         d.Solver = "cvodesstiff";    %d.Jacobian = @JacobiFun;

    tic
    %[~,Phi_t] = ode45(@rhs, 0:TimeLine.yMax/100:TimeLine.yMax, phi_ic(~bound), opts);
    %[~,Phi_t] = ode15s(@rhs, 0:Tmax/99:Tmax, phi_ic(~bound), opts);
    [~,Phi_t] = ode15s(@rhs, TimeLine.Pts.y, phi_ic(~bound), opts);
    
    % Sundials (4xslow for easy problem)
    %sol = solve(d,0,Tmax);
    %fh = solutionFcn(d,0,Tmax);
    %Phi_t = fh(0:Tmax/99:Tmax)'; %sol.Solution'; 
    
    % Euler
    %Phi_t = phi_ic(~bound);
    %Phi_t = Phi_t' + Tmax * rhs(0,Phi_t);
    toc

    %   rhs(1e-15, 0.5*ones(N-N_bound,1));      % A wee test that takes 0.047151 s

    %% Now process the solution for it to also contain the domain
    %Rho_t = zeros(1,N);     % rn let's just interpolate once
    % 
    d_bound = A_bound_inv * ( B_No_Flux(:,~bound) * Phi_t' - normal_grad(:,~bound) * dF(Phi_t')  );
    % Find a fixed point
    rho_bound = R(R(R(zeros(N_bound,1), d_bound), d_bound), d_bound);
    sprintf('%.2e', norm( rho_bound - R(rho_bound, d_bound), 'fro') )


    %Fixed_Point = @(u) deal(NR(u,d_bound), dNR(u,d_bound));             % Run just this one every time to set the parameter
    % rho_bound = fsolve(Fixed_Point, zeros(size(d_bound)), ...
    %                            optimoptions('fsolve', 'Display', 'none','Algorithm','trust-region','SpecifyObjectiveGradient', true));
    Rho_t = zeros(size(d_bound,2), N);
    Rho_t(:,bound)  = rho_bound';
    Rho_t(:,~bound) = Phi_t;



    %% Define ODE in time
    function dydt = rhs(t,rho_inner)
        %t

        %% Interior-to-boundary
        % Compute constant vector
        d_bound = A_bound_inv * ( B_No_Flux(:,~bound) * rho_inner - ...
                                   normal_grad(:,~bound) * dF(rho_inner) );

        %Fixed_Point = @(u) deal(NR(u,d_bound), dNR(u,d_bound));             % Run just this one every time to set the parameter
        % Compute u at bound
        %rho_bound = fsolve(Fixed_Point, zeros(N_bound,1), ...
        %                       optimoptions('fsolve', 'Display', 'none','Algorithm','trust-region','SpecifyObjectiveGradient', true));
        %
        %plot(rho_bound)

        % Question: Does it work if I directly apply a fixed point
        % iteration?
        x_0 = zeros(N_bound, size(rho_inner,2));         % Initial point
        rho_bound = R(R(R( x_0, d_bound), d_bound), d_bound);
        if norm( rho_bound - R(rho_bound, d_bound), 'fro' ) > 1e-10
            rho_bound = R(R(rho_bound, d_bound), d_bound);
        end
        % It does, most of the time it is 1e-14!
        

        %% Consolidate full domain information
        % This changes to vectorise code
        rho = zeros(N,size(rho_bound,2));               % [Was (N,1) before vectorisation]
        rho(bound,:)  = rho_bound;                      % [Remove ,: for unvectorised]
        rho(~bound,:) = rho_inner;                      % [Remove ,: for unvectorised]

        %% Compute dynamics
    
        % Compute the first derivative of mu and apply nonlocal operator
        flux = getFlux(rho);
        % Obtain Laplacian of mu and add the source term
        dydt = div * flux;
    
        % no-flux BCs: gives that j.n = 0 on boundary
        % We don't need this extra operator:    Solved using DAE and discarded at the end
        % dydt(bound,:) = normal * flux;        % [Remove ,: for unvectorised]

        %% Project to interior of domain
        dydt = dydt(~bound,:);                          % [Remove ,: for unvectorised]
    end

    function f = getFlux(phi)
        % the flux is the derivative of the difference between F' and the nonlocality
        f1 = 1 * dF(phi) - 1 * (Conv * phi + Conv_D .* phi);
        f  = grad * f1;
    end

    %% Jacobian: Not really useful

    function drhs = JacobiFun(t,rho_inner)
        %% Interior-to-boundary
        % Compute constant vector
        d_bound = A_bound_inv * ( B_No_Flux(:,~bound) * rho_inner - ...
                                   normal_grad(:,~bound) * dF(rho_inner) );
        % Fixed point
        x_0 = zeros(N_bound, 1);         % Initial point
        rho_bound = R(R(R( x_0, d_bound), d_bound), d_bound);
        if norm( rho_bound - R(rho_bound, d_bound), 'fro' ) > 1e-10
            rho_bound = R(R(rho_bound, d_bound), d_bound);
        end
        %% Consolidate full domain information
        rho = zeros(N,1);
        rho(bound)  = rho_bound;
        rho(~bound) = rho_inner;
        
        %% Full Jacobian
        % Spatial part
        z = rho;
        z(z >  1.0 - 1e-3) =  1.0 - 1e-3;
        z(z < -1.0 + 1e-3) = -1.0 + 1e-3;

        Jmu = scalarOperator(ddF(z) - Conv_D) - Conv;
        drhs = (Lap_Space * Jmu);
        drhs = drhs(~bound,~bound);
        % Boundary:     We don't need this for the inner points
        %drhs(bound,:) = normal * grad * Jmu;
    end
    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end

end
% Life lesson «If you have a code that works already: don't touch it.»