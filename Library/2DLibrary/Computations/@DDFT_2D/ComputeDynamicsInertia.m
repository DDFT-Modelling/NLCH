function ComputeDynamicsInertia(this)

%% Solves
% 
% $$\frac{\partial \varrho}{\partial t} = \nabla \cdot ( \varrho \nabla \mu )$$
% 
% in computational variable x: $\varrho = e^{(x-Vext)/(k_B T)}$
% 
% $$\frac{\partial x}{\partial t} = k_B T \Delta \mu + (\nabla x - \nabla V)\cdot \nabla \mu$$
% 

    % Initialization
    optsPhys    = this.optsPhys;
    optsNum     = this.optsNum;
        
    if(isfield(optsPhys,'sigmaS'))    
        R       = optsPhys.sigmaS/2;
    else
        R       = [];
    end
    
    if(isfield(optsPhys,'Fext'))    
        Fext    = optsPhys.Fext;
    else
        %Fext    = [0,0];
        Fext    = [];
    end
    
    if(isfield(optsPhys,'BCWall_U'))
        BCWall_U = optsPhys.BCWall_U;
    else
        BCWall_U = [];
    end
    
    if(isfield(optsPhys,'mS'))
        mS = optsPhys.mS;
    else
        mS = 1;
    end
    mInv        = mS.^(-1);
    Diff        = this.IDC.Diff;
    
    if(isfield(optsPhys,'viscosity'))
        doVisc = true;
    else
        doVisc = false;
    end
        
    if(isstruct(optsNum.plotTimes))
        tI         = optsNum.plotTimes.t_int;
        t_n        = optsNum.plotTimes.t_n;
        
        plotTimes   = tI(1)+(tI(2)-tI(1))*(0:1:(t_n-1))/(t_n-1);
    else
        plotTimes   = optsNum.plotTimes;
    end
    
    %!!!!!!!!!!!!!!!!!!!!!!
    div = [Diff.Dy1' Diff.Dy2'];
    
    kBT         = optsPhys.kBT;
    gammaS      = optsPhys.gammaS;
    nSpecies    = optsPhys.nSpecies;
    M           = this.IDC.M;
    Vext        = this.Vext;
    Vext_grad   = this.Vext_grad;
    IntMatrHI   = this.IntMatrHI;
    IntMatrFex  = this.IntMatrFex;
    Int_of_path = this.Int_of_path;
    Conv        = this.IntMatrV2;
    Ind         = this.IDC.Ind;        
    getFex      = str2func(['Fex_',optsNum.FexNum.Fex]);    
    doHI        = this.doHI;
    markVinf    = (Vext == inf);
    if(strcmp(this.IDC.polar,'polar'))
        polarShape = true;
    else
        polarShape = false;
    end
    
    subArea     = this.subArea;
    
	I           = eye(M);  
    eyes        = [I I];
    
    PtsCart = this.IDC.GetCartPts();
	y1S     = repmat(PtsCart.y1_kv,1,nSpecies); 
    y2S     = repmat(PtsCart.y2_kv,1,nSpecies);

    ythS    = repmat(this.IDC.Pts.y2_kv,1,nSpecies);
    
    
	x_ic = GetInitialCondition(this);
	mu   = this.mu;    
    
    fprintf(1,'Computing dynamics ...'); 

    if(isfield(optsNum,'PlotArea'))
        optsNumT = rmfield(optsNum,'PlotArea');
    else
        optsNumT = optsNum;
    end
    ignoreList = 'PlotAreaCart';
    
    [this.dynamicsResult,recEq,paramsEq] = DataStorage('Dynamics',...
                                    @ComputeDDFTDynamics,...
                                    v2struct(optsNumT,optsPhys),[],[],ignoreList); %true      
                        
    this.dynamicsResult.t = plotTimes;            
    this.FilenameDyn      = paramsEq.Filename;            
                     
    function data = ComputeDDFTDynamics(params,misc)        
       
        mMx              = ones(M,nSpecies);        
        mMx(markVinf)    = 0;
        mMuv             = ones(2*M,nSpecies);
        if(doVisc)
            mMuv([Ind.finite;Ind.finite],:)   = 0;
        else
            mMuv([Ind.finite1;Ind.finite2],:) = 0;
        end
        mM = [mMx;mMuv];
        mM = mM(:);

        x_ic = [x_ic;zeros(2*M,nSpecies)]; % pad with zero velocity
        
        opts    = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([ones(nSpecies,1);mM]));    
        [~,X_t] = ode15s(@dx_dt,plotTimes,[zeros(nSpecies,1);x_ic(:)],opts);   
                       
        nPlots    = length(plotTimes);

        accFlux   = X_t(:,1:nSpecies);
        X_t       = X_t(:,nSpecies+1:end)';
        X_t       = reshape(X_t,3*M,nSpecies,nPlots);
        
        Y_t       = X_t(1:M,:,:);
        UV_t      = X_t(M+1:3*M,:,:);
        X_t = Y_t;
         
        rho_t     = zeros(M,nSpecies,nPlots);
        flux_t    = zeros(2*M,nSpecies,nPlots);
        V_t       = zeros(M,nSpecies,nPlots);   

        for i = 1:length(plotTimes)
            rho_t(:,:,i)  = exp((Y_t(:,:,i)-Vext)/kBT);
            flux_t(:,:,i) = GetFlux([X_t(:,:,i);UV_t(:,:,i)],plotTimes(i));           
            V_t(:,:,i)    = Vext + getVAdd(y1S,y2S,plotTimes(i),optsPhys.V1);
        end
       
        data       = v2struct(X_t,UV_t,rho_t,mu,flux_t,V_t); %IntMatrFex
      %  data.shape = this.IDC;
        if(this.doSubArea) 
            data.Subspace = v2struct(accFlux); %subArea
        end
    end   
    function dxdt = dx_dt(t,x)
        
        % ignore first row of entries. This is mass in subsystem 
        x        = x(nSpecies+1:end);              
        
        x        = reshape(x,3*M,nSpecies);
        
        y  = x(1:M,:);
        u  = x(M+1:2*M,:);
        v  = x(2*M+1:3*M,:);
        uv = x(M+1:3*M,:);
        
        % convective term; C = v.grad
        C  = zeros(M,M,nSpecies);
        Cu = zeros(M,nSpecies);
        Cv = zeros(M,nSpecies);
        for iSpecies = 1:nSpecies
            C(:,:,iSpecies)  = sparse(1:M,1:M,u(:,iSpecies))*Diff.Dy1 + sparse(1:M,1:M,v(:,iSpecies))*Diff.Dy2;
            Cu(:,iSpecies)   = C(:,:,iSpecies)*u(:,iSpecies);
            Cv(:,iSpecies)   = C(:,:,iSpecies)*v(:,iSpecies);
        end
        
        mu_s     = GetExcessChemPotential(y,t,mu);        
        mu_s(markVinf) = 0;
        
        h_s1    = Diff.Dy1*y - Vext_grad(1:M,:);
        h_s2    = Diff.Dy2*y - Vext_grad(1+M:end,:);
                
        h_s1(markVinf)  = 0; 
        h_s2(markVinf)  = 0; 
        
        dydt     = -kBT*(Diff.div*uv) - (u.*h_s1  + v.*h_s2);
        duvdt    = - [Cu;Cv] - gammaS*uv - mInv.*(Diff.grad*mu_s);
                
        if(~isempty(Fext))
            duvdt = duvdt + [Fext(1)*ones(M,1);Fext(2)*ones(M,1)];
        end
                
        if(doVisc)
            rho_s    = exp((y-Vext)/kBT);
            [zeta,dzeta_drho] = BulkViscosity(rho_s,optsPhys);
            [eta,deta_drho]   = ShearViscosity(rho_s,optsPhys);
            
            etaR       = repmat(eta./rho_s,2,1);
            zetaR      = repmat(zeta./rho_s,2,1);
            dzeta_drho = repmat(dzeta_drho,2,1);
            deta_drho  = repmat(deta_drho,2,1);
            
            duvdt = duvdt + etaR.*(Diff.LapVec*uv) ...
                          + (zetaR + etaR/3).*(Diff.gradDiv*uv);
            
            if((max(abs(dzeta_drho)) ~= 0) || (max(abs(deta_drho)) ~= 0))
                dydy1      = (Diff.Dy1*y)/kBT;
                dydy2      = (Diff.Dy2*y)/kBT;                
                
                duvdt      = duvdt + (dzeta_drho - 2/3*deta_drho).*repmat(Diff.div*uv,2,1).*[dydy1;dydy2];                                
                dudy2dvdy1 = (Diff.Dy2*u + Diff.Dy1*v);
                
                visc1      = 2*(Diff.Dy1*u).*dydy1 + dudy2dvdy1.*dydy2;
                visc2      = dudy2dvdy1.*dydy1 + 2*(Diff.Dy2*v).*dydy2;
                duvdt      = duvdt + [visc1;visc2].*deta_drho;
            end       
        end
        
        if(doHI)
            error('HI not implemented yet');
%             rho_s    = exp((x-Vext)/kBT);
%             rho_s    = [rho_s;rho_s];
%             gradMu_s = Diff.grad*mu_s;
%             HI_s     = ComputeHI(rho_s,gradMu_s,IntMatrHI);            
%             dxdt     = dxdt + kBT*Diff.div*HI_s + eyes*( h_s.*HI_s );  
        end
        
        duvdt([Ind.finite1 ; false(M,1)],:)  = Ind.normalFinite1*u;
        duvdt([false(M,1)  ; Ind.finite2],:) = Ind.normalFinite2*v;
        
        if(doVisc) %set velocities parallel to the interface 
            duvdt([Ind.finite2 ; false(M,1)],:)   = u(Ind.finite2) - WallBC(t,Ind,BCWall_U); 
            duvdt([false(M,1) ; Ind.finite1],:)   = v(Ind.finite1);
        end
        
        dydt(markVinf)         = y(markVinf) - x_ic(markVinf);

        dxdt = [dydt;duvdt];
        
        dxdt = [(Int_of_path*GetFlux(x,t))';dxdt(:)];
    end
    function mu_s = GetExcessChemPotential(x,t,mu)
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R) + ...
                                Fex_Meanfield(rho_s,Conv,kBT);
        
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu(iSpecies);
        end

        mu_s = mu_s + x + getVAdd(y1S,y2S,t,optsPhys.V1);
    end   
    function flux = GetFlux(x,t)
        x        = reshape(x,3*M,nSpecies);        
        y  = x(1:M,:);
        uv = x(M+1:3*M,:);
        
        rho_s = exp((y-Vext)/kBT);       
        flux  = [rho_s;rho_s].*uv;                              
        if(polarShape)
            %then transform to cartesian corrdinates
            flux = GetCartesianFromPolarFlux(flux,ythS);
        end
    end

end