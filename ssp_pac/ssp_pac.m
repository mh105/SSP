function [regressed_pac, SS_tot, X_t_tot] = ssp_pac (modelOsc, slowID, fastID,Kappa_tot, Nlag,Nconst)
%% SSP_PAC computes phase amplitude coupling between a fast and a slower oscillation, (and resampled estimations) after oscillation decomposition
%   INPUTS : modelOsc, structure obtained fromn SSP_DECOMP
%            Kappa   : number of resampled time series for CI estimation
%            slowID  : corresponding to slower oscillation in x_t_n
%            fastID  : corresponding to faster oscillation in x_t_n
%            Nlag    : max lag used to resample oscillations time series
%            Nconst (optionnal) : 1: constraint on K_mod, >1 linear, constraint, 0: none

%   OUTPUTS : latent_pac : 3 x (Kappa+1) x Nwin  Set of 3D Vectors
%   Enconding Strength K^mod_t and phase Phi_mod_t_k of the oscillation for
%   each window and resampled series
%             SS_tot 3 x (Kappa+1) x Nwin sum of squares : [ESS  TSS  RSS]'
%             X_t_tot     : 2 x Nosc x nT x (Kappa+1) Resampled oscillations (might be large)


%% Graps Decomposition Parameter and Init
staPointId = modelOsc.startPointId;
endPointId = modelOsc.endPointId;
lengthWind   = max(endPointId-staPointId+1);
p = size(modelOsc.res{1,1}.x_t_n,1);
Fs = modelOsc.Fs;

Kappa1 = Kappa_tot(1,1);
if Kappa1==0
    Kappa2=0;
    
elseif Kappa1>0
    
    if size(Kappa_tot,1)<size(Kappa_tot,2)
        Kappa_tot=Kappa_tot';
    end
    
    if size(Kappa_tot,1)==1
        Kappa2 = 0;
    else
        Kappa2 = Kappa_tot(2,1);
    end
    
end

if nargin<5
    Nlag = inf;  % Max Lag used to resample the time series

end

if nargin<6
    Nconst = 1;
end

if nargout>2
    saveXt = 1;
    warning('Outputing X_t_tot might result in unnecessary memory allocation')
    X_t_tot= zeros(p,lengthWind+1, length(staPointId), Kappa1+1);
else
    saveXt=0;
end

% Used for regression
M_t = @(Phi) [ones(length(Phi),1) cos(Phi), sin(Phi)];


regressed_pac  = zeros(3,(1+Kappa1)*(1+Kappa2),length(staPointId));
SS_tot = zeros(3,1+Kappa1,length(staPointId));
multPhase_tot = zeros(1,length(staPointId));

for ii=1:length(staPointId) % Parallel loop on each time window 
    %% Resample time series
    disp(['Window: ',num2str(ii),'/' ,num2str(length(staPointId)), ' resampling'])
    
    
    % From harmonic/oscillation decomposition
    Yreg_t     = modelOsc.res{1,ii}.y;
    x_t_n = modelOsc.res{1,ii}.x_t_n(:,1:end);
    
    
    paramsTot=modelOsc.res{1,ii}.model_prams;
    Phi  =paramsTot.Phi;
    Q    =paramsTot.Q;
    R    =paramsTot.R;
    mu   =paramsTot.mu;
    sigma=paramsTot.sigma;
    Nosc = size(Phi,1)/2;
    
    % X_ii_n ~ (2 x Nosc) x Nwind+1 x Kappa
    X_ii_n = zeros(size(x_t_n,1)  , size(x_t_n,2) , Kappa1+1);
    X_ii_n(:,:,1) = x_t_n;
    
    % Resample if necessary
    if Kappa1>0
        %tic
        [x_t_n,Pstruct] = get_Pts_n_KF(Yreg_t,Phi,Q,R,mu,sigma,Nlag);
        %toc
        MU_tmp=x_t_n(:,1:end);
        MU=MU_tmp(:);
        PP=Pstruct.PP;
        
        % Make sure that PP is semi positive definite before sampling with it
        [CholDec,err] = cholcov(Pstruct.PP);
        if err==0
            X_ii_n_tmp=mvnrnd(MU, PP, Kappa1,CholDec)';
        else
            PP=nearest_spd(PP);
            X_ii_n_tmp=mvnrnd(MU, PP, Kappa1)';
        end
        X_ii_n(:,:,2:end) = reshape( X_ii_n_tmp , [size(x_t_n,1)  , size(x_t_n,2),Kappa1 ]);
        
        % Gather resampled series if necessary
        if saveXt
            X_t_tot(:,:,ii,:)=X_ii_n;
        end
    end
    
    
    %% Estimate Phase Amplitude Coupling
    disp(['Window: ',num2str(ii),'/' ,num2str(length(staPointId)), ' pac'])
    
    % Take oscillation corresponding to ID slow and fast
    w_i = zeros(1,Nosc);
    for nosc=1:Nosc
        [~,w_tmp]=get_rotmat_pam(Phi(2*(nosc-1)+1:2*nosc,2*(nosc-1)+1:2*nosc));
        w_i(1,nosc)=w_tmp;
    end
    
    [~,sortedID]= sort(abs(w_i));
    w1 = w_i(sortedID(slowID));
    w2 = w_i(sortedID(fastID));
    
    slowIDcur=sortedID(slowID);
    fastIDcur=sortedID(fastID);
    
    % (If w_i<0, it changes the sign of the modulation phase)
    multPhase_tot(1,ii)= sign(w1)*sign(w2);
    
    % Grasp slow oscilation phase and fast oscilation amplitude
    % Phase/Amp ~ 1 x Nwind+1 x Kappa || X_ii_n ~ (2 x Nosc) x Nwind+1 x Kappa
    phase_slow = atan2(X_ii_n(2*slowIDcur,2:end,:), X_ii_n(2*slowIDcur-1,2:end,:))*sign(w1);
    ampli_fast = sqrt(X_ii_n(2*fastIDcur,2:end,:).^2 + X_ii_n(2*fastIDcur-1,2:end,:).^2);
    
    latent_pac_buffer = zeros(3,(1+Kappa1)*(1+Kappa2));
    SS_buffer         = zeros(3,1+Kappa1);
    
    
    for kk = 1:Kappa1+1
        
        ampli_fast_cur = ampli_fast(1,:,kk)';
        phase_slow_cur = phase_slow(1,:,kk)';
        
        Yreg_t = ampli_fast_cur;       % New observation ~ n x 1
        Xreg_t = M_t (phase_slow_cur); % Regression Matrix ~ n x 3
        
        
        % Priors Params : a, nu, b, V
        A_mean_t = mean(Yreg_t);
        a_prior  = [A_mean_t;0; 0];
        nu_prior = 3;
        b_prior  = 1;
        V_prior = (1/A_mean_t) * diag([1, 4/pi, 4/pi])* b_prior * nu_prior/(nu_prior-2);
        
        
        
        % Posterior Params : a, nu, b, V
        [a_post,nu_post, V_post,b_post] = get_post_param(Yreg_t,Xreg_t,a_prior,nu_prior,V_prior,b_prior);
        
        
        % Define linear constraintss if necessary
        if Nconst>0
            theta_Const = linspace(-pi,pi, Nconst)';
            f = zeros(Nconst,1); 
            W = [ones(Nconst,1) -cos(theta_Const) -sin(theta_Const)];
        end
        
        
        
        latent_pac_ii = zeros(3,1+Kappa2);
        
        % Grasp estimate
        if Nconst==0
            latent_pac_ii(:,1) = a_post;
        elseif Nconst>0
            
            % Define an upper bound for A_t 
            f_slow = Fs*abs(w1)/(2*pi);
            T_so_hat = min(10,1/f_slow);
            T_wind = length(Yreg_t)/Fs;
            
            
            
            Cmax = 1/(1-T_so_hat/(T_wind*pi));
            %assert(Cmax>=1)
            A_0_max = max([a_post(1,1), sqrt(a_post(2,1)^2+a_post(3,1)^2) ,A_mean_t ]) * Cmax;
           
            
            
            % Grasp pdf
            dbeta=0.025 * A_0_max ; % Sample resolution
            betamLim = A_0_max;
            res = 20;  % Resolution multiplicator
            [pdf_post_const,~,~,~,beta_0,beta_1,beta_2] = mvt3_const_pdf(dbeta,betamLim ,a_post, nu_post, V_post, b_post, W,f);

    
            
            % Return the mode
            [~,IndexMode] = max(pdf_post_const(:));
            [beta_mode_id_0,beta_mode_id_1,beta_mode_id_2] = ind2sub(size(pdf_post_const),IndexMode);
            b0 = beta_0(beta_mode_id_0);
            b1 = beta_1(beta_mode_id_1);
            b2 = beta_2(beta_mode_id_2);
            
            latent_pac_ii(:,1) =[b0;b1;b2];
            
        end
        
        % Resample the linear regression
        if Kappa2>0
            if Nconst==0
                beta_samp_post = mvt3rnd_hs(a_post,(b_post * inv(V_post)),nu_post, Kappa2);
            elseif Nconst>0
                beta_samp_post  = sample_from_pdf3(beta_0,beta_1,beta_2,pdf_post_const,Kappa2,res);
         
            end 
            latent_pac_ii(:,2:end) = beta_samp_post; % ~ 3 x (Kappa2)

        end
        

        % Yreg_t ~ WindowLength || Xreg_t ~ WindowLength x 2 || Yhat_t ~ WindowLength x (Kappa2+1)
        Yhat_t = Xreg_t * latent_pac_ii;
   
        % Sums of square
        ESS = (Yhat_t(:,1) - mean(Yreg_t))'*(Yhat_t(:,1) - mean(Yreg_t)); % ~ 1 x (1+Kappa2)
        RSS = (Yreg_t-Yhat_t(:,1))'*(Yreg_t-Yhat_t(:,1))        ;% ~ 1 x (1+Kappa2)
        TSS = (Yreg_t - mean(Yreg_t))'*(Yreg_t - mean(Yreg_t)) ;% ~ 1 x 1
        
        IdK_cur = 1 + (kk-1)*(1+Kappa2);
        IdK_nm1 = 1 + (kk)*(1+Kappa2)-1;
        
        latent_pac_buffer(:,IdK_cur:IdK_nm1) = latent_pac_ii;
        SS_buffer(:,kk) = [ESS;TSS;RSS];
        
    end
    
    regressed_pac(:,:,ii) = latent_pac_buffer;
    SS_tot    (:,:,ii) = SS_buffer;
    
    
    
    
end

end








