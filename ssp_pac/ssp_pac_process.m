function processed_ssp_pac = ssp_pac_process(regressed_pac, SS_tot,doArFiltering,tVect,Nlinconst,name)
%% SSP_PAC_PROCESS gathers results from modulation linear regression (and resampling) and filter them if necessary
% INPUTS:
%         - regressed_pac~ 3 x Kappa+1 x Nwind where Kappa is the number of
%         drawn samples and Nwind the number of sample windows
%         - SS_tot sums (ESS,TSS,RSS) of square ~ 3 x Kappa1+1 x Nwind
%         - doArFiltering, binary variable, use or not a 2nd state space to
%         filter regressed_pac
%         - tVect, time vector
%         - Nlinconst, constrain or not the linear regression (so that 0<K_mod<1)
%
% OUTPUTS:
%         - processed_ssp_pac, a structure containing estimated parameters K_mod, Phi_mod,R2, filtered parameters
%         and credible intervals.
%         Where, for 1 window: A_fast_t = A_0 (1 + K_mod cos(Phi_mod + Phi_slow_t))


if nargin <5
    Nlinconst = 0;
end

if nargin<6
   name='' ;
end

processed_ssp_pac = struct();
%% Grasp Fit parameters
%disp('Grasp Params')
Kappa   = size(regressed_pac,2)-1;
Kappa1  = size(SS_tot,2)-1;
Kappa2  = (Kappa+1)/(1+Kappa1)-1;
dim_pac = size(regressed_pac,1);
Nwind   = length(tVect);

% Fraction used to construct CI
alpha_ci = 0.05;

%% Goodness of fit
%disp('Goodness of fit')
ESS = reshape(SS_tot(1,:,:), Kappa1+1, Nwind);
TSS = reshape(SS_tot(2,:,:), Kappa1+1, Nwind);
RSS = reshape(SS_tot(3,:,:), Kappa1+1, Nwind);

R2 =  ESS./TSS;
R2_cen = mean(R2,1 );
R2_inf = quantile(R2, alpha_ci,1);
R2_sup = quantile(R2, 1-alpha_ci,1);


%%  Construct CI on regression Vector
%disp('Construct CI on regression Vector')

% Define a distance between pac vector and Kth resampled series ||pac -pac_k||^2
dlatent_pac = regressed_pac(:,1,:)-regressed_pac ;
normDiff = sum(dlatent_pac.^2,1);
[~, best_latent] = sort(normDiff, 2);

% Keep (1-alpha)% of samples
Nkept = max(floor(Kappa*(1-alpha_ci)),1);
regressed_pac_kept = zeros(dim_pac,Nkept,Nwind);
for tt = 1:Nwind
    regressed_pac_kept(:, :,tt) = regressed_pac(:,best_latent(1,1:Nkept,tt),tt);
end


%% Get phase amplitude coupling params K, Phi and respective CI
%disp('Estimate K, Phi and CI')

[K_kept, Phi_kept] = get_pac_params(regressed_pac_kept);

% K_mod Kept
K_mean = mean(K_kept,1);
K_sup = max(K_kept, [],1);
K_inf = min(K_kept,[], 1);

% Phi_mod_kept
[Phi_mean,Phi_sup,Phi_inf] = get_phase_msi(Phi_kept);

% Bin number
Npb = 18;

% Modulogram and Modulation Index
[PAC, pbins, MI_KL_raw, MI_L1_raw] = get_pac_from_kphi(K_mean,Phi_mean,Npb);
[~, ~, MI_KL_raw_sup, MI_L1_raw_sup]   = get_pac_from_kphi(min(K_sup,1),Phi_mean,Npb); %  min(K_sup,1) we sometimes have small rounding issues for K
[~, ~, MI_KL_raw_inf, MI_L1_raw_inf]   = get_pac_from_kphi(min(K_inf,1),Phi_mean,Npb); %  min(K_inf,1) we sometimes have small rounding issues for K

% Gather Estimates (mis = mean - inf - sup)
K_mis     = [K_mean;K_inf;K_sup];
R2_mis    = [R2_cen;R2_inf;R2_sup];
Phi_mis   = [Phi_mean;Phi_inf;Phi_sup];
MI_KL_mis = [MI_KL_raw;MI_KL_raw_inf;MI_KL_raw_sup];
MI_L1_mis = [MI_L1_raw;MI_L1_raw_inf;MI_L1_raw_sup];

%% Filter Latent Pac

if doArFiltering
    disp('AR filtering...')
    
    % Choose a reference to filter pac parameters
    ref_kk = 1;% (-1 -> Raw | = 0 -> mean | = k -> kth sample)
    IC_str = 'BIC'; % = AIC or BIC
    [regressed_pac_filt,ar_coeff] = filter_latent_pac(regressed_pac_kept, ref_kk, IC_str);
    
    
    
    %% Get filtered phase amplitude coupling params K, Phi and respective CI
    %disp('Filtered K, Phi and CI')
    % Get phase amplitude coupling params K, Phi
    [K_kept_filt, Phi_kept_filt] = get_pac_params(regressed_pac_filt);
    
    % K_mod Kept
    K_mean_filt = mean(K_kept_filt,1);
    K_sup_filt  = max(K_kept_filt, [],1);
    K_inf_filt  = min(K_kept_filt,[], 1);
    
    % Phi_mod_kept
    [Phi_mean_filt,Phi_sup_filt,Phi_inf_filt] = get_phase_msi(Phi_kept_filt);
    
    % Modulogram and Modulation Index
    [PAC_filt, pbins, MI_KL_filt, MI_L1_filt] = get_pac_from_kphi(min(K_mean_filt,1),Phi_mean_filt,Npb);
    [~, ~, MI_KL_filt_sup, MI_L1_filt_sup] = get_pac_from_kphi(min(K_sup_filt,1),Phi_mean_filt,Npb);  %  min(K_sup,1) we sometimes have small rounding issues for K
    [~, ~, MI_KL_filt_inf, MI_L1_filt_inf] = get_pac_from_kphi(min(K_inf_filt,1),Phi_mean_filt,Npb);  %  min(K_sup,1) we sometimes have small rounding issues for K
    
    % Gather Estimates (mis = mean - inf - sup)
    K_mis_filt     = [K_mean_filt;K_inf_filt;K_sup_filt];
    Phi_mis_filt   = [Phi_mean_filt;Phi_inf_filt;Phi_sup_filt];
    MI_KL_mis_filt = [MI_KL_filt;MI_KL_filt_sup; MI_KL_filt_inf];
    MI_L1_mis_filt = [MI_L1_filt;MI_L1_filt_sup; MI_L1_filt_inf];
    
end
disp('AR filtered')

%% Gather and save estimates

% Processing Parameters
processed_ssp_pac.alpha_ci=alpha_ci;
processed_ssp_pac.Npb=Npb;
processed_ssp_pac.tVect=tVect;
processed_ssp_pac.name=name;
processed_ssp_pac.Nlinconst=Nlinconst;
processed_ssp_pac.pbins=pbins;

% Sample Numbers
processed_ssp_pac.Kappa  = Kappa;
processed_ssp_pac.Kappa1 = Kappa1;
processed_ssp_pac.Kappa2  = Kappa2;

% Goodness of fit estimate (mis = mean - inf - sup)
processed_ssp_pac.SS_tot=SS_tot;
processed_ssp_pac.R2=R2;
processed_ssp_pac.R2_mis=R2_mis;

% Modulograms
processed_ssp_pac.PAC=PAC;

% PAC raw estimates (mis = mean - inf - sup)
processed_ssp_pac.K_mis=K_mis;
processed_ssp_pac.Phi_mis=Phi_mis;
processed_ssp_pac.MI_KL_mis=MI_KL_mis;
processed_ssp_pac.MI_L1_mis=MI_L1_mis;

if doArFiltering
    
    % Modulograms
    processed_ssp_pac.PAC_filt=PAC_filt;
    
    % Filtering Parameters
    processed_ssp_pac.IC_str=IC_str;
    processed_ssp_pac.ref_kk=ref_kk;
    processed_ssp_pac.ar_coeff=ar_coeff;
    
    % PAC filtered estimates (mis = mean - inf - sup)
    processed_ssp_pac.K_mis_filt = K_mis_filt;
    processed_ssp_pac.Phi_mis_filt = Phi_mis_filt;
    processed_ssp_pac.MI_KL_mis_filt = MI_KL_mis_filt;
    processed_ssp_pac.MI_L1_mis_filt = MI_L1_mis_filt;
    
end


end