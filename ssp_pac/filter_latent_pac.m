function [latent_pac_filt,ar_coeff] = filter_latent_pac(latent_pac, ref_kk, IC_str)
%% Inputs : Latent pac ~ dim_pac x Nsamp    x Nwind
%                       2 or 3  x #Samples x number of windows
%           ref_kk : reference used for the filter
%                    =-1 -> Choose raw estimate
%                    =0  -> Choose mean estimate
%                    =k  -> Choose a given sample
%           IC_str : String 'AIC' or 'BIC' chose AR order method
%
% Outputs : latent_pac_filt ~ dim_pac x Nsamp    x Nwind
%               latent_pac filtered using an AR model which coefficient are
%               estimate with BURG algorithm on latent_pac(:,ref_kk,:)
%

dim_pac = size(latent_pac,1);
Nsamp   = size(latent_pac,2);
Nwind   = size(latent_pac,3);

if nargin<3
    IC_str = 'BIC';
end

if nargin<2
    ref_kk = 0;
end

if ref_kk ==-1
    % Choose raw estimate
    latent_pac_ref = reshape(latent_pac(:,1,:),dim_pac,Nwind);
elseif ref_kk ==0
    % Choose mean estimate
    latent_pac_ref = reshape(mean(latent_pac(:,:,:),2),dim_pac,Nwind);
else
    % Choose a given sample
    latent_pac_ref = reshape(latent_pac(:,ref_kk,:),dim_pac,Nwind);
end


ar_order_max = max(1,min(10,floor(Nwind/5)));
ll_tot= zeros(1,ar_order_max);


for ARp = 1:ar_order_max 
    %disp(['AR order test : ', num2str(ARp), '/', num2str(ar_order_max)])
    [~,~,~,~,mll] = mvyule_ar(latent_pac_ref,ARp);
    ll_tot(1,ARp) = -mll;
end

% Eeach AR param is a dim_pac x dim_pac matrix 
Nparam = (1:ar_order_max)*dim_pac*dim_pac;

AIC = 2*Nparam-2*ll_tot;
BIC = Nparam*log(Nwind)- 2*ll_tot;

% Use AIC or BIC to select optimal filter
if strcmp(IC_str,'AIC')
    [~,ARp_hat] = min(AIC,[],2);
elseif strcmp(IC_str,'BIC')
    [~,ARp_hat] = min(BIC,[],2);
else 
    error('Incorect String for Information criterion')
end

[Phi_hat,M,Q_hat,R_hat] = mvyule_ar(latent_pac_ref,ARp_hat);
%mu = zeros(dim_pac*ARp_hat,1);
mu = repmat([1;0.01;0], ARp_hat,1);
Sigma = 1*eye(dim_pac*ARp_hat,dim_pac*ARp_hat);
ar_coeff = Phi_hat(1:dim_pac,:);


gamma_kk = kalman_filter_mv(latent_pac_ref,M,Phi_hat,Q_hat,R_hat,mu,Sigma);
latent_pac_ref_filt = gamma_kk(1:dim_pac,2:end);


% Filter each sample (with same coefficients ?)
latent_pac_filt = zeros(dim_pac,Nsamp,Nwind);
parfor kk = 1:size(latent_pac,2)   
    %disp(['Filtering, sample # ', num2str(kk), '/' , num2str(Nsamp)])
    gamma_kk = kalman_filter_mv(squeeze(latent_pac(:,kk,:)),M,Phi_hat,Q_hat,R_hat,mu,Sigma); % Modify the squeeze to reshape
    latent_pac_filt(:,kk,:)=gamma_kk(1:dim_pac,2:end);
end





end








