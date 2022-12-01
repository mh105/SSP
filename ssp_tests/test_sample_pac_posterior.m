%% Test the second step of (d)SSP: linear regression of modulation parameters and posterior sampling 

addpath(genpath('./../ssp_pac'))

%% Generate PAC coupled Signal
Fs = 250;
Tlength_sec = 6;
n = Fs*Tlength_sec;
tt = 1:n;
f_so = 0.5;

% Modulation Amplitude
K_mod =0.9;

% Modulation Phase
Phi_mod = pi/3;

% Fast Amplitude nominal
A_0 =1;

noise_level = A_0* 0.05;
obs_nose =  10*normrnd(0,noise_level, 1,n);
phase_noise = 2*pi *  normrnd(0,noise_level, 1,n);

% Slow Phase
phi_t_so_0 = 2*pi*f_so/Fs*tt + phase_noise;

% Fast amplitude
A_al_t_norm =A_0*(1+ K_mod*cos(Phi_mod +phi_t_so_0 )) +obs_nose;

% Plot
figure; 
subplot(2,1,1) 
scatter(tt/Fs, wrapToPi(phi_t_so_0),10,'b', 'filled')
title('Slow Phase')
xlabel('[sec]')
ylim([-pi pi])

subplot(2,1,2) 
plot(tt/Fs, A_al_t_norm,'r')
title('Fast Amplitude')
xlabel('[sec]')

%% Define Regressors, prior and posterior parameters
y_t = A_al_t_norm';
X_t = [ones(length(phi_t_so_0),1)  cos(phi_t_so_0') sin(phi_t_so_0')];


% Prior and posterior Parameters : a, nu, b, V
NlinConst = 100;
dim_y = n;
dim_x1 = NlinConst;
dim_x2 = 2; %alpha = dim_x2;

theta_r = linspace(-pi,pi, NlinConst)';
f =  zeros(NlinConst,1); 
W = [ones(NlinConst,1) -cos(theta_r) -sin(theta_r)];

A_mean_t = mean(y_t);
a_prior  = [A_mean_t;0; 0];
nu_prior = 3;
b_prior  = 1;
V_prior = (1/A_mean_t) * diag([3,12,12])* b_prior * nu_prior/(nu_prior-2);

% Posterior and posterior Parameters : a, nu, b, V
[a_post,nu_post, V_post,b_post] = get_post_param(y_t,X_t,a_prior,nu_prior,V_prior,b_prior);

% Check posterior estimates
A_0_post = a_post(1,1);
K_post = sqrt(a_post(2,1).^2+a_post(3,1).^2)/a_post(1,1);
Phi_post = -atan2(a_post(3,1),a_post(2,1));
disp(['Esti : A_0 = ',num2str(A_0_post), ' K^{mod} = ', num2str(K_post), ' phi^{mod} = ', num2str(Phi_post)])
disp(['True : A_0 = ',num2str(A_0), ' K^{mod} = ', num2str(K_mod), ' phi^{mod} = ', num2str(Phi_mod)])
 
%% Draw samples and estimate prior/posterior pdf

T_so_hat = min(10,1/mean(diff(phi_t_so_0)*Fs/(2*pi)));
T_wind = length(phi_t_so_0)/Fs;
Cmax = 1/(1-T_so_hat/(T_wind*pi));
assert(Cmax>=1)
A_0_max = A_mean_t * Cmax;

dbeta=0.025 * A_0_max ; % Sample resolution
betamLim = A_0_max;

res = 100;  % Resolution multiplicator
Nsample = 2000; % Number of Samples

disp(['Pdf Pior'])
[pdf_prior,beta_0_grid,beta_1_grid,beta_2_grid,beta_0,beta_1,beta_2] = mvt3pdf_hs(dbeta,betamLim ,a_prior, nu_prior,V_prior, b_prior);
[pdf_prior_const] = mvt3_const_pdf(dbeta,betamLim ,a_prior, nu_prior, V_prior, b_prior, W,f);

disp(['Sample Post'])
[pdf_post] = mvt3pdf_hs(dbeta,betamLim ,a_post, nu_post,V_post, b_post);
[pdf_post_const] = mvt3_const_pdf(dbeta,betamLim ,a_post, nu_post, V_post, b_post, W,f);

disp(['Sample Unconstrained'])
beta_samp_prior = mvt3rnd_hs(a_prior,b_prior*inv(V_prior),nu_prior, Nsample);
beta_samp_post  = mvt3rnd_hs(a_post, b_post *inv(V_post) ,nu_post , Nsample);

disp(['Sample Constrained'])
[beta_samp_prior_const] = sample_from_pdf3(beta_0,beta_1,beta_2,pdf_prior_const,Nsample,res);
[beta_samp_post_const]  = sample_from_pdf3(beta_0,beta_1,beta_2,pdf_post_const,Nsample,res);

%% Deduce the PAC parameters K_mod and Phi_mod
disp(['Pac Params'])
K_hat_prior =sqrt( beta_samp_prior(2,:).^2 + beta_samp_prior(3,:).^2)./beta_samp_prior(1,:);
K_hat_post  =sqrt( beta_samp_post(2,:).^2  + beta_samp_post(3,:).^2) ./beta_samp_post(1,:);

K_hat_prior_const =sqrt( beta_samp_prior_const(2,:).^2 + beta_samp_prior_const(3,:).^2)./beta_samp_prior_const(1,:);
K_hat_post_const  =sqrt( beta_samp_post_const(2,:).^2 + beta_samp_post_const(3,:).^2)./beta_samp_post_const(1,:);

Phi_hat_prior = -atan2( beta_samp_prior(3,:),  beta_samp_prior(2,:));
Phi_hat_post  = -atan2( beta_samp_post(3,:),  beta_samp_post(2,:));

Phi_hat_prior_const = -atan2( beta_samp_prior_const(3,:),  beta_samp_prior_const(2,:));
Phi_hat_post_const  = -atan2( beta_samp_post_const(3,:),  beta_samp_post_const(2,:));


%% Plot: pdf, drawn samples, K_mod and Phi_mod histograms 

beta_0_test =A_mean_t;
[~,beta_0_id] = min(abs(beta_0-beta_0_test));

figure; 
subplot(4,4,1)
surf(squeeze(beta_1_grid(beta_0_id,:,:)),squeeze(beta_2_grid(beta_0_id,:,:)), squeeze(pdf_prior(beta_0_id,:,:)),'linestyle', 'none')
view([0 1 0])
xlim([-betamLim betamLim])
ylim([-betamLim betamLim])
title('Prior')
xlabel('\beta_1')
ylabel('\beta_2')
zlabel(['pdf(.| \beta_0 = ', num2str(beta_0_test) ,')'],'fontweight','bold')

subplot(4,4,2)
surf(squeeze(beta_1_grid(beta_0_id,:,:)),squeeze(beta_2_grid(beta_0_id,:,:)), squeeze(pdf_prior_const(beta_0_id,:,:)),'linestyle', 'none')
view([0 1 0])
xlim([-betamLim betamLim])
ylim([-betamLim betamLim])
title('Prior Constrained')

subplot(4,4,3)
surf(squeeze(beta_1_grid(beta_0_id,:,:)),squeeze(beta_2_grid(beta_0_id,:,:)), squeeze(pdf_post(beta_0_id,:,:)),'linestyle', 'none')
view([0 1 0])
xlim([-betamLim betamLim])
ylim([-betamLim betamLim])
title('Posterior')

subplot(4,4,4)
surf(squeeze(beta_1_grid(beta_0_id,:,:)),squeeze(beta_2_grid(beta_0_id,:,:)), squeeze(pdf_post_const(beta_0_id,:,:)),'linestyle', 'none')
view([0 1 0])
xlim([-betamLim betamLim])
ylim([-betamLim betamLim])
title('Posterior Constrained')

beta_1_renorm = linspace(-1,1,length(beta_1));
beta_2_renorm = linspace(-1,1,length(beta_2));

subplot(4,4,5)
histogram2(beta_samp_prior(2,:)./beta_samp_prior(1,:),beta_samp_prior(3,:)./beta_samp_prior(1,:),beta_1_renorm,beta_2_renorm,'Normalization', 'probability', 'facecolor', 'flat')
view(2)
xlabel('\beta_1 / \beta_0')
ylabel('\beta_2 / \beta_0')
title([num2str(Nsample), ' samples'])
xlim([-1 1])
ylim([-1 1])

subplot(4,4,6)
histogram2(beta_samp_prior_const(2,:)./beta_samp_prior_const(1,:),beta_samp_prior_const(3,:)./beta_samp_prior_const(1,:),beta_1_renorm,beta_2_renorm,'Normalization', 'probability', 'facecolor', 'flat')
view(2)
xlabel('\beta_1 / \beta_0')
ylabel('\beta_2 / \beta_0')
xlim([-1 1])
ylim([-1 1])



subplot(4,4,7)
histogram2(beta_samp_post(2,:)./beta_samp_post(1,:),beta_samp_post(3,:)./beta_samp_post(1,:),beta_1_renorm,beta_2_renorm,'Normalization', 'probability', 'facecolor', 'flat')
view(2)
xlabel('\beta_1 / \beta_0')
ylabel('\beta_2 / \beta_0')
xlim([-1 1])
ylim([-1 1])


subplot(4,4,8)
histogram2(beta_samp_post_const(2,:)./beta_samp_post_const(1,:),beta_samp_post_const(3,:)./beta_samp_post_const(1,:),beta_1_renorm,beta_2_renorm,'Normalization', 'probability', 'facecolor', 'flat')
view(2)
xlabel('\beta_1 / \beta_0')
ylabel('\beta_2 / \beta_0')
xlim([-1 1])
ylim([-1 1])


subplot(4,4,9)
histogram(K_hat_prior,'Normalization', 'probability');
ylabel(['K^{mod}'],'fontweight','bold')
title(['K^{mod}_{true} = ' , num2str(K_mod)])
xlim([0 1])


subplot(4,4,13)
histogram(Phi_hat_prior,'Normalization', 'probability');
ylabel(['\phi^{mod}'],'fontweight','bold')
title(['\phi^{mod}_{true} = ' , num2str(Phi_mod)])
xlim([-pi pi])


subplot(4,4,10)
histogram(K_hat_prior_const,'Normalization', 'probability');
xlim([0 1])


subplot(4,4,14)
histogram(Phi_hat_prior_const,'Normalization', 'probability');
xlim([-pi pi])

subplot(4,4,11)
histogram(K_hat_post,'Normalization', 'probability');
title(['K^{mod}_{mean} = ' , num2str(mean(K_hat_post))])
xlim([0 1])


subplot(4,4,15)
histogram(Phi_hat_post,'Normalization', 'probability');
title(['\phi^{mod}_{mean} = ' , num2str(mean(Phi_hat_post))])
xlim([-pi pi])

subplot(4,4,12)
histogram(K_hat_post_const,'Normalization', 'probability');
title(['K^{mod}_{mean} = ' , num2str(mean(K_hat_post_const))])
xlim([0 1])


subplot(4,4,16)
histogram(Phi_hat_post_const,'Normalization', 'probability');
title(['\phi^{mod}_{mean} = ' , num2str(mean(Phi_hat_post_const))])
xlim([-pi pi])





