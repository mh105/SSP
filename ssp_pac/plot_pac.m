%function plot_pac(tVect,pbins, PAC ,R2_mis, K_mis ,Phi_mis , MI_KL_mis,MI_L1_mis,Kappa1,Kappa2,alpha_ci)
function plot_pac(processed_ssp_pac, filtered)

%% Plots a summary of the parametric modulation estimation
%
% INPUTS: - tVect time vector, length:Nwind
%         - PAC : parametric phase amplitude coupling Nwind x Npb
%         - mis: mean, inferior and superior bound
%         - R2 variance explained
%         - Phi, K, modulation phase and strength
%         - MI modulation index
%         - L1: L1 norm
%         - KL: Kullback-Leibler
%         - Kappa1, Kappa2 : number samples drawn from the posterior
%         distributions
%         - alpha_ci, fraction used to construct Credible intervals
%         - filtered, boolean, plot the raw or filtered data

tVect   = processed_ssp_pac.tVect;
pbins   = processed_ssp_pac.pbins;
R2_mis  = processed_ssp_pac.R2_mis;
Kappa1 = processed_ssp_pac.Kappa1;
Kappa2 = processed_ssp_pac.Kappa2;
alpha_ci = processed_ssp_pac.alpha_ci;

if nargin<2
    filtered = 0;
end

if filtered
    PAC     = processed_ssp_pac.PAC_filt;
    K_mis   = processed_ssp_pac.K_mis_filt;
    Phi_mis = processed_ssp_pac.Phi_mis_filt;
    MI_KL_mis = processed_ssp_pac.MI_KL_mis_filt;
    MI_L1_mis = processed_ssp_pac.MI_L1_mis_filt;
else
    PAC     = processed_ssp_pac.PAC;
    K_mis   = processed_ssp_pac.K_mis;
    Phi_mis = processed_ssp_pac.Phi_mis;
    MI_KL_mis = processed_ssp_pac.MI_KL_mis;
    MI_L1_mis = processed_ssp_pac.MI_L1_mis;
end


R2_mean = R2_mis(1,:);
if size(R2_mis,1) == 3
    R2_inf = R2_mis(2,:);
    R2_sup = R2_mis(3,:);
end

K_mean = K_mis(1,:);
if size(K_mis,1) == 3
    K_inf = K_mis(2,:);
    K_sup = K_mis(3,:);
end

Phi_mean = Phi_mis(1,:);
if size(R2_mis,1) == 3
    Phi_inf = Phi_mis(2,:);
    Phi_sup = Phi_mis(3,:);
end

MI_KL_mean = MI_KL_mis(1,:);
if size(MI_KL_mis,1) == 3
    MI_KL_inf = MI_KL_mis(2,:);
    MI_KL_sup = MI_KL_mis(3,:);
end

MI_L1_mean = MI_L1_mis(1,:);
if size(MI_KL_mis,1) == 3
    MI_L1_inf = MI_L1_mis(2,:);
    MI_L1_sup = MI_L1_mis(3,:);
end


subplot(3,2,1 ); hold on
ciplot(R2_inf,R2_sup,tVect,[0.5 0.5 0.5])
plot(tVect,R2_mean)
box on
xlim([tVect(1) tVect(end)])
xlabel('[min]')
ylim([0 1])
title('R^2')

subplot(3,2,3); hold on
cmin = 0; cmax = 2/(2*pi);
imagesc(tVect ,pbins,PAC);
col = redbluemap; colormap(col);
colorbar('location', 'west')
set(gca, 'clim', [cmin cmax]);
axis tight
xlim([tVect(1) tVect(end)])
ylabel('(rad)')
xlabel('[min]')
title('PAC')
box on

subplot(3,2,2);hold on
ciplot(K_inf,K_sup,tVect,[0.5 0.5 0.5])
l1=plot(tVect,K_sup,'Color',[0.5 0.5 0.5]);
l2=plot(tVect,K_mean,'g');
ylim([0 1])
xlim([tVect(1) tVect(end)])
box on
xlabel('[min]')
legend([l1,l2], {[num2str(floor(100*(1-alpha_ci))), '% CI'],'Mean'})
title(['K^{mod}_t. ', num2str(Kappa1),'\times', num2str(Kappa2), 'samples'])

subplot(3,2,4);hold on
if (Kappa1+Kappa2)>0
    plot_ci_phase_hs(Phi_sup,Phi_inf, tVect,[0.5 0.5 0.5])
end
scatter(tVect,Phi_mean,10,'g','filled')
ylim([-pi, pi])
xlim([tVect(1) tVect(end)])
box on
xlabel('[min]')
title(['\phi^{mod}_t. ', num2str(Kappa1),'\times', num2str(Kappa2), 'samples'])

subplot(3,2,5 ); hold on
ciplot(MI_KL_inf,MI_KL_sup,tVect,[0.5 0.5 0.5])
plot(tVect,MI_KL_mean,'k')
box on
xlim([tVect(1) tVect(end)])
xlabel('[min]')
title('MI_{KL}')

subplot(3,2,6); hold on
ciplot(MI_L1_inf,MI_L1_sup,tVect,[0.5 0.5 0.5])
plot(tVect,MI_L1_mean,'k')
box on
xlim([tVect(1) tVect(end)])
xlabel('[min]')
title('MI_{L_1}')



end