addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../ssp_pac'))
addpath(genpath('../usual_pam'))


f_so_transdn=1;
ta_sec2=(1:Fs*6)/Fs;
a_param=1;
K_true_transdn = 0;
Phi_true_transdn = 0;
pac_true_transdn =  get_pac_from_kphi(K_true_transdn,Phi_true_transdn,Npb);

x_t_tmp=cos(2*pi*f_so_transdn*ta_sec2);
x_t_transdn= x_t_tmp+a_param*x_t_tmp.^2;

noiseAmp = 0.1;
y_t_transdn= x_t_transdn + normrnd(0, noiseAmp, [1, length(ta_sec2)] ) ;

figure; hold on
plot(ta_sec2, x_t_tmp)
plot(ta_sec2, y_t_transdn)
legend('x_t','y_t')


%% Standard Processing 
% Bandpass of the solution
slow_freq_transdn = [0.8 1.2]; fast_freq_transdn = [0.9 3.1]; order=100;
[x_filt_slow_transdn,  tail1] =quickbandpass(y_t_transdn,Fs,slow_freq_transdn, order);   %Bandpass at the phase frequency
[x_filt_alpha_transdn, tail2] =quickbandpass(y_t_transdn,Fs,fast_freq_transdn, order);

% Coupling
dPhi= 2*pi/Npb;
pbins=linspace(-pi,pi,Npb);
pbins2=linspace(-pi,pi,Npb-1);

x_al_hilb = hilbert(x_filt_alpha_transdn(1:end-tail1));
x_so_hil  = hilbert(x_filt_slow_transdn(1:end-tail2));

amp_al_hilb_transdn   =abs(x_al_hilb);
phase_so_hilb_transdn = angle(x_so_hil);
pac_standard_transdn=phaseamp(amp_al_hilb_transdn,phase_so_hilb_transdn,pbins);
pac_standard_transdn= pac_standard_transdn/sum(pac_standard_transdn)/dPhi;



%% State-Space Oscillator Decomposition



em_its=400;
harmTest = 8;

f_tot_transdn = [1.5, 2];
a_tot_transdn = [0.99 0.99];
sigma2_tot_transdn = [0.1 0.1];
% Include or not the alpha componemt in the state space oscillator model
modelOsc_harm = cell(1,harmTest);
parfor harmtest_cur = 1:(harmTest)
    disp([' Harm#', num2str(harmtest_cur), '/' , num2str(harmTest)])
    
    NNharm_transdn = [harmtest_cur,1];
    
    % 1 osc + harm
    init_params=struct();
    init_params.f_init     =  f_tot_transdn(1,1:1)
    init_params.a_init     =  a_tot_transdn(1,1:1)
    init_params.sigma2_init=  sigma2_tot_transdn(1,1:1)
    init_params.R        =  4;
    init_params.NNharm      = NNharm_transdn(1,1:1);
    modelOsc_1osc=ssp_decomp(y_t_transdn',Fs,1,length(ta_sec2),em_its,eps, init_params,'',0);
    
    
    
    % 2 osc + harm
    init_params=struct();
    init_params.f_init     =  f_tot_transdn(1,1:2)
    init_params.a_init     =  a_tot_transdn(1,1:2)
    init_params.sigma2_init=  sigma2_tot_transdn(1,1:2)
    init_params.R        =  4;
    init_params.NNharm      =  NNharm_transdn(1,1:2);
    modelOsc_2osc=ssp_decomp(y_t_transdn',Fs,1,length(ta_sec2),em_its,eps, init_params,'',0);
    
    modelOsc_harm{:,harmtest_cur}={modelOsc_1osc;modelOsc_2osc};
    
end




%% Compute Information criterion to Select the best model
lltot = zeros(size(modelOsc_harm{1,1},1),harmTest);
Kparam = zeros(size(modelOsc_harm{1,1},1),harmTest);

for harmtest_cur = 1:(harmTest)
    for j = 1:size(modelOsc_harm{1,1},1)   
        lltot(j,harmtest_cur) = modelOsc_harm{1,harmtest_cur}{j}.res{1,1}.ll - ( 0.5*length(y_t_transdn) *log(2*pi)); 
        Kparam_tmp =   sum(3+2*(modelOsc_harm{1,harmtest_cur}{j}.init_params.NNharm-1)) +1;
        Kparam(j,harmtest_cur) = Kparam_tmp;
    end  
end





AIC_transdn = -2*lltot + 2 * Kparam;
BIC_transdn = -2*lltot + 2 * Kparam * log(length(y_t_transdn));

dAIC_transdn = AIC_transdn- min(min(AIC_transdn));
dBIC_transdn = BIC_transdn- min(min(BIC_transdn));

pdAIC = exp(-0.5*dAIC_transdn)./sum(sum(exp(-0.5*dAIC_transdn)));
pdBIC = exp(-0.5*dBIC_transdn)./sum(sum(exp(-0.5*dBIC_transdn)));

IC = dAIC_transdn;
[mIC, locMax_1_tmp] = min(dAIC_transdn);
[~, locMax_2]   = min(mIC);

locX = locMax_1_tmp(locMax_2);
locY = locMax_2;

best_model_transdn = modelOsc_harm{1,locY}{locX};
NNharm_transdn     = best_model_transdn.init_params.NNharm;



x_t_n_transdn   = best_model_transdn.res{1,1}.x_t_n;
ta_sec2 = best_model_transdn.res{1,1}.ta;

x_t_n_slow_transdn =  sum(x_t_n_transdn(1+2*((1:NNharm_transdn(1))-1),2:end),1);
x_t_n_fast_transdn = zeros(size(x_t_n_slow_transdn));
K_ssp_transdn = 0;
Phi_ssp_transdn = 0;

if length(NNharm_transdn)>1
    x_t_n_fast_transdn = x_t_n_transdn(1+2*sum(NNharm_transdn(1)),:);
    Kappa_tot = [10;10];
    slowID= 1; fastID=2;
    Nlag = Inf; NlinConst = 100;
    [regressed_pac, SS_tot, X_t_tot] = ssp_pac (best_model_transdn, Kappa_tot, slowID, fastID, Nlag,NlinConst);
    [K_ssp_transdn, Phi_ssp_transdn] = get_pac_params(regressed_pac);
    
end

pac_ssp_transdn =  get_pac_from_kphi(K_ssp_transdn,Phi_ssp_transdn,Npb);




%%


figure;
subplot(3,6,13); hold on
plot(ta_sec2, y_t_transdn-mean(y_t_transdn),'k')
plot(ta_sec2, x_t_transdn-mean(x_t_transdn),'b')
box on
axis tight
ylim0 = ylim;

subplot(3,6,14); hold on
plot(ta_sec2, x_filt_slow_transdn,'b')
plot(ta_sec2, x_filt_alpha_transdn,'r')
ylim(ylim0)
box on


subplot(3,6,15); hold on
plot(pbins, pac_true_transdn,'k','linewidth',2)
plot(pbins2, pac_standard_transdn,'g','linewidth',2)
axis tight
ylim1 = ylim;
box on


subplot(3,6,16); hold on
plot(1:harmTest, dAIC_transdn(1,:),'color',[1 0.6 0], 'linewidth', 2)
plot(1:harmTest, dBIC_transdn(1,:),'-.','color',[1 0.6 0], 'linewidth', 2)

plot(1:harmTest, dAIC_transdn(2,:),'color',[0.5 0.5 0.5], 'linewidth', 2)
plot(1:harmTest, dBIC_transdn(2,:),'-.','color',[0.5 0.5 0.5], 'linewidth', 2)
axis tight;box on


subplot(3,6,17); hold on
plot(ta_sec2,x_t_n_slow_transdn,'b')
if length(NNharm_transdn)>1
   plot(ta_sec2,x_t_n_fast_transdn,'r') 
end
ylim(ylim0);  box on


subplot(3,6,17); hold on
plot(ta_sec2,x_t_n_slow_transdn,'b')
if length(NNharm_transdn)>1
   plot(ta_sec2,x_t_n_fast_transdn,'r') 
end
ylim(ylim0);  box on


subplot(3,6,18); hold on
plot(pbins, pac_true_transdn,'k','linewidth',2) 
if length(NNharm_transdn)>1

plot(pbins, pac_ssp_transdn,'g','linewidth',2)
else
   
   plot(pbins, pac_ssp_transdn,':g','linewidth',2)
end

ylim(ylim1);  box on


