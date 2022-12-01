addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../ssp_pac'))

%% Generate 1 slow oscillation, 1 modulated fast during Nwind windows of length N.
Fs=250; % Hz
N=3*Fs; % Window length
Nwind=24;    % Window Number
tt= (1:Nwind*N);

StaPointId = ((1:Nwind)-1)*N +1;
EndPointId = StaPointId+N-1;

f_so = 1;  w_so = f_so *2*pi /Fs; 
f_al = 10; w_al = f_al *2*pi /Fs;
Phi_so  = w_so*tt;


lambda=0.05* 1/Nwind;
t_0= (Nwind/2)*N;


Phi_mod = pi/3;

x_so = 10 * cos(Phi_so);
x_al = 2 * cos(w_al*tt) .* (1 + 0.9 * cos(-Phi_mod+Phi_so) );

sigma_noise= 2;
y_t = x_so +  x_al + normrnd(0,sigma_noise,[1,Nwind*N]);
y_t=y_t';

figure
subplot(2,1,1); hold on
plot(tt/Fs, x_so)
plot(tt/Fs, x_al)
axis tight

subplot(2,1,2)
plot(tt/Fs, y_t,'k')
axis tight

%% Initialize and fit oscillator model
convergenceTolerance=eps;
em_its=100;
doPlot=0;
namecur='test_fit';

init_params=struct();
init_params.f_init     =  [2, 6];
init_params.a_init     =  [0.98, 0.98];
init_params.sigma2_init=  [1,1];
init_params.R          =  1;
init_params.NNharm     =  [1, 1];

modelOsc=ssp_decomp(y_t,Fs,StaPointId,EndPointId,em_its,convergenceTolerance, init_params,namecur,doPlot);

%% Plot a summary of the decomposition
crange = [-20 20];
plot_summary(modelOsc, crange)


%% Use a linear regression and a second state space to estimate and time constrain Phase Amplitude Coupling
Kappa_tot = [10,10]; slowID = 1; fastID = 2;
[processed_ssp_pac,regressed_pac,~, X_t_tot] = ssp_pac_main(modelOsc,slowID,fastID,Kappa_tot);

%% Plots a summary of the pac decomposition
% Non filtered
figure
plot_pac(processed_ssp_pac,0)

% Filtered
figure
plot_pac(processed_ssp_pac,1)























%% Plot 1 Window of interest
WOI = 1;
Fs = modelOsc.Fs;
ta_woi = modelOsc.startPointId(WOI):modelOsc.endPointId(WOI);
yt_woi = y_t(ta_woi);

length_wind = size(X_t_tot,2)-1;
number_wind = size(X_t_tot,3);
kappa_1 =   size(X_t_tot,4);
kappa_tot = size(regressed_pac,2);
kappa_2 = kappa_tot/(kappa_1); 

x_slow = reshape(X_t_tot(2*slowID-1,2:end,WOI,:),length_wind,kappa_1);
x_fast = reshape(X_t_tot(2*fastID-1,2:end,WOI,:),length_wind,kappa_1);

phase_slow = atan2(X_t_tot(2*slowID,2:end,WOI,:), X_t_tot(2*slowID-1,2:end,WOI,:));
phase_slow = reshape(phase_slow,length_wind,kappa_1);

ampli_fast = sqrt(X_t_tot(2*fastID,2:end,WOI,:).^2 +  X_t_tot(2*fastID-1,2:end,WOI,:).^2);
ampli_fast = reshape(ampli_fast,length_wind,kappa_1);

regressed_pac_woi = reshape(regressed_pac(:,:,WOI),3,kappa_tot);

M_t = @(Phi) [ones(length(Phi),1) cos(Phi), sin(Phi)];

regressed_amp = zeros(length_wind,kappa_tot);

for kk1 = 1:kappa_1
    Xreg_t = M_t (phase_slow(:,kk1));
    
    regressed_pac_woi_kk2 = regressed_pac_woi(:,(kk1-1)*kappa_2+1:(kk1)*kappa_2);

    regressed_amp_kk2 = Xreg_t*regressed_pac_woi_kk2;
    regressed_amp(:,(kk1-1)*kappa_2+1:(kk1)*kappa_2) = regressed_amp_kk2;

end

MIS = @(X) [mean(X,2),min(X,[],2),max(X,[],2)];

[x_slow_mean, x_slow_inf, x_slow_sup] = get_mis(MIS(x_slow));
[x_fast_mean, x_fast_inf, x_fast_sup] = get_mis(MIS(x_fast));
[ampli_fast_mean, ampli_fast_inf, ampli_fast_sup] = get_mis(MIS(ampli_fast));
[regressed_amp_mean, regressed_amp_inf, regressed_amp_sup] = get_mis(MIS(regressed_amp));
[~,phase_slow_sup,phase_slow_inf] = get_phase_msi(phase_slow');

figure;
subplot(3,2,1); hold on
plot(ta_woi/Fs,yt_woi, 'k'); box on
axis tight
ylim0=ylim;
title('Signal')
ylabel('[a.u]')

subplot(3,2,2); hold on
ciplot(regressed_amp_inf,regressed_amp_sup,ta_woi/Fs,[0.5 0.5 0.5],1)
plot(ta_woi/Fs,regressed_amp_mean, 'g'); box on
ylim([0 ylim0(2)/2])
title('Regressed Amplitude')
ylabel('[a.u]')

subplot(3,2,3); hold on
ciplot(x_slow_inf,x_slow_sup,ta_woi/Fs,[0.5 0.5 0.5],1)
plot(ta_woi/Fs,x_slow_mean, 'b'); box on
ylim(ylim0)
title('Slow Oscillation')
ylabel('[a.u]')

subplot(3,2,4); hold on
ciplot(x_fast_inf,x_fast_sup,ta_woi/Fs,[0.5 0.5 0.5],1)
plot(ta_woi/Fs,x_fast_mean, 'r'); box on
ylim(ylim0)
title('Fast Oscillation')
ylabel('[a.u]')

subplot(3,2,5); hold on
plot_ci_phase_hs(phase_slow_sup, phase_slow_inf, ta_woi/Fs,[0.5 0.5 0.5])
scatter(ta_woi/Fs,phase_slow(:,1),10, 'b','filled'); box on
ylim([-pi pi])
title('Slow Phase')
xlabel('[sec]')
ylabel('[rad]')

subplot(3,2,6); hold on
ciplot(ampli_fast_inf,ampli_fast_sup,ta_woi/Fs,[0.5 0.5 0.5],1)
plot(ta_woi/Fs,ampli_fast_mean, 'r'); box on
ylim([0 ylim0(2)/2])
title('Fast Amplitude')
xlabel('[sec]')
ylabel('[a.u]')



