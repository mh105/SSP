
% Generate Data
addpath(genpath('./../'))

Fs = 250;        % SamplingFrequency
Nwindow=5000;    % Number of Window
lengthWindSec=6; % Window length ins second
lengthTot=Nwindow*lengthWindSec*Fs;
t_tot=1:lengthTot;
starind=t_tot(1:lengthWindSec*Fs:end-lengthWindSec*Fs+1)';
stopind=t_tot(lengthWindSec*Fs:lengthWindSec*Fs:end)';

Nstep=6; % Noise level number
Tstep=linspace(1,lengthTot,Nstep+1)';
StepFunction=sum(heaviside(t_tot-Tstep(1:end-1)),1)/Nstep;


K_th  =0.9*ones(1,lengthTot);
Phi_th=(pi/3)*ones(1,lengthTot);

w_so_th = (2*pi/Fs)*linspace(0.5,0.5 ,lengthTot);
w_al_th = (2*pi/Fs)*floor(10*linspace(12,12 ,lengthTot))/10;

SNR=40*StepFunction;

phi_so_th=   w_so_th.*t_tot;
phi_so_th_n= w_so_th.*t_tot+w_so_th.*normrnd(0,SNR);%+noise

phi_al_th= w_al_th.*t_tot;
phi_al_th_n= w_al_th.*t_tot+w_so_th.*normrnd(0,SNR);%+noise

x_so_th=10*cos(phi_so_th_n);
x_al_th=5*cos(phi_al_th_n);
x_al_th_mod=x_al_th.*(1+K_th.*cos(phi_so_th_n+Phi_th));

y_tot_th= x_so_th+ x_al_th_mod +normrnd(0,SNR);%+noise obs

%figure; plot(t_tot/Fs, y_tot_th)

%% Plot Samples from each noise level

StepIndexes = find(diff(StepFunction));
StepSize    = floor(mean(diff(StepIndexes)));
timeWind    = 4*Fs;
SampleStart=StepIndexes+timeWind;
SampleEnd=SampleStart+timeWind;
tickVal = unique(SNR);
figure
ylimTmp=[-120 120];

for jj=1:length(SampleStart)
    subplot(3,length(SampleStart),jj ); hold on
    ta=t_tot(SampleStart(jj):SampleEnd(jj))/Fs; ta=ta-ta(1);
    Ycurr=y_tot_th(SampleStart(jj):SampleEnd(jj));
    x_so_cur=x_so_th(SampleStart(jj):SampleEnd(jj));
    x_al_cur=x_al_th_mod(SampleStart(jj):SampleEnd(jj));
    plot(ta,  Ycurr, 'color', 'k')
    scalemuVolt=20; scaleTime=1;
    axis off;title(['\sigma^{2}_N = ',num2str(floor(tickVal(jj+1)))]); ylim([ylimTmp]);
    if jj==1
    plot([ta(1); ta(1)], [min(Ycurr); min(Ycurr)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(Ycurr); min(Ycurr)], '-k', 'LineWidth', 3)
    text(ta(1),min(Ycurr)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end
    
    subplot(3,length(SampleStart),jj+ length(SampleStart)*2); hold on
    plot(ta,  x_so_cur, 'color', 'b')
    plot(ta,  x_al_cur, 'color', 'R')
    ylimTmp2=ylimTmp/4;
    ylim([ylimTmp2]);
    if jj==1
    plot([ta(1); ta(1)], [min(x_so_cur); min(x_so_cur)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so_cur); min(x_so_cur)], '-k', 'LineWidth', 2)
    text(ta(1),min(x_so_cur)-(ylimTmp2(2)-ylimTmp2(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end
    axis off;
  
    subplot(3,length(SampleStart),jj+ length(SampleStart)); hold on
    plot(ta,  Ycurr-x_al_cur-x_so_cur, 'color', [0.5 0.5 0.5])
    ylim([ylimTmp]);

    axis off;
    if jj==1
    plot([ta(1); ta(1)], [min(Ycurr); min(Ycurr)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(Ycurr); min(Ycurr)], '-k', 'LineWidth', 3)
    text(ta(1),min(Ycurr)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end

end
set(gcf, 'Position', [ 0.6667    0.6667  853.3333  469.3333])
set(gcf,'PaperOrientation','landscape');
%print(gcf, '-dpdf', ['/Users/apple/Desktop/figure4supalt'])



%% Osc decomposition
dpath = '/media/hdhs/data/snr_dss_vs_standard/';
save([dpath, 'snr_dss_vs_standard_simulated_signal'])

em_its=500;
init_params=struct();
init_params.f_init     =  [5, 9];
init_params.a_init     =  [0.99, 0.99];
init_params.sigma2_init=  [1, 1];
init_params.R          =  100;
init_params.NNharm      =  [1,1];

modelOsc=dss_decomp(y_tot_th',Fs,starind,stopind,em_its,eps, init_params,'',0,0);
%save([dpath, 'snr_dss_vs_standard_modelOsc'], 'modelOsc', '-v7.3')

%% PAC estimation
slowID = 1; fastID = 2; Kappa = 0;
[K_mod_t_0, Phi_mod_t_0] = dss_pac(modelOsc, Kappa, slowID, fastID);
%save([dpath, 'snr_dss_vs_standard_PAC'], 'K_mod_t_0', 'Phi_mod_t_0', '-v7.3')






