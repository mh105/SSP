addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../ssp_pac'))
addpath(genpath('../usual_pam'))

%% Create modulated sin
f_so=2.5;
f_al=10;
Fs=250;
ta=(1:Fs*6)/Fs;
K_mod_th=0.7;
phi_0=pi/3;
x_t=cos(2*pi*f_so*ta);
PhisoNoise=0*2*pi*0.2*rand(1, length(ta));

obsNoise=2*rand(1, length(ta));

x_so_th=4*cos(2*pi*f_so*ta);
x_al_th=cos(2*pi*f_al*ta);

y_t_tot_nonnoisy=x_so_th+x_al_th.*(1+a.*cos(2*pi*f_so*ta+phi_0));
y_t_tot=x_so_th+x_al_th.*(1+a.*(1+0.3*rand(1, length(ta))).*cos(2*pi*f_so*ta+phi_0)+PhisoNoise)+ obsNoise;

figure; hold on
plot(ta, x_t)
plot(ta, y_t_tot)
legend('x_t','y_t')

%% "Badly" Bandpass of the solution
phasefreq = [2 4]; ampfreq = [9 11]; order=100;
[x_slow, tail1] =quickbandpass(y_t_tot,Fs,phasefreq, order);   %Bandpass at the phase frequency
[x_alpha, tail2]=quickbandpass  (y_t_tot,Fs,ampfreq, order);

% phasefreq = [2 4]; ampfreq = [5 15]; order=100;
% [x_slow_good, tail1] =quickbandpass(y_t_tot,Fs,phasefreq, order);   %Bandpass at the phase frequency
% [x_alpha_good, tail2]=quickbandpass  (y_t_tot,Fs,ampfreq, order);

%% Use Oscillation Decomposition
em_its=350;
init_params=struct();
init_params.f_init     =  [4, 10];
init_params.a_init     =  [0.99, 0.99];
init_params.sigma2_init=  [1, 5];
init_params.R          =  5;
init_params.NNharm      =  [1,1];
modelOsc=ssp_decomp(y_t_tot',Fs,1,length(ta),em_its,eps, init_params,'',1);

Phi=modelOsc.res{1,1}.model_prams.Phi;
Q=modelOsc.res{1,1}.model_prams.Q;
q_so=Q(1,1);
q_al=Q(3,3);
[a_so,w_so]=get_rotmat_pam(Phi(1:2,1:2));
[a_al,w_al]=get_rotmat_pam(Phi(3:4,3:4));

f_so=w_so*Fs/(2*pi);
f_al=w_al*Fs/(2*pi);

x_so=modelOsc.res{1,1}.x_t_n(1,2:end);
x_al=modelOsc.res{1,1}.x_t_n(3,2:end);
phi_so=atan2(modelOsc.res{1,1}.x_t_n(2,2:end),modelOsc.res{1,1}.x_t_n(1,2:end));
phi_al=atan2(modelOsc.res{1,1}.x_t_n(4,2:end),modelOsc.res{1,1}.x_t_n(3,2:end));
A_al=x_al./cos(phi_al);
tafit = modelOsc.res{1,1}.ta;




%% Parametric PAC
Kappa_tot = [200;200];
slowID= 1;
fastID= 2;
Nlag = Inf; NlinConst = 100;

[regressed_pac, SS_tot, X_t_tot] = ssp_pac (modelOsc, Kappa_tot, slowID, fastID, Nlag,NlinConst);
[K, Phi, A_0] = get_pac_params(regressed_pac);

kappa_tot = size(K,1);
Npb = 18;
pac_param = zeros(Npb,kappa_tot);
pbins = linspace(-pi,pi,Npb);

for kk=1:kappa_tot
    [pac, pbins, mi_kl, mi_l1] = get_pac_from_kphi(K(kk,1),Phi(kk,1),Npb);
    pac_param(:,kk)=pac;
end

pac_param_inf  = min(pac_param,[],2);
pac_param_sup  = max(pac_param,[],2);
pac_param_mean = mean(pac_param,2);


%% Standard PAC
pbins2=linspace(-pi,pi,Npb-1);

x_al_hilb= hilbert(x_alpha(1:end-tail1));
env_Al   =abs(x_al_hilb);
x_so_hil=hilbert(x_slow(1:end-tail2));
phase_so= angle(x_so_hil);


dPhi = 2*pi/Npb;
pa_usual=phaseamp(env_Al,phase_so,pbins);
pa_usual=pa_usual/sum(pa_usual)/dPhi;


% Check Spectral Content
specparams.Fs = Fs;
specparams.tapers = [4 2];

[Sso_mt_osc,f1]=mtspectrumc(x_so,specparams);
[SAl_mt_osc,f1]=mtspectrumc(x_al,specparams);
[Sso_mt_filt,f2]=mtspectrumc(x_so_hil,specparams);
[SAl_mt_filt,f2]=mtspectrumc(x_al_hilb,specparams);
[SyMt,fyt]=mtspectrumc(y_t_tot,specparams);
[SyMt_nonnoisy ,fyt]=mtspectrumc(y_t_tot_nonnoisy,specparams);
[SsoMt_nonnoisy,fyt]=mtspectrumc(x_so_th,specparams);
[Salt_nonnoisy ,fyt]=mtspectrumc(x_al_th,specparams);

% [Sso_mt_osc,f1]=mtspectrumc(x_so,specparams);
% [SAl_mt_osc,f1]=mtspectrumc(x_al,specparams);
% [Sso_mt_filt,f2]=mtspectrumc(x_so_hil,specparams);
% [SAl_mt_filt,f2]=mtspectrumc(x_al_hilb,specparams);
% [SyMt,fyt]=mtspectrumc(y_t_tot,specparams);
% [SyMt_nonnoisy ,fyt]=mtspectrumc(y_t_tot_nonnoisy,specparams);
% [SsoMt_nonnoisy,fyt]=mtspectrumc(x_so_th,specparams);
% [Salt_nonnoisy ,fyt]=mtspectrumc(x_al_th,specparams);


%% Plots 
tafit=modelOsc.res{1,1}.ta;
tspan=4*Fs:5*Fs;
scaleTime=1;
scalemuVolt=2;
ylim1 = -50;
ylim2 = 20;
xlim1 = 0;
xlim2 = 20;
figure(1);

% Filtered Signal
p1=subplot(2,3,1); hold on
plot(ta(tspan), y_t_tot(tspan), 'color', 'k')
plot(tafit(tspan),x_slow(tspan), 'linewidth', 1.3, 'color', 'b')
plot(tafit(tspan),x_alpha(tspan), 'linewidth', 1.3,'color', 'r')
title('Signal')
set(gca,'xtick',[],'ytick',[])
box on;axis tight
ylimTmp=ylim;ymin1=ylim;
plot([ta(tspan(1)); ta(tspan(1))], [ymin1(1); ymin1(1)+scalemuVolt], '-k',  [ta(tspan(1)); ta(tspan(1+scaleTime*Fs))], [ymin1(1);ymin1(1)], '-k', 'LineWidth', 2)
text(ta(tspan(1)),ymin1(1)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
ylabel('Filtered','FontWeight','bold', 'Fontsize', 11)

% Decomposed Signal
p2=subplot(2,3,4); hold on
plot(ta(tspan), y_t_tot(tspan), 'color', 'k')
plot(tafit(tspan),x_so(tspan), 'linewidth', 1.3, 'color', 'b')
plot(tafit(tspan),x_al(tspan), 'linewidth', 1.3, 'color', 'r')
linkaxes([p1 p2])
axis tight
set(gca,'xtick',[],'ytick',[])
box on;axis tight
ylimTmp=ylim;ymin1=ylim;
plot([ta(tspan(1)); ta(tspan(1))], [ymin1(1); ymin1(1)+scalemuVolt], '-k',  [ta(tspan(1)); ta(tspan(1+scaleTime*Fs))], [ymin1(1);ymin1(1)], '-k', 'LineWidth', 2)
text(ta(tspan(1)),ymin1(1)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
ylabel('Decomposed','FontWeight','bold', 'Fontsize', 11)


% Multitpaper PSD on filter
s1=subplot(2,3,2); hold on
rectangle('Position', [phasefreq(2) ylim1 ampfreq(1)-phasefreq(2)  ylim2-ylim1], 'FaceColor', [0.7 0.7 0.7], 'lineStyle','none')
rectangle('Position', [ampfreq(2) ylim1 xlim2-ampfreq(2)  ylim2-ylim1], 'FaceColor', [0.7 0.7 0.7], 'lineStyle','none')
plot(fyt , 10*log10(SyMt), 'color', 'k', 'linewidth', 1.3)
%plot(f2, 10*log10(Sso_mt_filt+SAl_mt_filt), 'color', 'g', 'linewidth', 1.3)

plot(f2, 10*log10(Sso_mt_filt),'-.', 'color', 'b', 'linewidth', 1.8)
plot(f2, 10*log10(SAl_mt_filt),'-.', 'color', 'r', 'linewidth', 1.8)

ylabel('dB')
set(gca,'xtick',[ 0 ,10, 20],'ytick',[ -40 ,-20, 0])
ylim([ylim1 ylim2])
xlim([xlim1 xlim2])
box on;
title('Multitaper PSD')

% Multitpaper PSD on decomposed
s2=subplot(2,3,5); hold on
plot(fyt , 10*log10(SyMt), 'color', 'k', 'linewidth', 1.3)
set(gca,'xtick',[ 0 ,10, 20],'ytick',[ -40 ,-20, 0])
%plot(f1, 10*log10(SAl_mt_osc+Sso_mt_osc), 'color', 'g', 'linewidth', 1.3)

plot(f1, 10*log10(Sso_mt_osc),'-.', 'color', 'b', 'linewidth', 1.8)
plot(f1, 10*log10(SAl_mt_osc),'-.', 'color', 'r', 'linewidth', 1.8)
linkaxes([s1,s2])
ylim([ylim1 ylim2])
xlim([xlim1 xlim2])
xlabel('(Hz)')
ylabel('dB')
box on;

phi_0=pi/3;
pac_true = get_pac_from_kphi(K_mod_th,phi_0,Npb);

% PAC standard
pam1=subplot(2,3,3); hold on

plot(pbins, pac_true, 'color', 'k', 'linewidth',  2)
plot(pbins2, pa_usual, 'color', [0 0.5 0], 'linewidth',  2)
set(gca,'xtick',[-pi,0, pi],'XTickLabel',{'-\pi','0','\pi'})
axis tight
title('PAM')
axis tight
legend('True','Estimate')
set(gca,'xtick',[-pi,0, pi],'XTickLabel',{'-\pi','0','\pi'},'ytick',[0.5, 1, 1.5],'YTickLabel',{'0.5','1','1.5'})
box on;

% Parametric PAC
pam2=subplot(2,3,6); hold on
ciplot(pac_param_inf,pac_param_sup,pbins,[0.75 0.75 0.75],1)


pp1=plot(pbins, pac_true, 'color', 'k', 'linewidth', 2);
pp2=plot(pbins,pac_param_inf,'color', [0.75 0.75 0.75], 'linewidth', 2);
pp3=plot(pbins,pac_param_mean, 'g', 'linewidth', 2);
plot(pbins,pac_true, 'k', 'linewidth', 2)


xlabel('(rad)')
legend([pp1,pp2,pp3],{'True','Resampled','Estimate'})
set(gca,'xtick',[-pi,0, pi],'XTickLabel',{'-\pi','0','\pi'},'ytick',[0.5, 1, 1.5],'YTickLabel',{'0.5','1','1.5'})
box on;

linkaxes([pam1 pam2])
xlim([-pi pi])
ylim([0 0.3])

h=gcf;
set(h,'Position',[50 50 3*250 2*250]);
set(h,'PaperOrientation','landscape');
print(gcf, '-dpdf', 'figure5_tmp.pdf')


%% Plots Figure 6

figure(2)

f61=subplot(2,1,1);cla; hold on
plot(f_y, 10*log10(H_i_osc(1,:)), 'color', [0 0 1], 'linewidth', 1.3)
plot(f_y, 10*log10(H_i_osc(2,:)), 'color', [1 0 0], 'linewidth', 1.3)
xlabel('(Hz)')
ylabel('dB')
title('Parametric PSD')

f62=subplot(2,1,2);cla; hold on
plot(fyt , 10*log10(SyMt), 'color', 'k', 'linewidth', 1.3)
plot(f2, 10*log10(Sso_mt_osc), 'color', [0 0 1], 'linewidth', 1.3)
plot(f2, 10*log10(SAl_mt_osc), 'color', [1 0 0], 'linewidth', 1.3)
xlabel('(Hz)')
ylabel('dB')
legend('raw','filtered')
title('Multitaper PSD')
linkaxes([f61 f62])
xlim([0 35])
ylim([-50 20])

h=gcf;
set(h,'Position',[50 50 800 550]);
set(h,'PaperOrientation','landscape');
print(gcf, '-dpdf', 'figure6_tmp.pdf')
