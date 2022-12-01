addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../ssp_pac'))
addpath(genpath('../usual_pam'))

% Initializarion
Fs = 250;
ta_sec = 0:1/Fs:6;
doPlot = 0;
Npb=18;


%% Generate Van der Pol oscillator (and modulated alpha)

% VdP parameters
w0 =5; Mu = 5;
odeVDP = @(t,y) vanderpoldemo2(t,y,Mu,w0);

% Grasp Vdp Oscilator
y_tmp = ode2(odeVDP,ta_sec,[2 0]);
x_slow_vdp = y_tmp(:,1)' -mean(y_tmp(:,1))';

% Modulation
K_true_vdp = 0.9;
f_fast_vdp = 10; ampfast =0.3; amp_noise = 0.15;
x_fast_vdp  = ampfast*cos(2*pi*ta_sec*f_fast_vdp) .*  (1+ K_true_vdp * x_slow_vdp/max(x_slow_vdp) );
pac_true_vdp = get_pac_from_kphi(K_true_vdp,0,18);

% Signal Tot
y_t_mod_vdp    = x_fast_vdp' + x_slow_vdp' + normrnd(0,amp_noise, [length(ta_sec),1]);
y_t_no_mod_vdp = x_slow_vdp' + normrnd(0,amp_noise, [length(ta_sec),1]);

figure
subplot(2,1,1); hold on
plot(ta_sec,y_t_mod_vdp,'k', 'color', [0.5 0.5 0.5]);
plot(ta_sec,y_t_no_mod_vdp, 'color', [0.85 0.6 0]);
axis tight; box on
ylim0 = ylim;

subplot(2,1,2); hold on
plot(ta_sec,x_slow_vdp,'b');
plot(ta_sec,x_fast_vdp,'r');
ylim(ylim0); box on
%% State-Space Oscillator Decomposition

em_its=300;
harmTest = 8;

f_tot_vdp = [2 13];
a_tot_vdp = [0.98 0.98];
sigma2_tot_vdp = [0.01 0.01];
best_models_vdp =cell(1,2);

% Include or not the alpha componemt in the state space oscillator model
for include_fast_osc = [0,1]
    % Launch a parfor loop on the number of harmonic
    
    if include_fast_osc
        y = y_t_mod_vdp;
    else
        y =y_t_no_mod_vdp;
    end
    
    modelOsc_harm = cell(1,harmTest);
    
    parfor harmtest_cur = 1:(harmTest)
        disp(['Vdp : ', 'With Alpha=', num2str(include_fast_osc), ' Harm#', num2str(harmtest_cur), '/' , num2str(harmTest)])
        
        NNharm = [harmtest_cur,1];
        
        % 1 osc + harm
        init_params=struct();
        init_params.f_init     =  f_tot_vdp(1,1:1)
        init_params.a_init     =  a_tot_vdp(1,1:1)
        init_params.sigma2_init=  sigma2_tot_vdp(1,1:1)
        init_params.R        =  5;
        init_params.NNharm      = NNharm(1,1:1);
        modelOsc_1osc=ssp_decomp(y,Fs,1,length(y),em_its,eps, init_params,'',doPlot);
        
        % 2 osc + harm
        init_params=struct();
        init_params.f_init     =  f_tot_vdp(1,1:2)
        init_params.a_init     =  a_tot_vdp(1,1:2)
        init_params.sigma2_init=  sigma2_tot_vdp(1,1:2)
        init_params.R        =  5;
        init_params.NNharm      =  NNharm(1,1:2);
        modelOsc_2osc=ssp_decomp(y,Fs,1,length(y),em_its,eps, init_params,'',doPlot);
        
        modelOsc_harm{:,harmtest_cur}={modelOsc_1osc;modelOsc_2osc};
        
    end
    
    
    
    
    % Compute Information criterion to Select the best model
    lltot = zeros(size(modelOsc_harm{1,1},1),harmTest);
    Kparam = zeros(size(modelOsc_harm{1,1},1),harmTest);
    
    for harmtest_cur = 1:(harmTest)
        for j = 1:size(modelOsc_harm{1,1},1)
            lltot(j,harmtest_cur) = modelOsc_harm{1,harmtest_cur}{j}.res{1,1}.ll - ( 0.5*length(y) *log(2*pi));
            
            Kparam_tmp = 3*  sum(modelOsc_harm{1,harmtest_cur}{j}.init_params.NNharm) +1;
            Kparam(j,harmtest_cur) = Kparam_tmp;
        end
    end
    
    AIC = -2*lltot + 2 * Kparam;
    BIC = -2*lltot + 2 * Kparam * log(length(y));
    
    dAIC = AIC- min(min(AIC));
    dBIC = BIC- min(min(BIC));
    
    pdAIC = exp(-0.5*dAIC)./sum(sum(exp(-0.5*dAIC)));
    pdBIC = exp(-0.5*dBIC)./sum(sum(exp(-0.5*dBIC)));
    
    IC = dAIC;
    [mIC, locMax_1_tmp] = min(dAIC);
    [~, locMax_2]   = min(mIC);
    
    
    locX = locMax_1_tmp(locMax_2);
    locY = locMax_2;
    
    modelOsc_osc_best = modelOsc_harm{1,locY}{locX};
    NNharm = modelOsc_osc_best.init_params.NNharm;
    
    best_models_vdp{1,include_fast_osc+1}.modelOsc_osc_best = modelOsc_osc_best;
    best_models_vdp{1,include_fast_osc+1}.NNharm=NNharm;
    best_models_vdp{1,include_fast_osc+1}.dAIC=dAIC;
    best_models_vdp{1,include_fast_osc+1}.dBIC=dBIC;
      
end


%% Use a linear regression to estimate Coupling and resample data.
for include_fast_osc = [0,1]
    disp(['Vdp : ', 'With Alpha=', num2str(include_fast_osc)])
    Kappa_tot = [200;200];
    modelOsc_osc_best = best_models_vdp{1,include_fast_osc+1}.modelOsc_osc_best;
    slowID= 1;
    fastID=modelOsc_osc_best.init_params.NNharm(1,1)+1;
    
    Nlag = 0; NlinConst = 100;
    NNharm = modelOsc_osc_best.init_params.NNharm;
    
    if length(NNharm) >1
        [regressed_pac, SS_tot, X_t_tot] = ssp_pac (modelOsc_osc_best, Kappa_tot, slowID, fastID, Nlag,NlinConst);
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
        
        best_models_vdp{1,include_fast_osc+1}.pac_param_inf=pac_param_inf;
        best_models_vdp{1,include_fast_osc+1}.pac_param_sup=pac_param_sup;
        best_models_vdp{1,include_fast_osc+1}.pac_param_mean=pac_param_mean;
        
    else
        
        best_models_vdp{1,include_fast_osc+1}.pac_param_inf=get_pac_from_kphi(0,0,Npb);
        best_models_vdp{1,include_fast_osc+1}.pac_param_sup=get_pac_from_kphi(0,0,Npb);
        best_models_vdp{1,include_fast_osc+1}.pac_param_mean=get_pac_from_kphi(0,0,Npb);
        
    end
    
    x_t_n_tot  = modelOsc_osc_best.res{1,1}.x_t_n(1:2:end, 2:end);
    x_t_n_slow = sum(x_t_n_tot(1:NNharm(1,1),:),1);
    if length(NNharm) >1
        x_t_n_fast = x_t_n_tot(NNharm(1,1)+1,:);
    else
        x_t_n_fast = zeros(size(x_t_n_slow));
        
        
    end
    ta_sec = modelOsc_osc_best.res{1,1}.ta;
    
    best_models_vdp{1,include_fast_osc+1}.x_t_n_slow=x_t_n_slow;
    best_models_vdp{1,include_fast_osc+1}.x_t_n_fast=x_t_n_fast;
    
end

%% Standard Processing
standard_processing_vdp = cell(1,2);
for include_fast_osc = [0,1]
    
    if include_fast_osc
        y = y_t_mod_vdp;
    else
        y =y_t_no_mod_vdp;
    end
    
    % Bandpass of the solution and Usual PAC
    slow_freq_VdP = [0.1 1.5]; fast_freq_VdP = [6 14]; order=100;
    [x_slow_filt, tail1] =quickbandpass(y',Fs,slow_freq_VdP, order);   %Bandpass at the phase frequency
    [x_fast_filt, tail2]=quickbandpass  (y',Fs,fast_freq_VdP, order);
    
    A_fast   = abs(hilbert(x_fast_filt));
    Phi_slow = angle(hilbert(x_slow_filt));
    
    pbins2 = linspace(-pi, pi, Npb-1);
    UsualPAM = phaseamp(A_fast,Phi_slow,pbins)  ;
    UsualPAM = UsualPAM ./ sum(UsualPAM) * (Npb/(2*pi));
    
    standard_processing_vdp{1,include_fast_osc+1}.x_slow_filt=x_slow_filt;
    standard_processing_vdp{1,include_fast_osc+1}.x_fast_filt=x_fast_filt;
    standard_processing_vdp{1,include_fast_osc+1}.UsualPAM=UsualPAM;
    
end



%%



%%


%% Plot PAC

figure;
withOrWithout = {'no','with'};

for include_fast_osc = [0,1]
    
    if include_fast_osc
        y = y_t_mod_vdp;
    else
        y =y_t_no_mod_vdp;
    end
    
    best_models_vdp{1,include_fast_osc+1}.modelOsc_osc_best = modelOsc_osc_best;
    NNharm = best_models_vdp{1,include_fast_osc+1}.NNharm;
    dAIC=best_models_vdp{1,include_fast_osc+1}.dAIC;
    dBIC=best_models_vdp{1,include_fast_osc+1}.dBIC;
    pac_param_inf=best_models_vdp{1,include_fast_osc+1}.pac_param_inf;
    pac_param_sup=best_models_vdp{1,include_fast_osc+1}.pac_param_sup;
    pac_param_mean=best_models_vdp{1,include_fast_osc+1}.pac_param_mean;
    x_t_n_slow=best_models_vdp{1,include_fast_osc+1}.x_t_n_slow;
    x_t_n_fast=best_models_vdp{1,include_fast_osc+1}.x_t_n_fast;
    
    
    x_slow_filt=standard_processing_vdp{1,include_fast_osc+1}.x_slow_filt;
    x_fast_filt=standard_processing_vdp{1,include_fast_osc+1}.x_fast_filt;
    UsualPAM=standard_processing_vdp{1,include_fast_osc+1}.UsualPAM;
    
    
    
    p0=subplot(2,6,6*include_fast_osc+1) ;hold on
    plot(ta_sec,y, 'k')
    plot(ta_sec,x_slow_vdp, 'b')
    if include_fast_osc ==1
        plot(ta_sec,x_fast_vdp , 'r')
    end
    ylabel('a.u')
    xlabel('sec')
    title(['Simulated Signal'])
    %: ',withornot{1,putAlpha+1}, ' fast component'
    legend('Signal Tot', 'Slow','Fast')
    box on
    
    subplot(2,6,6*include_fast_osc+4) ;hold on
    plot(1:(harmTest) ,dAIC(1,:), 'color', [1 0.6 0], 'linewidth', 2)
    plot(1:(harmTest) ,dBIC(1,:), '-.', 'color', [1 0.6 0], 'linewidth', 2)
    plot(1:(harmTest) ,dAIC(2,:), 'color', [0.5 0.5 0.5], 'linewidth', 2)
    plot(1:(harmTest) ,dBIC(2,:), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
    xlabel('Number of Harmonics')
    ylabel('IC')
    legend('AIC : no fast', 'BIC : no fast', 'AIC : with fast', 'BIC : with fast')
    withornot={'without', 'with'};
    axis tight
    box on
    title(['IC'])
    
    p1=subplot(2,6,6*include_fast_osc+2) ;hold on
    plot(ta_sec,x_slow_filt , 'b')
    plot(ta_sec,x_fast_filt , 'r')
    ylabel('a.u')
    xlabel('sec')
    title(['Filtered'])
    box on
    
    p2=subplot(2,6,6*include_fast_osc+5) ;hold on
    plot(ta_sec,x_t_n_slow , 'b')
    if length(NNharm) >1
        plot(ta_sec,x_t_n_fast , 'r')
    end
    ylabel('a.u')
    xlabel('sec')
    title(['Best Osc/Harm'])
    box on
    linkaxes([p0,p1,p2])
    ylim([-2.1 2.1])
    axis tight
    
    p3=subplot(2,6,6*include_fast_osc+3); hold on
    hold on
    l31=plot(pbins, UsualPAM, 'g', 'linewidth', 2);
    title(['Standard PAC'])
    if include_fast_osc
        l32=plot(pbins, pac_true_vdp, 'Color', 'k', 'linewidth', 2);
    else
        l32=plot(pbins, (1/(2*pi)) * ones(1, length(pbins)), 'Color', 'k', 'linewidth', 2);
    end
    box on
    xlabel('rad')
    ylim([0 0.35])
    xlim([-pi pi])
    legend([l31, l32],{'True','Estimate'})
    
    p4=subplot(2,6,6*include_fast_osc+6); hold on
    %truePAC=(1/(2*pi))*(1+(sin(pi/Npb)/(pi/Npb)) * K_mod *cos(pbins+0));
    
    if length(NNharm)>1
        ciplot(pac_param_inf,pac_param_sup,pbins,[0.75 0.75 0.75],1)
        l43 =plot(pbins,pac_param_inf,'color', [0.75 0.75 0.75], 'linewidth', 2);
        l42 =plot(pbins,pac_param_mean, 'g', 'linewidth', 2);
        plot(pbins,pac_true_vdp, 'k', 'linewidth', 2)
        
    end
    ylim([0 0.35])
    xlim([-pi pi])
    
    
    
    if include_fast_osc
        l41=plot(pbins, pac_true_vdp, 'Color', 'k', 'linewidth', 2);
    else
        l41=plot(pbins, (1/(2*pi)) * ones(1, length(pbins)), 'Color', 'k', 'linewidth', 2);
    end
    
    if length(NNharm)>1
        l42=plot(pbins, pac_param_mean, 'g', 'linewidth', 2);
    else
        l42=plot(pbins, pac_param_mean,':', 'color','g', 'linewidth', 2);
    end
    
    if length(NNharm)>1
        legend([l41,l42,l43],{'True','Estimate', 'Resampled'})
    else
        legend([l41,l42],{'True','Estimate'})
    end
    
    box on
    linkaxes([p3,p4])
    ylim([0 0.5])
    xlim([-pi pi])
    %legend('True','Boostrap','Estimate')
    xlabel('rad')
    title(['Parametric  PAC'])
    
    
end



h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [-1.7 0 13.8 2*2.3]);
%print([spath, 'bestOscVsFilt_', withornot{1,putAlpha+1}, '_fast'], '-dpdf')


%% Side function
function dydt = vanderpoldemo2(t,y,epsilon,w_0)
dydt = [y(2); epsilon*w_0*(1-y(1)^2)*y(2)-w_0^2*y(1)];
end
