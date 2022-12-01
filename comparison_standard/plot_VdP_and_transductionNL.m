%transduction_nonlinearities.m
%VanDerPole.m

figure;

for include_fast_osc = [0,1]
    
    if include_fast_osc
        y_vdp = y_t_mod_vdp;
    else
        y_vdp =y_t_no_mod_vdp;
    end
    
    modelOsc_osc_best = best_models_vdp{1,include_fast_osc+1}.modelOsc_osc_best;
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
    pac_std_VdP=standard_processing_vdp{1,include_fast_osc+1}.UsualPAM;
    
    % VdP Simulated Signal
    subplot(3,6,6*include_fast_osc+1) ;hold on
    plot(ta_sec,y_vdp, 'k')
    plot(ta_sec,x_slow_vdp, 'b')
    if include_fast_osc ==1
        plot(ta_sec,x_fast_vdp , 'r')
    end
    axis tight
    ylim0=ylim;
    ylabel('[a.u]')
    legend('Signal Tot', 'Slow','Fast')
    box on
    
    % VdP information criterion
    subplot(3,6,6*include_fast_osc+4) ;hold on
    plot(1:(harmTest) ,dAIC(1,:), 'color', [1 0.6 0], 'linewidth', 2)
    plot(1:(harmTest) ,dBIC(1,:), '-.', 'color', [1 0.6 0], 'linewidth', 2)
    plot(1:(harmTest) ,dAIC(2,:), 'color', [0.5 0.5 0.5], 'linewidth', 2)
    plot(1:(harmTest) ,dBIC(2,:), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
    legend('AIC : no fast', 'BIC : no fast', 'AIC : with fast', 'BIC : with fast')
    axis tight
    box on

    % VdP Filtered Signal
    subplot(3,6,6*include_fast_osc+2) ;hold on
    plot(ta_sec,x_slow_filt , 'b')
    plot(ta_sec,x_fast_filt , 'r')
    box on
    ylim(ylim0)
    legend(['Slow:', num2str(slow_freq_VdP(1)),'-',num2str(slow_freq_VdP(2)),'Hz'], ['Fast:', num2str(fast_freq_VdP(1)),'-',num2str(fast_freq_VdP(2)),'Hz'])
    
    % VdP state space decomposed
    subplot(3,6,6*include_fast_osc+5) ;hold on
    plot(ta_sec,x_t_n_slow , 'b')
    if length(NNharm) >1
        plot(ta_sec,x_t_n_fast , 'r')
    end
    box on
    ylim(ylim0)
    
    % Standard PAC
    subplot(3,6,6*include_fast_osc+3); hold on
    l31=plot(pbins2, pac_std_VdP, 'g', 'linewidth', 2);
    if include_fast_osc
        l32=plot(pbins, pac_true_vdp, 'Color', 'k', 'linewidth', 2);
    else
        l32=plot(pbins, (1/(2*pi)) * ones(1, length(pbins)), 'Color', 'k', 'linewidth', 2);
    end
    box on
    axis tight
    ylim1 = ylim;
    xlim([-pi pi])
    legend([l31, l32],{'True','Estimate'})
    
    % SSP PAC
    subplot(3,6,6*include_fast_osc+6); hold on
    if length(NNharm)>1
        ciplot(pac_param_inf,pac_param_sup,pbins,[0.75 0.75 0.75],1)
        l43 =plot(pbins,pac_param_inf,'color', [0.75 0.75 0.75], 'linewidth', 2);
        l42 =plot(pbins,pac_param_mean, 'g', 'linewidth', 2);
        plot(pbins,pac_true_vdp, 'k', 'linewidth', 2)
        
    end
    ylim(ylim1)
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
    ylim(ylim1)
    xlim([-pi pi])
    
end

% Signal with non linear trasduction
subplot(3,6,13); hold on
plot(ta_sec2, y_t_transdn-mean(y_t_transdn),'k')
plot(ta_sec2, x_t_transdn-mean(x_t_transdn),'b')
box on
axis tight
ylim0 = ylim;
xlabel('[sec]')
ylabel('[a.u]')
legend('Signal Tot', 'Transducted')

% Filtered Signal with non linear trasduction
subplot(3,6,14); hold on
plot(ta_sec2, x_filt_slow_transdn,'b')
plot(ta_sec2, x_filt_alpha_transdn,'r')
legend(['Slow:', num2str(slow_freq_transdn(1)),'-',num2str(slow_freq_transdn(2)),'Hz'], ['Fast:', num2str(fast_freq_transdn(1)),'-',num2str(fast_freq_transdn(2)),'Hz'])

ylim(ylim0)
box on
xlabel('[sec]')

% Standard PAC
subplot(3,6,15); hold on
plot(pbins, pac_true_transdn,'k','linewidth',2)
plot(pbins2, pac_standard_transdn,'g','linewidth',2)
axis tight
ylim1 = ylim;
box on
xlabel('[rad]')
legend({'True','Estimate'})

% Information criteria
subplot(3,6,16); hold on
plot(1:harmTest, dAIC_transdn(1,:),'color',[1 0.6 0], 'linewidth', 2)
plot(1:harmTest, dBIC_transdn(1,:),'-.','color',[1 0.6 0], 'linewidth', 2)

plot(1:harmTest, dAIC_transdn(2,:),'color',[0.5 0.5 0.5], 'linewidth', 2)
plot(1:harmTest, dBIC_transdn(2,:),'-.','color',[0.5 0.5 0.5], 'linewidth', 2)
axis tight;box on
legend('AIC : no fast', 'BIC : no fast', 'AIC : with fast', 'BIC : with fast')
xlabel('#Slow Harmonics')

% State space decomposed signal with non linear trasduction
subplot(3,6,17); hold on
plot(ta_sec2,x_t_n_slow_transdn,'b')
if length(NNharm_transdn)>1
   plot(ta_sec2,x_t_n_fast_transdn,'r') 
end
ylim(ylim0);  box on
xlabel('[sec]')


% SSP pac
subplot(3,6,18); hold on
plot(pbins, pac_true_transdn,'k','linewidth',2) 
if length(NNharm_transdn)>1

plot(pbins, pac_ssp_transdn,'g','linewidth',2)
else
   
   plot(pbins, pac_ssp_transdn,':g','linewidth',2)
end
xlim([-pi pi])
ylim(ylim1);  box on
xlabel('[rad]')
legend({'True','Estimate'})

h=gcf; savePic =1;
h.Renderer='Painters';
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position', [68         210        1107         524]);
if savePic
    print(h, '-dpdf', ['VdP_and_transduction_nonlinearities' ,'pdf'])
end