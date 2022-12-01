addpath(genpath('../'))


Fs=250;


dpath ='/media/hdhs/data/eeg_propofol_data/dss_lin_pac_3D/';
bhvr_path = '/media/hdhs/data/eeg_propofol_data/behavioral_data/';

%subjectId={'02','04','13', '03'};
subjectId={'04'};

xlim1 = 30;
xlim2 = 140;

cramge = 0.25/pi;
cmin = 1/(2*pi)-cramge;
cmax = 1/(2*pi)+cramge;

mimin = 0;
mimax=0.05;

for subjectCur=1:length(subjectId)
    %% Load and Grasp behavioral data
    name_bhvr = [bhvr_path, 'eeganes_bhvr_summary_sub',subjectId{subjectCur}];
    load([name_bhvr,'.mat'])
    
    tprop_min = bhvr_summary.prop_time/(bhvr_summary.Fs_behav*60);
    tburs_min = bhvr_summary.t_burst/(bhvr_summary.Fs_behav*60);
    tverb_min = bhvr_summary.t_verb/(bhvr_summary.Fs_behav*60);
    
    %% Load and Grasp dsslin processed data
    name_scur = [dpath, 'dss_lin_pac_new_ss_filtering_sub',subjectId{subjectCur} ];
    load([name_scur,'.mat'])
    
    % 'Distance
    MI_type = 'KL';
    
    % Processing Parameters
    tVect = dsslin_pac.tVect;
    pbins_dss = dsslin_pac.pbins;
    Kappa1 = dsslin_pac.Kappa1;
    Kappa2 = dsslin_pac.Kappa2;
    alpha_ci = dsslin_pac.alpha_ci;
    
    % Goodness of fit
    R2_mis=dsslin_pac.R2_mis;
    
    % Raw estimates
    K_mis = dsslin_pac.K_mis;150
    Phi_mis = dsslin_pac.Phi_mis;
    PAC_raw = dsslin_pac.PAC_raw;
    MI_KL_mis = dsslin_pac.MI_KL_mis;
    MI_L1_mis = dsslin_pac.MI_L1_mis;
    
    % Filtered estimates
    K_mis_filt = dsslin_pac.K_mis_filt;
    Phi_mis_filt = dsslin_pac.Phi_mis_filt;
    PAC_filt = dsslin_pac.PAC_filt;
    MI_KL_mis_filt = dsslin_pac.MI_KL_mis_filt;
    MI_L1_mis_filt = dsslin_pac.MI_L1_mis_filt;
    
    % R2
    [R2_mean, R2_inf, R2_sup] = get_mis(R2_mis);
    xlim([tVect(1) tVect(end)])
    % Ks
    [K_mean, K_inf, K_sup] = get_mis(K_mis);
    [K_mean_filt, K_inf_filt, K_sup_filt] = get_mis(K_mis_filt);
    
    %Phis
    [Phi_mean, Phi_inf, Phi_sup] = get_mis(Phi_mis);
    [Phi_mean_filt,Phi_inf_filt, Phi_sup_filt] = get_mis(Phi_mis_filt);
    
    % MIs
    [MI_KL_mean, MI_KL_inf, MI_KL_sup] = get_mis(MI_KL_mis);
    [MI_L1_mean, MI_L1_inf, MI_L1_sup] = get_mis(MI_L1_mis);
    
    [MI_KL_mean_filt, MI_KL_inf_filt, MI_KL_sup_filt] = get_mis(MI_KL_mis_filt);
    [MI_L1_mean_filt, MI_L1_inf_filt, MI_L1_sup_filt] = get_mis(MI_L1_mis_filt);
    
    %% Grasp standard processed data
    
    load([dpath,'standard_pac_sub', subjectId{subjectCur}])
    Fs = standard_pac.pac{1,1}.Fs;
    
    plottedLength = [6,120];
    Nplots = length(plottedLength);
    corresLengthId = zeros(1,Nplots);
    
    for twsec = 1:Nplots
        corresLengthId(1,twsec) = find(standard_pac.length_wind_sec_tot == plottedLength(1,twsec) );
    end
    
    Npb2 = standard_pac.Npb;
    pbins2 = linspace(-pi,pi,Npb2-1);
    
    
    %% Plot params
   
    
    
    
    if strcmp(MI_type, 'KL')
        MI_plot_inf  = MI_KL_inf;
        MI_plot_sup  = MI_KL_sup;
        MI_plot_mean = MI_KL_mean;
        MI_plot_inf_filt  = MI_KL_inf_filt;
        MI_plot_sup_filt  = MI_KL_sup_filt;
        MI_plot_mean_filt = MI_KL_mean_filt;
    elseif strcmp(MI_type, 'L1')
        MI_plot_inf  = MI_L1_inf;
        MI_plot_sup  = MI_L1_sup;
        MI_plot_mean = MI_L1_mean;
        MI_plot_inf_filt  = MI_L1_inf_filt;
        MI_plot_sup_filt  = MI_L1_sup_filt;
        MI_plot_mean_filt = MI_L1_mean_filt;
    else
        error('Invalid MI_type')
    end
    
    
    
    %% Plots
    
    figmain=figure;
    
    % Propofol Concentration
    subplot(Nplots+3, 2, 2);hold on
    plot(tprop_min,bhvr_summary.prop_conc,'Color', 'k')
    ylabel('ug/ mL')
    box on 
    xlim([xlim1 xlim2])
    
    % Probability of response
    subplot(Nplots+3, 2, 1);hold on
    %ciplot(bhvr_summary.p05_burst,bhvr_summary.p95_burst,tburs_min,[0 0 1], 0.5); hold on
    %p1=plot(tburs_min,bhvr_summary.p_burst, 'linewidth',2, 'color' , 'b' );
    ciplot(bhvr_summary.p05_verb,bhvr_summary.p95_verb,tverb_min,[0.5 0.5 0.5],1); hold on
    p2=plot(tverb_min,bhvr_summary.p_verb, 'linewidth',2, 'color' , 'k' );
    ylabel('P(response)')
    %legend([p1,p2], {'Burst', 'Verbal'})
    box on ;xlim([xlim1 xlim2])
    
    % Plot Standard Processing
    for twsec_tmp = 1:Nplots
        
        twsec = corresLengthId(1,twsec_tmp);
        curpac = standard_pac.pac{1,twsec};
        twind_sec_cur = standard_pac.length_wind_sec_tot(1,twsec);
        tvecmin = 0.5*(curpac.startPointID+curpac.startPointID)/(60*Fs);
        pac_cur = curpac.pac;
        
        if strcmp(MI_type, 'KL')
            MI_cur = curpac.MI_KL;
            MI_perm_cur =  curpac.MI_KL_perm;
        elseif strcmp(MI_type, 'L1')
            MI_cur = curpac.MI_L1;
            MI_perm_cur =  curpac.MI_L1_perm; 
        else
            error('Invalid MI_type')
        end
        
        thr_mi = quantile(MI_perm_cur, 1-alpha_ci, 2);
        sig_mi = MI_cur>= thr_mi;
        
        % Standard modulogram
        subplot(Nplots+3, 2,  2*twsec_tmp+1)
        imagesc(tvecmin ,pbins2,pac_cur');
        col = redbluemap; colormap(col);
        colorbar('location', 'west')
        set(gca, 'clim', [cmin cmax]);
        axis tight
        xlim([xlim1 xlim2])
        ylabel([' Wind :',num2str(twind_sec_cur), 's'], 'FontWeight', 'bold');
        box on ;xlim([tVect(1) tVect(end)])
        %title(['Sub ', subjectId{subjectCur}, ' Wind ',num2str(twind_sec_cur), 's'])
        
         % Standard MI with significance
        subplot(Nplots+3, 2, 2*twsec_tmp+2); hold on
        p1 = plot(tvecmin,MI_cur, 'Color', 'k'); hold on
        p2 = scatter(tvecmin(find(sig_mi)) ,quantile(MI_cur,0.95)* ones(1,sum(sig_mi)), 10, [0.5 0.5 0.5],'*');
        ylabel(['MI ' , MI_type])
        box on ;xlim([tVect(1) tVect(end)])
        if twsec_tmp== 1
            legend([p1,p2], {'MI', ['p<', num2str(alpha_ci), '. rp= ' , num2str(size(MI_perm_cur,2))]});
        end
        xlim([xlim1 xlim2])
    end
    
    % PAM colorbar
    Xnum = Nplots+3;
    Ynum = 2;
    pos1=get(subplot(Xnum,Ynum,1+Ynum), 'Position');
    pos2=get(subplot(Xnum,Ynum,(Xnum-2)*Ynum), 'Position');
    widtBar=pos2(1)/50;
    middleFig=(pos1(1)+pos2(1)+pos2(3))/2;%-widtBar/1000;
    highBar= abs(pos1(2)+pos1(4)-(pos2(2)));
    cax = colorbar('Position', [middleFig pos2(2)  widtBar highBar]); %caxis([0.5 1.5]);
    set(gca, 'clim', [cmin, cmax]); %0.2-1.8
    set(cax  , 'Ticks', [0 0.1 0.15 0.2 0.25]);
    ylabel(cax,'Alpha power Distribution');
    uu=colormap(cax,col);
    
    % Plot DSS Raw
    % DSS Raw Modulogram
    subplot(Nplots+3, 2, 2*Nplots+3);hold on
    imagesc(tVect ,pbins_dss,PAC_raw);
    col = redbluemap; colormap(col);
    colorbar('location', 'west')
    set(gca, 'clim', [cmin cmax]);
    axis tight
    ylabel('DSS Lin Raw', 'FontWeight', 'bold');
    box on ;xlim([xlim1 xlim2])
    
    % DSS Raw MI
    subplot(Nplots+3, 2, 2*Nplots+4);hold on
    ciplot(MI_plot_inf,MI_plot_sup,tVect,[0.5 0.5 0.5])
    m1 = plot(tVect,MI_plot_sup, '-', 'Color', [0.5 0.5 0.5]);
    m2 = plot(tVect,MI_plot_mean, '-', 'Color', 'k');
    ylabel(['MI ' , MI_type])
    legend([m1,m2], {[num2str(floor(100*(1-alpha_ci))), '% CI, rs=', num2str(Kappa1),'\times', num2str(Kappa2)],'Mean'})
    box on ;xlim([xlim1 xlim2])
    ylim([mimin mimax])
    
    
    % Plot DSS Filtered
    % DSS Filtered Modulogram
    subplot(Nplots+3, 2, 2*Nplots+5);hold on
    imagesc(tVect ,pbins_dss,PAC_filt);
    col = redbluemap; colormap(col);
    colorbar('location', 'west')
    set(gca, 'clim', [cmin cmax]);
    axis tight
    ylabel('DSS Lin Filt', 'FontWeight', 'bold');
    box on ;xlim([xlim1 xlim2])
    xlabel('[min]')
    
    % DSS Filtered MI
    subplot(Nplots+3, 2, 2*Nplots+6);hold on
    ciplot(MI_plot_inf_filt,MI_plot_sup_filt,tVect,[0.5 0.5 0.5])
    plot(tVect,MI_plot_mean_filt, '-', 'Color', 'k')
    xlim([tVect(1) tVect(end)])
    ylabel(['MI ' , MI_type])
    box on ;xlim([xlim1 xlim2])
    ylim([mimin mimax])
    xlabel('[min]')
    
    
    % Save Figure
    savePic=1;
    figmain.Renderer='Painters';
    set(figmain, 'position', [274    66   902   984])
    %set(figmain,'PaperOrientation','landscape');
    if savePic
        print(figmain, '-dpdf', ['pac_comparison_sub',subjectId{subjectCur} ,'pdf'])
    end
end







