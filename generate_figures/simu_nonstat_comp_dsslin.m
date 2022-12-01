%% Load Data
addpath(genpath('./../'))
lpath = '/media/hdhs/data/non_stationnary_simulations/';
load([lpath, 'non_stationnary_sim_signal_tot_final_oscdecomp.mat'])
load([lpath, 'non_stationnary_sim_signal_tot_final_oscdecomp_dss_decomp.mat'])



%% Compare Processed
saveFig = 1;

cramge = 0.25/pi;
cmin = 1/(2*pi)-cramge;
cmax = 1/(2*pi)+cramge;
Fs = 250;
mimax = 0.1;
%for nwidth= 1 :length(delta_f_gen)
%    for cur_amp = 1:length(sigma_slow)

lengthSigToPlot = [1200,300,120];
lengthWinToPlot = [6,120];

for nwidth= 1
    for cur_amp = 1
        
        figmain=figure;
        for tsignal = 1:length(lengthSigToPlot)
            
            SigLenId = find(T_tot_sec_tot   ==lengthSigToPlot(tsignal));
            
            tt = signal_tot{1,SigLenId}.tt/(60*Fs);
            K_cur =signal_tot{1,SigLenId}.Kmod;
            K_cur_clean=signal_tot{1,SigLenId}.Kmod_clean;
            P_cur=signal_tot{1,SigLenId}.Phimod;
            P_cur_clean= signal_tot{1,SigLenId}.Phimod_clean;
            y_t_cur = squeeze(signal_tot{1,SigLenId}.y_tot(nwidth,cur_amp,:));
            
            % Plot sample signal
            if tsignal == 1
                tsamplePlot = (locKmaxId-1.5*Fs):(locKmaxId+1.5*Fs);
                tlimmax =  [tt(1), tt(end)];
                
                subplot(2+length(lengthWinToPlot),(length(lengthSigToPlot)+1)*2, 2); hold on
                plot(tt(tsamplePlot),y_t_cur(tsamplePlot),'k');
                axis tight
                box on
                
                scalemuVolt = 10;
                scaleTime = 1;
                tt1 = tt(tsamplePlot(1));
                tt2 = tt(tsamplePlot(scaleTime*Fs));
                
                ylimTmp=ylim;
                plot([tt1; tt1], [min(y_t_cur(tsamplePlot)); min(y_t_cur(tsamplePlot))+scalemuVolt], '-k',  [tt1; tt2], [ min(y_t_cur(tsamplePlot));  min(y_t_cur(tsamplePlot))], '-k', 'LineWidth', 2)
                text(tt1,min(y_t_cur(tsamplePlot))-(ylimTmp(2)-ylimTmp(1))/10 ,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
                set(gca, 'Visible', 'off')
                title('dsdssd')
                
                subplot(2+length(lengthWinToPlot),(length(lengthSigToPlot)+1)*2, 2+(length(lengthSigToPlot)+1)*2); hold on
                title(['\color{blue}\Delta f=', num2str(delta_f_gen(1,nwidth )), ' \sigma_{s}=' , num2str(sigma_slow(1,cur_amp)), ' \sigma_{f}=' , num2str(sigma_fast(1,cur_amp))])
            end
            
            
            % Plot modulation parameters
            subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1),2+2*tsignal); hold on
            plot(tt,K_cur_clean, 'k','linewidth',1.5)
            xlim(tlimmax)
            xlabel('[min]')
            box on
            
            subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1),1+2*tsignal); hold on
            plot(tt,P_cur_clean, 'k','linewidth',1.5)
            xlim(tlimmax)
            ylim([0 pi+0.1])
            box on
            xlabel('[min]')
            
            
            
            % Plot standard PAC
            for twindow = 1:length(lengthWinToPlot)
                
                WinLenId = find(pac_wind_sec_tot==lengthWinToPlot(twindow));
                
                % Grasp Parameters
                Nwind_cur = floor(T_tot_sec_tot(1,SigLenId)/pac_wind_sec_tot(1,WinLenId));
                tvecmin = (1:Nwind_cur)*pac_wind_sec_tot(1,WinLenId)/60;
                pac_cur   =  signal_tot{1,SigLenId}.pac_std{1,WinLenId}.pac_standard(:,:,nwidth,cur_amp);
                
                % Modulation Index And Permutted NI
                mi_kl_standard      = signal_tot{1,SigLenId}.pac_std{1,WinLenId}.mi_kl_standard(:,nwidth,cur_amp);
                mi_kl_standard_perm = signal_tot{1,SigLenId}.pac_std{1,WinLenId}.mi_kl_standard_perm(:,:,nwidth,cur_amp);
                
                % Point at which MI is considered significant
                alpha_ci = 0.05;
                th_kl = quantile(mi_kl_standard_perm, 1-alpha_ci, 2);
                sig_kl = mi_kl_standard> th_kl;
                
                % Plot Modulograms
                subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1) ,3 +(tsignal-1)*2 + (twindow-1)*2*(length(lengthSigToPlot)+1) + 2*(length(lengthSigToPlot)+1)*1); hold on
                imagesc(tvecmin ,pbins2,pac_cur');
                col = redbluemap; colormap(col);
                set(gca, 'clim', [cmin cmax]);
                axis tight
                
                
                if tsignal == 1
                    ylabel('[rad]')
                    colorbar('location', 'west')
                end
                
                % Plot MI,' || Signal Tot : ',num2str(T_tot_sec_tot(1,tsignal)), 's'
                subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1) ,4 +(tsignal-1)*2 + (twindow-1)*2*(length(lengthSigToPlot)+1) + 2*(length(lengthSigToPlot)+1)*1); hold on
                plot(tvecmin,mi_kl_standard, 'color', 'k')
                scatter(tvecmin(find(sig_kl)), 0.8*mimax* ones(1,length(find(sig_kl))), 40, 'm', '*')
                box on
                ylabel('MI')
                xlim([0 tvecmin(end)])
                ylim([0 mimax])
                
                if tsignal == 1 && twindow == 1
                    legend(['MI_{KL}. '],['p<', num2str(alpha_ci)]) %, 'rp=', num2str(size(mi_kl_standard_perm,2))
                end
                
            end
            
            
            % Plot Parametric DSS PAC
            dsslin_pac = dss_lin_pac_tot{SigLenId,nwidth,cur_amp} ;
            dsslin_pac_processed = dsslin_pac.dsslin_pac_processed;
            MI_KL_mean_raw  = dsslin_pac_processed.MI_KL_mis(1,:);
            MI_KL_inf  = dsslin_pac_processed.MI_KL_mis(2,:);
            MI_KL_sup  = dsslin_pac_processed.MI_KL_mis(3,:);
            
            MI_KL_mean_filt = dsslin_pac_processed.MI_KL_mis_filt(1,:);
            MI_KL_inf_filt  = dsslin_pac_processed.MI_KL_mis_filt(2,:);
            MI_KL_sup_filt  = dsslin_pac_processed.MI_KL_mis_filt(3,:);
            
            tVectmindss = dsslin_pac_processed.tVect/60;
            
            
            % Plot Filteres Modulogram
            subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1) ,3 +(tsignal-1)*2 + length(lengthWinToPlot)*2*(length(lengthSigToPlot)+1) + 2*(length(lengthSigToPlot)+1)*1); hold on
            imagesc(tVectmindss ,dsslin_pac_processed.pbins,dsslin_pac_processed.PAC_filt);
            col = redbluemap; colormap(col);
            if tsignal == 1
                ylabel('[rad]')
                colorbar('location', 'west')
            end
            set(gca, 'clim', [cmin cmax]);
            axis tight
            xlim([tVectmindss(1) tVectmindss(end)])
            ylabel('[rad]')
            xlabel('[min]')
            box on
            
            % Plot Filtered MI
            subplot(2+length(lengthWinToPlot),2*(length(lengthSigToPlot)+1) ,4+(tsignal-1)*2 + length(lengthWinToPlot)*2*(length(lengthSigToPlot)+1) + 2*(length(lengthSigToPlot)+1)*1); hold on
            ciplot(MI_KL_inf_filt,MI_KL_sup_filt,tVectmindss, [0.5 0.5 0.5]); hold on;
            p1 = plot(tVectmindss,MI_KL_mean_filt,'k');
            p2 = plot(tVectmindss,MI_KL_sup_filt,'color',[0.5 0.5 0.5]);
            box on
            xlim([tVectmindss(1) tVectmindss(end)])
            xlabel('[min]')
            ylim([0 mimax])
            ylabel('MI')
            
            
            
            if tsignal==1
                legend([p1,p2], {'MI_{KL}','CI'}) % ['rs=', num2str(dsslin_pac_processed.Kappa1),'x' , num2str(dsslin_pac_processed.Kappa2)]
            end
            
        end
        
        
        saveFig =1;
        set(gcf,'Position',[-9         356        1771         652])
        figmain.Renderer='Painters';
        orient(figmain,'landscape')
        pathDrop='';
        if saveFig
            namecur = ['non_stat_mod_comparison_dF',num2str(delta_f_gen(1,nwidth )),'Hz_sigId',num2str(cur_amp)];
            print(gcf, '-dpdf', [pathDrop,namecur,'fig1_aug'])
        end
    end
end




%% Plot signal Sample

figmain=figure;
for nwidth= 1:length(delta_f_gen)
    for cur_amp = 1:length(sigma_slow)
       
        tsignal = 2;
        
        
        SigLenId = find(T_tot_sec_tot   ==lengthSigToPlot(tsignal));
        
        tt = signal_tot{1,SigLenId}.tt/(60*Fs);
        K_cur =signal_tot{1,SigLenId}.Kmod;
        K_cur_clean=signal_tot{1,SigLenId}.Kmod_clean;
        P_cur=signal_tot{1,SigLenId}.Phimod;
        P_cur_clean= signal_tot{1,SigLenId}.Phimod_clean;
        y_t_cur = squeeze(signal_tot{1,SigLenId}.y_tot(nwidth,cur_amp,:));
        x_s_cur =  squeeze(signal_tot{1,SigLenId}.x_slow(nwidth,cur_amp,:));
        x_f_cur =  squeeze(signal_tot{1,SigLenId}.x_fast(nwidth,cur_amp,:));
        noi_cur = y_t_cur - x_s_cur - x_f_cur;
        
        
        tsamplePlot = (locKmaxId-1.5*Fs):(locKmaxId+1.5*Fs);
        
        scalemuVolt = 10;
        scaleTime = 1;
        tt1 = tt(tsamplePlot(1));
        tt2 = tt(tsamplePlot(scaleTime*Fs));
        
        subplot(3,length(sigma_slow)*length(delta_f_gen) , length(sigma_slow) * (nwidth-1) + cur_amp ); hold on
        plot(tt(tsamplePlot),y_t_cur(tsamplePlot),'k');
        
      
        
        axis tight
        box on
        
        if (nwidth == 1) && (cur_amp ==1)
        ylimTmp=ylim;
        plot([tt1; tt1], [min(y_t_cur(tsamplePlot)); min(y_t_cur(tsamplePlot))+scalemuVolt], '-k',  [tt1; tt2], [ min(y_t_cur(tsamplePlot));  min(y_t_cur(tsamplePlot))], '-k', 'LineWidth', 2)
        text(tt1,min(y_t_cur(tsamplePlot))-(ylimTmp(2)-ylimTmp(1))/10 ,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
        
        end
        ylim(ylimTmp)
        
        
        set(gca, 'Visible', 'off')
        
       
        
        subplot(3,length(sigma_slow)*length(delta_f_gen) , length(sigma_slow) * (nwidth-1) + cur_amp +length(sigma_slow)*length(delta_f_gen)); hold on
        plot(tt(tsamplePlot),noi_cur(tsamplePlot),'color', [0.5 0.5 0.5]);
        axis tight
        ylim(ylimTmp)
        set(gca, 'Visible', 'off')

        subplot(3,length(sigma_slow)*length(delta_f_gen) , length(sigma_slow) * (nwidth-1) + cur_amp +length(sigma_slow)*length(delta_f_gen)*2); hold on
        plot(tt(tsamplePlot),x_f_cur(tsamplePlot),'r');
        plot(tt(tsamplePlot),x_s_cur(tsamplePlot),'b');
        
        axis tight
        ylim(ylimTmp)
        poscur =get(gca, 'Position');
        set(gca, 'Visible', 'off')
         text(mean([tt1,tt1]), ylimTmp(2)    ,['\color{blue}\Delta f=', num2str(delta_f_gen(1,nwidth )), ' \sigma_{s}=' , num2str(sigma_slow(1,cur_amp)), ' \sigma_{f}=' , num2str(sigma_fast(1,cur_amp))],'HorizontalAlignment', 'left')
        
    end
end


saveFig =1;
set(gcf, 'position', [196         415        1120         488])

figmain.Renderer='Painters';
orient(figmain,'landscape')
pathDrop='';
if saveFig
    namecur = ['non_stat_mod_sample_signal'];
    print(gcf, '-dpdf', [pathDrop,namecur,'fig1_aug'])
end


