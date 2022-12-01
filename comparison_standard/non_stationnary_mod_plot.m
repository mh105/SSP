%% Load Data
addpath(genpath('./../'))
lpath = '/media/hdhs/data/non_stationnary_simulations/';
load([lpath, 'non_stationnary_sim_signal_tot_final_oscdecomp.mat'])
load([lpath, 'non_stationnary_sim_signal_tot_final_oscdecomp_dss_decomp.mat'])



%% Compare
saveFig = 1;

cramge = 0.25/pi;
cmin = 1/(2*pi)-cramge;
cmax = 1/(2*pi)+cramge;


for nwidth= 1 :length(delta_f_gen)
    for cur_amp = 1:length(sigma_slow)
        
        figure
        for tsignal = 1:length(T_tot_sec_tot)
            %% Plot Standard Processing
            for twindow = 1:length(pac_wind_sec_tot)
              
                % Grasp Parameters
                Nwind_cur = floor(T_tot_sec_tot(1,tsignal)/pac_wind_sec_tot(1,twindow));
                tvecsec = (1:Nwind_cur)*pac_wind_sec_tot(1,twindow);
                pac_cur   =  signal_tot{1,tsignal}.pac_std{1,twindow}.pac_standard(:,:,nwidth,cur_amp);
                
                % Modulation Index And Permutted NI
                mi_kl_standard      = signal_tot{1,tsignal}.pac_std{1,twindow}.mi_kl_standard(:,nwidth,cur_amp);
                mi_kl_standard_perm = signal_tot{1,tsignal}.pac_std{1,twindow}.mi_kl_standard_perm(:,:,nwidth,cur_amp);
                
                % Point at which MI is considered significant
                alpha_ci = 0.05;
                th_kl = quantile(mi_kl_standard_perm, 1-alpha_ci, 2);
                sig_kl = mi_kl_standard> th_kl;
                
                % Plot Modulograms
                subplot(2 + length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(twindow-1)*length(T_tot_sec_tot) + 2*tsignal-1  )
                imagesc(tvecsec ,pbins2,pac_cur');
                col = redbluemap; colormap(col);
                colorbar('location', 'west')
                set(gca, 'clim', [cmin cmax]);
                axis tight
                %xlim([tvecmin(1) tvecmin(end)])
                if tsignal==1
                    ylabel([' Std Wind:',num2str(pac_wind_sec_tot(1,twindow)), 's'],'fontweight','bold')
                else
                    ylabel('(rad)')
                end
                
                if twindow==1
                    title(['Signal Tot : ',num2str(T_tot_sec_tot(1,tsignal)), 's'    ])
                end
                
                % Plot MI,' || Signal Tot : ',num2str(T_tot_sec_tot(1,tsignal)), 's'
                subplot(2+ length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(twindow-1)*length(T_tot_sec_tot) + 2*tsignal); hold on
                plot(tvecsec,mi_kl_standard, 'color', 'k')
                scatter(tvecsec(find(sig_kl)), quantile(mi_kl_standard, 0.95)* ones(1,length(find(sig_kl))), 20, 'g', '*') 
                box on
                ylabel('MI_{KL}')
                ylim([0 0.2])
                
                
                if (tsignal==length(T_tot_sec_tot)) && (twindow==1)
                    title(['\color{blue}\Delta f=', num2str(delta_f_gen(1,nwidth )), ' \sigma_{s}=' , num2str(sigma_slow(1,cur_amp)), ' \sigma_{f}=' , num2str(sigma_fast(1,cur_amp))])
                end
                
                
                if (tsignal==1) && (twindow==1)
                    legend(['MI_{KL}. ', 'rp=', num2str(size(mi_kl_standard_perm,2))],['p<', num2str(alpha_ci)])
                end
                
            end
            
            %% PLot SSP
            
            
            dsslin_pac = dss_lin_pac_tot{tsignal,nwidth,cur_amp} ;
            dsslin_pac_processed = dsslin_pac.dsslin_pac_processed;
            MI_KL_mean_raw  = dsslin_pac_processed.MI_KL_mis(1,:);
            MI_KL_inf  = dsslin_pac_processed.MI_KL_mis(2,:);
            MI_KL_sup  = dsslin_pac_processed.MI_KL_mis(3,:);
            
            MI_KL_mean_filt = dsslin_pac_processed.MI_KL_mis_filt(1,:);
            MI_KL_inf_filt  = dsslin_pac_processed.MI_KL_mis_filt(2,:);
            MI_KL_sup_filt  = dsslin_pac_processed.MI_KL_mis_filt(3,:);
            
            
            % Plot Raw Modulogram
            subplot(2+ length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(length(pac_wind_sec_tot)+1-1)*length(T_tot_sec_tot) + 2*tsignal-1)
            imagesc(dsslin_pac_processed.tVect ,dsslin_pac_processed.pbins,dsslin_pac_processed.PAC_raw);
            col = redbluemap; colormap(col);
            colorbar('location', 'west')
            set(gca, 'clim', [cmin cmax]);
            axis tight
            xlim([dsslin_pac_processed.tVect(1) dsslin_pac_processed.tVect(end)])
            ylabel('(rad)')
            box on
            
            if tsignal==1
                ylabel([' DSS lin Raw'],'fontweight','bold')
            end
            
            
            % Plot Filteres Modulogram
            subplot(2+ length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(length(pac_wind_sec_tot)+2-1)*length(T_tot_sec_tot) + 2*tsignal-1)
            imagesc(dsslin_pac_processed.tVect ,dsslin_pac_processed.pbins,dsslin_pac_processed.PAC_filt);
            col = redbluemap; colormap(col);
            colorbar('location', 'west')
            set(gca, 'clim', [cmin cmax]);
            axis tight
            xlim([dsslin_pac_processed.tVect(1) dsslin_pac_processed.tVect(end)])
            ylabel('(rad)')
            xlabel('[min]')
            box on
            
            if tsignal==1
                ylabel([' DSS lin Filtered'],'fontweight','bold')
            end
            
            
            % Plot Raw MI
            subplot(2+ length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(length(pac_wind_sec_tot)+1-1)*length(T_tot_sec_tot) + 2*tsignal)
            ciplot(MI_KL_inf,MI_KL_sup,dsslin_pac_processed.tVect, [0.5 0.5 0.5]); hold on;
            p1=plot(dsslin_pac_processed.tVect,MI_KL_mean_raw,'k');
            p2=plot(dsslin_pac_processed.tVect,MI_KL_sup,'color', [0.5 0.5 0.5]);
            
            box on
            xlim([dsslin_pac_processed.tVect(1) dsslin_pac_processed.tVect(end)])
            ylim([0 0.2])
            ylabel('MI_{KL}')
            
            if tsignal==1
                legend([p1,p2], {'MI_{KL}',['rs=', num2str(dsslin_pac_processed.Kappa1),'x' , num2str(dsslin_pac_processed.Kappa2)]})
            end
            
            
            % Plot Filtered MI
            subplot(2+ length(pac_wind_sec_tot),2*length(T_tot_sec_tot), 2*(length(pac_wind_sec_tot)+2-1)*length(T_tot_sec_tot) + 2*tsignal)
            ciplot(MI_KL_inf_filt,MI_KL_sup_filt,dsslin_pac_processed.tVect, [0.5 0.5 0.5]); hold on;
            
            plot(dsslin_pac_processed.tVect,MI_KL_mean_filt,'k')
            box on
            xlim([dsslin_pac_processed.tVect(1) dsslin_pac_processed.tVect(end)])
            xlabel('[min]')
            ylim([0 0.2])
            ylabel('MI_{KL}')
            
            
            
        end

        namecur = ['non_stat_mod_comparison_dF',num2str(delta_f_gen(1,nwidth )),'Hz_sigId',num2str(cur_amp)];


        if saveFig
            set(gcf, 'position', [1          18        1887        1056])
            saveas(gcf,[lpath,namecur,'.png' ])
        end
    end
end















