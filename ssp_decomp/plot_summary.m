function plot_summary(modelOsc,crange , WOI)
%% Plots a quick summary of the oscillation decomposition model returned by SSP_DECOMP
%   modelOsc is a strucure estimated with ssp_decomp
%   crange = [cmin cmax] color range for the PSD plot
%   WOI : optionnal parameter allowing to plot the psd for a giving window
addpath(genpath('../spectral_analysis'))  % Needs Chronux toolbox to estimate multitaper PSD.

staPointId = modelOsc.startPointId;
endPointId = modelOsc.endPointId;
Nwind  = length(staPointId);
Fs     = modelOsc.Fs;

f_y=0.01:0.01:120;




[f_tot, a_tot, q_tot, R_tot, s_tot_param, s_param , Nosc_ind,NNharm, s_mt,s_mt_tot, fmt] = gather_params(modelOsc,f_y);


%% Scatter plot fitted model parameter for each window
colorCur = [0         0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0         0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

figure

% Oscillations peak frequency
subplot(2,2,1); hold on
for nosc=1:Nosc_ind
    for nosc_r=1:NNharm(1,nosc)
        curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nosc_r) ;
        scatter(staPointId/(60*Fs), f_tot(curId,:), 20,colorCur(nosc,:))
    end
end
title('Freq'); xlabel('[min]')

% Scaling parameter a in (0,1)
subplot(2,2,2); hold on
for nosc=1:Nosc_ind
    for nosc_r=1:NNharm(1,nosc)
        curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nosc_r) ;
        scatter(staPointId/(60*Fs), a_tot(curId,:), 20,colorCur(nosc,:))
    end
end
title('a'); xlabel('[min]')

% Process Noise Covariance
subplot(2,2,3); hold on
for nosc=1:Nosc_ind
    for nosc_r=1:NNharm(1,nosc)
        curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nosc_r) ;
        scatter(staPointId/(60*Fs), q_tot(curId,:), 20,colorCur(nosc,:))
    end
end
title('\sigma^2'); xlabel('[min]')

% Observation Noise
subplot(2,2,4); hold on
scatter(staPointId/(60*Fs), R_tot(1,:), 20,'k')
title('Freq'); xlabel('[min]')
title('R'); xlabel('[min]')

%% Plot multitpaer and parametric Spectra
figure

% Color scale

if nargin > 1
    cmin=crange(1,1);
    cmax=crange(1,2);
else
    cmax=60;
    cmin=-30;
end



% Parametric Spectra
p=subplot(2,1,1);
imagesc(staPointId/(60*Fs), f_y, 10*log10(s_tot_param))
colormap(jet)
axis xy;
ylabel('Freq [Hz]')
ylim([0 60])
title('Parametric PSD [dB]')
xlabel(['Time [min]'])
caxis([cmin cmax])
colorbar

% Multitaper Spectra
p2=subplot(2,1,2);
imagesc(staPointId/(60*Fs), fmt, 10*log10(s_mt))
colormap(jet)
axis xy;
ylabel('Freq [Hz]')
ylim([0 60])
title('Multitaper [dB]')
xlabel(['Time [min]'])
caxis([cmin cmax])
colorbar

linkaxes([p,p2])
ylim([0 60])


Nosc     = sum(NNharm);    % Total number of oscillation
Nosc_ind = length(NNharm);
%% Plot psd for one window of interest

if nargin > 2
    figure; hold on
    plot( fmt,  10*log10(s_mt(:,WOI)))
    plot( f_y, 10*log10(s_tot_param(:,WOI)))
    
    for nosc=1:Nosc_ind
        for nn_osc_r=1:  NNharm(nosc)
            curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r) ;
            
            plot( f_y, 10*log10(s_param (:,WOI, curId)), 'Color', colorCur( nosc,:), 'linewidth', 2)
            plot( fmt, 10*log10(s_mt_tot(:,WOI, curId)), 'Color', colorCur( nosc,:), 'linewidth', 2)
            
        end
    end
    
end


%% Plot temporal signal

sumHarmn=1; % Plot all harmonics separately or together
x_t_n_tot=zeros(2*Nosc,(endPointId(end)-staPointId(1)));
y_tot=zeros(1,(endPointId(end)-staPointId(1)));
ta=(0:endPointId(end)-staPointId(1))/Fs;

for tt=1:Nwind
    x_t_n_cur=modelOsc.res{1,tt}.x_t_n(:,2:end);
    x_t_n_tot(:, (staPointId(tt):endPointId(tt))-staPointId(1)+1)=x_t_n_cur;
    y_tot(1,(staPointId(tt):endPointId(tt))-staPointId(1)+1)=modelOsc.res{1,tt}.y';
end


figure;
p=subplot(Nosc_ind+1,1,1);hold on
plot(ta/60, y_tot, 'color','k')
plot(ta/60, sum(x_t_n_tot(1:2:end,:),1), 'color','g')

for nosc=1:Nosc_ind
    pcur=subplot(Nosc_ind+1,1,nosc+1);hold on
    
    if sumHarmn
        curId= 2*((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + 1) -1;
        curR=curId:2:curId+2*(NNharm(1,nosc)-1);
        x_r_sum1= sum(x_t_n_tot ( curR,:),1);
        x_r_sum2= sum(x_t_n_tot( curR+1,:),1);
        
        plot(ta/60, x_r_sum1, 'linewidth',2, 'color','b')
        plot(ta/60, x_r_sum2, 'linewidth',2, 'color',[0.6 0.6 1])
        
        
    else
        
        for nosc_r=1:NNharm(1,nosc)
            curId= 2*((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nosc_r) -1;
            
            plot(ta/60, x_t_n_tot(curId,:), 'linewidth',2, 'color','b')
            plot(ta/60, x_t_n_tot(curId+1,:), 'linewidth',2, 'color',[0.6 0.6 1])
            
        end
    end
    
    p=[p, pcur];
end

linkaxes(p)

end
