function [f_tot, a_tot, sigma2_tot,R_hat, contrib_tot]=get_init_osc(signal,startPointID,endPointID,Fs,TW,Ktapers,Nosc_max, doRedress,doPlot)
%% GET_INIT_OSC returns initialization parameters for oscillation decomposition (according to Matsuda and al. (2017) model)
%
%
%   inputs :    signal (total (n windows) or just 1 window) (nT x 1)
%               start/end PointID : ID of the windows starting and ending points (only if full signal, else startPointID= 1, endPointID = length(signal))
%               TW : timewidth product, used to compute Multitaper spemctrogra
%               Ktapers : number of tapers used for Multitaper spemctrogra
%               Nosc_max : number of oscillation to be fitted (eg : 6, not nececessarily the number used for decomposition)
%               doRedress : 0 or 1, redress the psd before fitting the oscillations
%               freq,    the frequency point at which psd_raw is calculated
%               Fs,      sampling fre
%               doPlot : plot a summary of the fit
%
%
%  outpurs :    f_tot       : Center frequencies  of the oscillations
%               a_tot       : In ]0,1[, amplitude of the oscillations
%               sigma2_tot  : Noise covariances   of the oscillations
%               R_hat       : Estimated observation noise
%               contrib_tot : Euristic estimation of the contribution of each oscillation. 
%                   (Can be used to sort oscillation and keep only the most
%                   important in the decomposition.)
%
%
%   The model used to redress the PSD is described in : Haller M, Donoghue T ... and Voytek B (2018)
%   Chronux toolbox is needed here
%   Does not work to initialize oscillations harmonics YET.


% Model used to redress psd (in dB)
psd_deviation_model = @ (b,k,X,F) b - log(1 + (F./k).^X );

% Signal parameters
Twind = mean(endPointID-startPointID)/Fs; % Window lengh in sec
Nwind = length(startPointID);             % Window number
fRes  = 2*TW/Twind;                       % Mt spectral resolution

% Initialize outputs
f_tot       = zeros(Nwind,Nosc_max);
a_tot       = zeros(Nwind,Nosc_max);
sigma2_tot  = zeros(Nwind,Nosc_max);
R_hat         = zeros(Nwind,1);
contrib_tot   = zeros(Nwind,Nosc_max);

% Multitaper parameters
specparams.Fs     = Fs;
specparams.tapers = [TW Ktapers]; 
specparams.fpass  = [0 100];

for windId = 1:Nwind
    %disp(['Wind :', num2str(windId), '/', num2str(Nwind)])
    
    % Grasp window #windId
    curSignal  = signal(startPointID(windId):endPointID(windId));
    [psd_raw,freq]=mtspectrumc(curSignal,specparams);
    psd_raw_dB=10*log10(psd_raw);
    
    % Estimate R_hat : observation noise
    rangeEnd = ceil(0.99*length(psd_raw_dB));
    R_hat_dB = mean(psd_raw_dB(rangeEnd:end));
    R_hat_t    = Fs*10^(R_hat_dB/10);
    
    % Remove PSD 1/f^X component according to 'psd_deviation_model' - Haller M, Donoghue T ... and Voytek B (2018)
    if doRedress
        initial_guess_psd_param = [psd_raw_dB(1),80,2];
        [redressedPSD,redressed_param,redressed_param_0] =resdressPSD(freq, psd_raw_dB, psd_deviation_model, initial_guess_psd_param);
        psd=redressedPSD;
    else
        psd = psd_raw_dB;
    end
    
    % Fit Nosc_max oscillation based on psd estimate
    initial_guess_osc_param = [0.99,0.1]; % initial guess common to each osccillation [a sigma2]
    [f_tot_t, a_tot_t, sigma2_tot_t] = get_osc_init_from_psd(freq, psd, initial_guess_osc_param, Nosc_max,Fs,fRes);
    
    % if the psd has been redressed before fit, remove psd_deviation_model(f_i) contribution
    if doRedress
        for nosc_cur = 1 : length(f_tot_t)
            [~,loc_f_nosc] = min(abs(freq-f_tot_t(1,nosc_cur)));
           
            b_fit = redressed_param(1);
            k_fit = redressed_param(2);
            X_fit = redressed_param(3);
            dev = psd_deviation_model(b_fit,k_fit,X_fit,freq);
            dev_nosc = dev(loc_f_nosc);
            
            sigma2_tot_nosc_prev = sigma2_tot_t(1,nosc_cur);
            sigma2_tot_t(1,nosc_cur) = sigma2_tot_nosc_prev * 10^(dev_nosc/10);
        end
    end
    
    % Get estimated oscillations parametric psd (See Matsuda (2017))
    [H_res, H_ires] = get_theoretical_psd(freq,Fs,f_tot_t,a_tot_t,sigma2_tot_t);

    
    % Estimate contribution of each oscillation to the psd (in ~ dB.Hz)
    for nosc = 1 : Nosc_max
        not_osc = find( (1:Nosc_max) ~=nosc );
        psd_no_nosc_dB= sum(H_ires(not_osc,:),1);
        contrib_tot(windId,nosc) = 10*log10(sum(abs(psd_no_nosc_dB-H_res )));
    end
   
    f_tot      (windId,:) = f_tot_t;
    a_tot      (windId,:) = a_tot_t;
    sigma2_tot (windId,:) = sigma2_tot_t;
    R_hat        (windId,1) = R_hat_t;
    

    % Plots the fitted spectrum
    % Plot fit steps
    if doPlot && (Nwind <= 3)
        figure;  
        % First pass fit
        subplot(4,2,1);hold on
        plot(freq, psd_raw_dB, 'k')
        if doRedress
            plot(freq, psd_deviation_model(redressed_param_0(1), redressed_param_0(2),redressed_param_0(3), freq),'r')
        else 
            plot(freq, 0*ones(1,length(freq)),'r')
        end
        xlim([0 freq(end)])
        title('1st pass fit')
        
        % Highlights non oscillatory components
        subplot(4,2,3);hold on
        plot(freq, psd_raw_dB,'k')
        if doRedress
            residuals = psd_deviation_model(redressed_param_0(1), redressed_param_0(2),redressed_param_0(3), freq)-psd_raw_dB';
            threshold = quantile(residuals, 0.2);
            keptId = find(residuals> threshold);
            plot(freq(keptId), psd_raw_dB(keptId),'b')
        end
        xlim([0 freq(end)])
        title('Non oscillatory components')
        
        % Second pass fit
        subplot(4,2,2);hold on
        plot(freq, psd_raw_dB,'k')
        if doRedress
            plot(freq,  psd_deviation_model(redressed_param(1), redressed_param(2),redressed_param(3), freq),'r')
        else 
            plot(freq, 0*ones(1,length(freq)),'r')
        end
        xlim([0 freq(end)]);
        title('2nd pass fit')
        
        % Plot redressed psd
        subplot(4,2,4);hold on
        plot(freq, psd, 'b')
        title('Redressed psd')
        xlim([0 freq(end)])
        
        % Plot oscillations and preliminary fit
        subplot(2,1,2); hold on
        plot(freq, psd_raw_dB,'color', 'k', 'linewidth',1.5 )
        plot(freq, 10*log10(H_res  +  R_hat_t/Fs) ,'color', 'b', 'linewidth',1.5 )
        plot(freq, 10*log10(H_ires  +  R_hat_t/Fs))
        
        if doRedress
            plot(freq,dev,'color', 'r', 'linewidth',1.5 )
        end
        title(['Fitted oscillations window #',num2str(windId), '/' ,num2str(Nwind)])
        
    end  
end

if doPlot && (Nwind > 3)
   figure;
   % Oscillation amplitude
   subplot(2,2,1); hold on
   for nosc =1:Nosc_max
       scatter(startPointID/Fs, a_tot(:,nosc))
   end
   title('a_{tot}')
   % Oscillation frequency
   subplot(2,2,3); hold on
   for nosc =1:Nosc_max
       scatter(startPointID/Fs, f_tot(:,nosc))
   end
   title('f_{tot}')
   % Oscillation contribution
   subplot(2,2,2); hold on
   for nosc =1:Nosc_max
       scatter(startPointID/Fs, contrib_tot(:,nosc))
   end
   title('Oscillation contribution') 
   % Oscillation process noise
   subplot(2,2,4); hold on
   for nosc =1:Nosc_max
       scatter(startPointID/Fs, sigma2_tot(:,nosc))
   end
   scatter(startPointID/Fs,R_hat,'k')
   title('sigma^2_{tot}')
end

end



%% Redress the power spectral density by fitting : h(freq) =  b - log (1 + (freq/k)^X )
%  Original idea from FOOOF algorithm : Haller M, Donoghue T ...Voytek B (2018)
function [redressedPSD, param_hat,param_hat_0] = resdressPSD(freq, raw_psd, deviation_model, initial_guess)

% Distances use for 1st and 2nd pass fits
dist     = @ (param, target, f_target)sum( abs(deviation_model(param(1),param(2), param(3), f_target) - target')  ,2);
dist0    = @ (param, target, f_target)sum( abs(deviation_model(initial_guess(1),param, initial_guess(3), f_target) - target')  ,2);

% First Pass fit : We fix b_0 and X_0 and only fit k
target_0    = raw_psd;
freq_0      = freq;
dist_0      = @ (param) dist0(param, target_0, freq_0);
options     = optimoptions('fmincon','Display','off'); lb=[10]; ub =[200];
param_hat_0 = fmincon(dist_0,initial_guess(2),[],[],[],[],lb ,ub, [],options);
param_hat_0 = [initial_guess(1),param_hat_0,initial_guess(3)];

% Residual are thresholded to remove PSD components which are believed to come from actual oscillations
residuals = deviation_model(param_hat_0(1), param_hat_0(2),param_hat_0(3), freq)-target_0';
threshold = quantile(residuals, 0.2);
keptId = find(residuals> threshold);

% Second Pass Fit : We fit the PSD's points associated with non oscillatory components
target_1  = raw_psd(keptId);
freq_1    = freq(keptId);
dist_1    = @ (param) dist(param, target_1, freq_1);
options   = optimoptions('fmincon','Display','off'); lb=[-10^5 10 0]; ub =[10^5 200 20];
param_hat_1 = fmincon(dist_1,initial_guess,[],[],[],[],lb ,ub, [],options);

% We remove non oscillatory fit from PSD and save deviation parameters
fitDeviation = deviation_model(param_hat_1(1), param_hat_1(2),param_hat_1(3), freq)';
redressedPSD = raw_psd - fitDeviation;
param_hat    = param_hat_1;

end



%% From (Redressed) PSD, remove Nosc oscillations according to ARMA(2,1) theoretical spectra
function [f_nosc_t, a_nosc_t, sigma2_nosc_t] = get_osc_init_from_psd(freq, psd_redressed, osc_param_guess,Nosc,Fs,fRes)

Fsf= 1/mean(diff(freq));
width_th=0.5;                 % Minimal width of a psd peak          .

% Initialization of the oscillations parameters (equirepartion of f_i on [0 100] and smaal a_i, sig2_i) 
osc_param_init=[0; 0.1;0.001].*ones(3,Nosc)+[1;0;0].*ones(3,Nosc).*(1:Nosc)*100/Nosc;

% Grasp and plot Current psd
psd_redressed_nosc = psd_redressed;

% From current psd, fit the strongest oscillation, then remove it
for nosc=1:Nosc
    
    % Find the prohiminent peaks in the psd
    % (Possible conflict with chronux toolbox 'findpeaks'). Be careful
    [psd_peaks_0,loc_psd_peaks_0] = findpeaks(psd_redressed_nosc, Fsf, 'MinPeakDistance',fRes);
    [psd_peaks,loc_psd_peaks ]    = findpeaks(psd_redressed_nosc, Fsf, 'MinPeakDistance',fRes,'MinPeakWidth',width_th);
    
    % Add the first peak if necessary
    peak_0 = psd_peaks_0(1);
    locp_0 = loc_psd_peaks_0(1);
    if psd_redressed_nosc(1)>peak_0
        peak_0=psd_redressed_nosc(1);
        locp_0 =1/Fsf;
    end
    
    if not(ismember(locp_0,loc_psd_peaks))
        psd_peaks=[peak_0;psd_peaks];
        loc_psd_peaks=[locp_0;loc_psd_peaks];
    end
    
    % Choose the biggest peak which is not too close from a previously fittted oscillation
    [~,loc_sorted_nosc] = sort(psd_peaks,'descend');
    osc_id_tmp=1;
    while (nosc>1) && (osc_id_tmp<= length(psd_peaks)) && (min(abs(loc_psd_peaks(loc_sorted_nosc(osc_id_tmp))-osc_param_init(1,1:(nosc-1))))<max(fRes,1)) 
        osc_id_tmp=osc_id_tmp+1;
    end
    
    % Break if no spectral peak
    if osc_id_tmp> length(psd_peaks)
        break
    end
   
    f_nosc_cur = max(loc_psd_peaks(loc_sorted_nosc(osc_id_tmp)), 0.1);
    
    % parametric psd of ARMA(2,1) which will be fitted to the current psd
    osc_psd  = @ (osc_param) max(10*log10(get_theoretical_psd(freq,Fs,osc_param(3),osc_param(1),osc_param(2))),0);
    
    % Width / Interval arround which we fit an oscillation to the current psd
    widthGen = 2*fRes; isInFitInterval = abs(f_nosc_cur -freq)<widthGen;
    widthFit = 5*fRes; RemoFitInterval = abs(f_nosc_cur -freq)<widthFit;

    % Distance used for the fit
    %dist_osc = @ (osc_param) sum(isInFitInterval.*abs(osc_psd(osc_param) - psd_redressed').^2);
    dist_osc = @ (osc_param) sum(isInFitInterval.*abs(osc_psd(osc_param) - psd_redressed_nosc').^2);

    % Constrained minimization over [a , sigma2]
    lb=[0.6,0.1,max(f_nosc_cur-fRes,0)]; ub=[0.99999, 2000, f_nosc_cur+fRes];
    options = optimoptions('fmincon','Display','off');
    [osc_param_hat,~]=fmincon(dist_osc,[osc_param_guess,f_nosc_cur] ,[],[],[],[], lb,ub, [],options);
    
    % Store results
    %osc_param_init  (1,nosc)  = f_nosc_cur;
    osc_param_init   (1,nosc)  =osc_param_hat(1,3);
    osc_param_init(2 :3,nosc) = osc_param_hat(1,1:2);
    
    % Remove fitted oscillation from the current psd
    removeFromPSD = osc_psd(osc_param_hat).*RemoFitInterval;
    
    % Plot fitted oscillation to check step by step if necessary
    livePlot=0;
    if livePlot
        figure;hold on
        plot(freq, osc_psd(osc_param_hat))
        xlim([0 , freq(end)])
        plot(freq,psd_redressed_nosc , 'k' )
        drawnow()
        pause()
    end
    
    psd_redressed_nosc = psd_redressed_nosc- removeFromPSD';
    
end

% Sort oscillation by frequency
[f_nosc_t, oscID] = sort(osc_param_init(1,:));
a_nosc_t          = osc_param_init(2,oscID);
sigma2_nosc_t     = osc_param_init(3,oscID);

end


