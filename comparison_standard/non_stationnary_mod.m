%% Use a non stationnary modulated Signal to compare Standard Analysis and SSP Lin
%%  Paths
addpath(genpath('./../'))
spath = './result_non_stationnary_simulations/';


%% Simulate signal

% Sampling Freq
Fs = 250; 

% Amplitude Modulation
Kmodmin = 0;
Kmodmax = 1;

% Modulation Phase
Phimodmin = 0;
Phimodmax = pi;

% Time vector and SignalLength
get_tt = @(t_tot_sec) (1:t_tot_sec*Fs);
T_tot_sec_1 = 20*60;
T_tot_sec_2 = 10*60;
T_tot_sec_3 = 5*60;
T_tot_sec_4 = 2*60;

T_tot_sec_tot = [T_tot_sec_1,T_tot_sec_2,T_tot_sec_3,T_tot_sec_4];
signal_tot = cell(1,length(T_tot_sec_tot));

% Slow component width
delta_f_gen = [1,3];

% Components + Noise amplitude
sigma_slow = [1, 0.5];
sigma_fast = [0.8 2];
sigma_noise = [1 1];

% Components main frequency peak
f_slow = 1;
f_fast = 10;
f_fast_band = [8,12];
Norder = floor(1.65*Fs./delta_f_gen)*2+1;

% Generate Signal for different total length
for tsignal = 1:length(T_tot_sec_tot)
    disp([num2str(tsignal), '/',num2str(length(T_tot_sec_tot))])
    
    % Get Time Vector
    signal_tot{1,tsignal}.T_tot_sec = T_tot_sec_tot(1,tsignal);
    tt = get_tt(T_tot_sec_tot(1,tsignal));
    
    % Noisy modulation
    K_cur_clean = get_K_mod  (tt,Kmodmin,Kmodmax);
    P_cur_clean = get_Phi_mod(tt,Phimodmin,Phimodmax);
    
    K_cur = max(K_cur_clean +  (Kmodmax-Kmodmin)/3*normrnd(0,1,Fs*T_tot_sec_tot(1,tsignal),1),0);
    P_cur = P_cur_clean +  (Phimodmax-Phimodmin)/3*normrnd(0,1,1,Fs*T_tot_sec_tot(1,tsignal));
    
    % Save params
    signal_tot{1,tsignal}.tt=tt;
    signal_tot{1,tsignal}.Kmod = K_cur;
    signal_tot{1,tsignal}.Kmod_clean = K_cur_clean;
    signal_tot{1,tsignal}.Phimod = P_cur;
    signal_tot{1,tsignal}.Phimod_clean = P_cur_clean;
    
        
    sig_tmp = normrnd(0,1,1,Fs*T_tot_sec_tot(1,tsignal));
    orferfilt = 400;
    
    x_slow = zeros(length(delta_f_gen),length(sigma_slow) ,Fs*T_tot_sec_tot(1,tsignal));
    x_fast = zeros(length(delta_f_gen),length(sigma_slow) ,Fs*T_tot_sec_tot(1,tsignal));
    y_tot  = zeros(length(delta_f_gen),length(sigma_slow) ,Fs*T_tot_sec_tot(1,tsignal));
    
    % Generate slow and fast component by filtering out gaussian noise
    for nwidth = 1 :length(delta_f_gen)
        slow_bandwidth = [max(f_slow - delta_f_gen(1,nwidth),0), f_slow + delta_f_gen(1,nwidth)];
        
        tt_filt = linspace(-Norder(1,nwidth)/2, Norder(1,nwidth)/2, Norder(1,nwidth));
        
        b_t_slow_1 = blackman(Norder(1,nwidth))'.*cos(2*pi*f_slow*tt_filt/Fs);
        b_t_slow_2 = blackman(Norder(1,nwidth))'.*sin(2*pi*f_slow*tt_filt/Fs);
        x_1_slow = conv(sig_tmp,b_t_slow_1,'same');
        x_2_slow = conv(sig_tmp,b_t_slow_2,'same');
        
        b_t_fast_1 = blackman(Norder(1,nwidth))'.*cos(2*pi*f_fast*tt_filt/Fs);
        X_1_fast = conv(sig_tmp,b_t_fast_1,'same');
        
        phi_slow = atan2(x_2_slow,x_1_slow);
        x_fast_tmp = X_1_fast .* (1+ K_cur'.*cos(phi_slow+P_cur));
                
        for cur_amp = 1:length(sigma_slow)
            x_slow(nwidth,cur_amp,:) = sigma_slow(1,cur_amp)* x_1_slow/std(x_1_slow);
            x_fast(nwidth,cur_amp,:) = sigma_fast(1,cur_amp)* x_fast_tmp/std(x_fast_tmp);
            y_tot (nwidth,cur_amp,:) = squeeze(x_slow(nwidth,cur_amp,:))' + squeeze(x_fast(nwidth,cur_amp,:))' + sigma_noise(1,cur_amp)*normrnd(0,1,1,Fs*T_tot_sec_tot(1,tsignal));
        end
        
    end
    
    
    signal_tot{1,tsignal}.x_slow=x_slow;
    signal_tot{1,tsignal}.x_fast=x_fast;
    signal_tot{1,tsignal}.y_tot=y_tot;
    
    signal_tot{1,tsignal}.delta_f_gen=delta_f_gen;
    signal_tot{1,tsignal}.sigma_slow=sigma_slow;
    signal_tot{1,tsignal}.sigma_fast=sigma_fast;
    signal_tot{1,tsignal}.sigma_noise=sigma_noise;
    
end

save([spath, 'non_stationnary_sim_signal_tot_final'],'-v7.3')


%% Check Signal

% figure
% p_tot=[];
% for nwidth = 1 :length(delta_f_gen)
%
%
%      for cur_amp = 1:length(sigma_slow)
%       x_s = squeeze(signal_tot{1,end}.x_slow(nwidth,cur_amp,:))   ;
%       x_f = squeeze(signal_tot{1,end}.x_fast(nwidth,cur_amp,:))   ;
%       y_t = squeeze(signal_tot{1,end}.y_tot(nwidth,cur_amp,:))   ;
%
%       p=subplot(length(delta_f_gen) ,length(sigma_slow),(nwidth-1) *length(sigma_slow)+cur_amp  ); hold on
%
%       %plot(x_s,'b')
%       %plot(x_f,'r')
%       plot(y_t,'k')
%       title(['\Delta f^{gen}_{slow} = ', num2str(delta_f_gen(1,nwidth )), ' \sigma_{s} = ' , num2str(sigma_slow(1,cur_amp)), ' \sigma_{f} = ' , num2str(sigma_fast(1,cur_amp)) ])
%
%       p_tot = [p_tot,p];
%
%      end
% end
%
% linkaxes([p_tot])



%% Perform Standard PAC Analysis

% Epoch length for processing
pac_wind_sec_tot = [6,30,60,120];
lengthWind_tot = pac_wind_sec_tot*Fs;

% Corresponding Filter Order
filter_order_tot = [100, 200,300,400];

% Fast and slow band definition
slowfreq = [0.1,1];
fastfreq = [8,12];

% Number of bins
Npb_standard = 18;
pbins2 = linspace(-pi,pi, Npb_standard);

% Permutation Analysis Parameters
doPerm = 1;
perm_range_sec = pac_wind_sec_tot/2;
Nperm = 200;

% Loop on Total Signal Length
for tsignal = 1:length(T_tot_sec_tot)    
    
    signal_tot{1,tsignal}.pac_std = cell(1,length(lengthWind_tot));
    signal_tot{1,tsignal}.pac_wind_sec_tot_std=pac_wind_sec_tot;
    ta = signal_tot{1,tsignal}.tt;
    N_tot=floor(length(ta)./lengthWind_tot);
    
    % Loop on the epoch length
    for twindow = 1:length(pac_wind_sec_tot)
        
        Nwind = N_tot(1,twindow);
        Lwind_sec = pac_wind_sec_tot(1,twindow);
        filter_order = filter_order_tot(1,twindow);
        
        pac_standard   = zeros(Nwind,Npb_standard-1,length(delta_f_gen),length(sigma_slow));
        mi_kl_standard = zeros(Nwind,length(delta_f_gen),length(sigma_slow));
        mi_l1_standard = zeros(Nwind,length(delta_f_gen),length(sigma_slow));
        
        mi_kl_standard_perm = zeros(Nwind,Nperm,length(delta_f_gen),length(sigma_slow));
        mi_l1_standard_perm = zeros(Nwind,Nperm,length(delta_f_gen),length(sigma_slow));
        
        
        % loop on slow component Spectral Width
        for nwidth = 1 :length(delta_f_gen)
            % loop on components amplitude variable
            for cur_amp = 1:length(sigma_slow)
                
                % Grasp Signal
                signal = squeeze(signal_tot{1,tsignal}.y_tot(nwidth,cur_amp,:))';
                startPointID=1+(0:Nwind-1)*Lwind_sec*Fs;
                endPointID=startPointID+Lwind_sec*Fs-1;
                
                pac_t_cur = zeros(Nwind,Npb_standard-1);
                MI_kl_t_cur = zeros(Nwind,1);
                MI_l1_t_cur = zeros(Nwind,1);
                
                tail=ceil(filter_order/2);
                endPointId_notail = endPointID-tail;
                phase_slow = zeros(1,sum(endPointId_notail-startPointID+1));
                ampli_fast = zeros(1,sum(endPointId_notail-startPointID+1));
                
                disp(['Total Length :', num2str(tsignal), '/',num2str(length(T_tot_sec_tot)), ' | Df :' ,num2str(nwidth), '/', num2str(length(delta_f_gen)), ' | sig : ', num2str(cur_amp), '/', num2str(length(sigma_slow))])
                         
                %% Estimate PAC - MI
                for tt = 1:Nwind
                    disp(['Total Length :', num2str(tsignal), '/',num2str(length(T_tot_sec_tot)), ' | Df :' ,num2str(nwidth), '/', num2str(length(delta_f_gen)), ' | sig : ', num2str(cur_amp), '/', num2str(length(sigma_slow)), ' | Wnd : ', num2str(tt), '/', num2str(Nwind)])
                
                    % Filter Signal
                    signal_cur = detrend(signal(1,startPointID(tt):endPointID(tt)));
                    [slow_tmp, tail_slow] = quickbandpass(signal_cur, Fs, slowfreq,filter_order);
                    [fast_tmp, tail_fast] = quickbandpass(signal_cur, Fs, fastfreq,filter_order);
                    
                    % Grasp Slow phase and fast amplitude
                    x_fast_hilb = hilbert(fast_tmp(1:end-tail_fast));
                    x_slow_hilb = hilbert(slow_tmp(1:end-tail_slow));
                    amp   = abs(x_fast_hilb);
                    phase = angle(x_slow_hilb);
                    
                    % Save them for permutation analysis
                    phase_slow (startPointID(tt):endPointId_notail(tt)) =phase;
                    ampli_fast (startPointID(tt):endPointId_notail(tt)) = amp;
                    
                    % Compute coupling
                    pa=phaseamp(amp,phase,pbins2);
                    pa = pa/(sum(pa)*(2*pi/(Npb_standard-1)));
                    pac_t_cur (tt,:) = pa;
                    
                    % Compute Modulation Indices
                    MI_L1_tmp = sum(abs(pa-(1/(2*pi)))*(2*pi/(Npb_standard-1)));
                    MI_KL_tmp = sum(pa.*log2(pa.*(2*pi))*(2*pi/(Npb_standard-1)));
                    
                    MI_l1_t_cur(tt,1)=MI_L1_tmp;
                    MI_kl_t_cur(tt,1)=MI_KL_tmp;
                end
                
                pac_standard(:,:,nwidth,cur_amp)=pac_t_cur;
                mi_kl_standard(:,nwidth,cur_amp) = MI_kl_t_cur;
                mi_l1_standard(:,nwidth,cur_amp) = MI_l1_t_cur;
                
                %% Do permutation test for significance
                if doPerm
                    
                    % First and last window used for permutation
                    Nwind_min = sum((1:Nwind)*Lwind_sec<perm_range_sec(1,twindow))+2;
                    Nwind_max = sum((Nwind*Lwind_sec-(1:Nwind)*Lwind_sec)>perm_range_sec(1,twindow));
                    Dt = 2*perm_range_sec(1,twindow)*(rand(1,Nperm)-1/2);
                    
                    MI_L1_perm = inf+zeros(Nwind,Nperm);
                    MI_KL_perm = inf+zeros(Nwind,Nperm);
                    
                    % loop on windows
                    for tt = Nwind_min:1:Nwind_max
                        disp(['Total Length :', num2str(tsignal), '/',num2str(length(T_tot_sec_tot)), ' | Df :' ,num2str(nwidth), '/', num2str(length(delta_f_gen)), ' | sig : ', num2str(cur_amp), '/', num2str(length(sigma_slow)), ' | Wnd : ', num2str(tt), '/', num2str(Nwind), '| Perm'])
                        phase_cur = phase_slow (startPointID(tt):endPointId_notail(tt));
                        
                        parfor perm_cur = 1:Nperm
                            
                            shifted_start = startPointID(tt)-floor(Dt(1,perm_cur)*Fs);
                            shifted_end   = endPointId_notail(tt)-floor(Dt(1,perm_cur)*Fs);
                            
                            
                            amp_perm =  ampli_fast(shifted_start:shifted_end);
                            pa = phaseamp(amp_perm,phase_cur,pbins2);
                            pa = pa/(sum(pa)*(2*pi/(Npb_standard-1)));
                            
                            MI_L1_perm_tmp = sum(abs(pa-(1/(2*pi)))*(2*pi/(Npb_standard-1)));
                            MI_KL_perm_tmp = sum(pa.*log2(pa.*(2*pi))*(2*pi/(Npb_standard-1)));
                            
                            MI_L1_perm(tt,perm_cur) = MI_L1_perm_tmp;
                            MI_KL_perm(tt,perm_cur) = MI_KL_perm_tmp;
                            
                        end
                        
                    end    
                    
                     mi_kl_standard_perm (:,:,nwidth,cur_amp) = MI_KL_perm;
                     mi_l1_standard_perm (:,:,nwidth,cur_amp) = MI_L1_perm;
                end
                
               
                
            end
            
        end
        
        signal_tot{1,tsignal}.pac_std{1,twindow}.pac_standard=pac_standard;
        signal_tot{1,tsignal}.pac_std{1,twindow}.mi_kl_standard=mi_kl_standard;
        signal_tot{1,tsignal}.pac_std{1,twindow}.mi_l1_standard=mi_l1_standard;
        
        signal_tot{1,tsignal}.pac_std{1,twindow}.mi_kl_standard_perm=mi_kl_standard_perm;
        signal_tot{1,tsignal}.pac_std{1,twindow}.mi_l1_standard_perm=mi_l1_standard_perm;

        
        
    end
end


save([spath, 'non_stationnary_sim_signal_tot_final'],'-v7.3')

%% Apply Oscillation decomposition and preliminary PAC

lengthWind_ssp = 6;
% Initialize and fit
em_its=1000;
doPlot=0;
namecur='test_fit';
init_params=struct();
init_params.f_init     =  [0.5, 8];
init_params.a_init     =  [0.98, 0.98];
init_params.sigma2_init=  [1,0.1];
init_params.R          =  10;
init_params.NNharm     =  [1, 1];
convergenceTolerance=eps;



for tsignal = 1:length(T_tot_sec_tot)
    disp([num2str(tsignal), '/',num2str(length(T_tot_sec_tot))])
    signal_tot{1,tsignal}.decomp = cell(length(delta_f_gen),length(sigma_slow));
    signal_tot{1,tsignal}.pac = cell(length(delta_f_gen),length(sigma_slow));
    
    
    
    Nwind=floor(T_tot_sec_tot(1,tsignal)/lengthWind_ssp);
    
    startPointID=1+(0:Nwind-1)*lengthWind_ssp*Fs;
    endPointID=startPointID+lengthWind_ssp*Fs-1;
    
    for nwidth = 1 :length(delta_f_gen)
        for cur_amp = 1:length(sigma_slow)
            
            disp(['Width Gen : ' , num2str(nwidth),'/', num2str(length(delta_f_gen) ), ' Amp Gen : ' , num2str(cur_amp),'/', num2str(length(sigma_slow) )])
            
            signal = squeeze(signal_tot{1,tsignal}.y_tot(nwidth,cur_amp,:))';
            modelOsc=ssp_decomp(signal',Fs,startPointID,endPointID,em_its,convergenceTolerance, init_params,namecur,doPlot);
            
            signal_tot{1,tsignal}.decomp{nwidth,cur_amp}=modelOsc;
            
            
            slowID= 1; fastID=2;
            Kappa = [0;0];
            Nlinconst = 100;Nlag =inf;
            %[latent_pac, SS_tot] =ssp_pac(modelOsc, Kappa, slowID, fastID, Nlag,Nlinconst);
            [latent_pac, SS_tot] =ssp_pac(modelOsc, Kappa, slowID, fastID, Nlag,Nlinconst);
            
            signal_tot{1,tsignal}.pac{nwidth,cur_amp} = latent_pac;
            signal_tot{1,tsignal}.SS_tot{nwidth,cur_amp} = SS_tot;
        end
    end
    
end



save([spath, 'non_stationnary_sim_signal_tot_final_oscdecomp'],'-v7.3')

%% Apply PAC with resampling
for tsignal = 1:length(T_tot_sec_tot)
    disp([num2str(tsignal), '/',num2str(length(T_tot_sec_tot))])
    
    for nwidth = 1 :length(delta_f_gen)
        for cur_amp = 1:length(sigma_slow)
            
            disp(['Width Gen : ' , num2str(nwidth),'/', num2str(length(delta_f_gen) ), ' Amp Gen : ' , num2str(cur_amp),'/', num2str(length(sigma_slow) )])
            
            modelOsc = signal_tot{1,tsignal}.decomp{nwidth,cur_amp};
            slowID= 1; fastID=2;
            Kappa = [200;200];
            Nlinconst = 100; Nlag =inf;
            [latent_pac, SS_tot] =ssp_pac(modelOsc, Kappa, slowID, fastID, Nlag,Nlinconst);
            
            name_cur = ['ns_simsig_paclin_resamp_ttot', num2str(T_tot_sec_tot(tsignal)), 's_dF', num2str(delta_f_gen(nwidth)), '_ampId', num2str(cur_amp),'2' ];
            save([spath, name_cur], 'latent_pac', 'SS_tot' , '-v7.3')
        end
    end
    
end

%% Process SSP with Resampling

ssp_lin_pac_tot = cell(length(T_tot_sec_tot) , length(delta_f_gen),length(sigma_slow) );
for tsignal = 1:length(T_tot_sec_tot)
    disp([num2str(tsignal), '/',num2str(length(T_tot_sec_tot))])
    
    for nwidth = 1 :length(delta_f_gen)
        for cur_amp = 1:length(sigma_slow)
            name_cur = ['ns_simsig_paclin_resamp_ttot', num2str(T_tot_sec_tot(tsignal)), 's_dF', num2str(delta_f_gen(nwidth)), '_ampId', num2str(cur_amp) ];
            load([spath, name_cur])
            
            ssplin_pac = struct();
            
            ssplin_pac.T_tot_sec   = T_tot_sec_tot(1,tsignal);
            ssplin_pac.delta_f_gen = delta_f_gen(1,nwidth);
            ssplin_pac.sigma_slow  = sigma_slow(1,cur_amp);
            ssplin_pac.sigma_fast  = sigma_fast(1,cur_amp);
            ssplin_pac.sigma_noise = sigma_noise(1,cur_amp);
            
            lengthWind_ssp = 6;
            Nwind=floor(T_tot_sec_tot(1,tsignal)/lengthWind_ssp);
            
            startPointId=1+(0:Nwind-1)*lengthWind_ssp*Fs;
            endPointId=startPointId+lengthWind_ssp*Fs-1;
            
            tVect = 0.5 * (startPointId + endPointId)/Fs;
            
            doArFiltering = 1;
            Nlinconst = 100; name = ''; Nlag = Inf;
            
            ssplin_pac_processed = ssp_pac_process(latent_pac,SS_tot,doArFiltering,tVect,Nlinconst,Nlag,name);
            ssplin_pac.ssplin_pac_processed = ssplin_pac_processed;
            
            ssp_lin_pac_tot{tsignal,nwidth,cur_amp} = ssplin_pac;
            
            
        end
    end
    
end

save([spath, 'non_stationnary_sim_signal_tot_final_oscdecomp_ssp_decomp'], 'ssp_lin_pac_tot', '-v7.3')


%% Side functions

function Kmod = get_K_mod(tt,Kmodmin,Kmodmax)

T_tot = length(tt);
t1 = floor(T_tot/4);
t2 = floor(T_tot/2);
t3 = floor(3*T_tot/4);
t4 = T_tot;

Kmod=zeros(T_tot,1);
Kmod(1:t1) = Kmodmin + (Kmodmax-Kmodmin)*((1:t1)-1)/(t1-1);
Kmod((t1+1):t2) = Kmodmax + (Kmodmin-Kmodmax)*(((t1+1):t2)-t1)/(t2-t1);
Kmod((t2+1):(2*t2)) = Kmod(1:t2);


end

function Phimod = get_Phi_mod(tt,Phimodmin,Phimodmax)


T_tot = length(tt);
t2 = floor(T_tot/2);

lambda = 20/t2;
Phimod = Phimodmin+ (Phimodmax-Phimodmin) * ( 1./(1 +  exp(-lambda*(tt- t2)  )));

end

