function modelOsc=ssp_decomp(signal,Fs,startPointID,endPointID,em_its,convergenceTolerance, init_params,fname,doPlot, doInitAuto)
%% SSP_DECOMP computes modelOSC : Oscillation decomposition of an EEG signal.
%                                      Adapted from Matsuda and al. 2017 and extended to harmonic estimation
%
%                  INPUTS   :          - signal : (nT) x 1 n windows of length T
%                                      - Fs     : sampling frequency [Hz]
%                                      - statPointId : starting point of each time window
%                                      - endPointId  :   endind point of each time window
%                                      - em_its : number of iteration of the EM algorithm
%                                      - init_params is a structure used to initialize the EM.
%                                            if doInitAuto = 1, it must contain :
%                                                  -f_init : osc frequencies [Hz]  1 x #indpt oscillations
%                                                  -a_init : osc scaling param     1 x #indpt oscillations
%                                                  -sigma2_init osc noise          1 x #indpt oscillations
%                                                  -R : observation noise covariance
%                                                  -NNharm : 1 x #indpt oscillations an array containing the number of harmonics for each independant oscillation
%                                            else, NNharm is sufficient and we use get_init_osc
%                                      - Convergence tolerence : convergence threshold to stop EM
%                                      - fname : a string with your filename
%                                      - doPlot : 0 no plot of EM evolution 
%                                                 1 plot without live parametric psd.       (not compatible with parfor)
%                                                 2 plot with    live parametric psd.       (not compatible with parfor)        
%                                                 3 plot with    disp iteration + last plot (not compatible with parfor)
%
%
%  OUTPUTS : modelOsc, a structure containing :
%                   - the main EM parameters
%                   - res : 1 cell per time window with :
%                           - ta : time vector (in second)
%                           - y  : initial eeg window
%                           - x_t_n : the Kalman filtered estimate of the hiddden components
%                           - P_t_n and P_t_tmin1 the associated covariances 
%                           - ll : final log likelihood of the EM estimates
%                           - model_prams : a struct containing the oscillation parameters (see Matsuda 2017)
%
%  In all what follows : f denotes a frequency in Hz and w an associated radial frequency

modelOsc=struct();
modelOsc.Fs=Fs;
modelOsc.em_its=em_its;
modelOsc.convergenceTolerance=convergenceTolerance;
modelOsc.name=fname;
modelOsc.startPointId =startPointID;
modelOsc.endPointId =endPointID;
modelOsc.res=cell(1, length(startPointID));
modelOsc.init_params=init_params;
res=cell(1, length(startPointID));

if nargin < 10
    doInitAuto = 0;
end

if nargin < 9
    doPlot =0;
end

% We use parallel computing when no plot is required
if doPlot > 0
    parforArg = 0;
else
    parforArg = Inf;
end


parfor (windID=1:length(startPointID),parforArg )
%for windID=1:length(startPointID) 
    if length(startPointID)>2
        disp(['Decomposition window: ', num2str(windID),'/', num2str(length(startPointID))])
    end
    % Isolate the tt^th eeg window
    y  = signal(startPointID(windID):endPointID(windID)); y=y-mean(y);
    ta =       (startPointID(windID):endPointID(windID))/Fs;
    

    if doInitAuto % automatically initialize EM parameters
        doRedress = 1;
        Twind = length(ta)/Fs; fRes= 1;
        TW     = (fRes/2)*Twind;
        Ktapers_max = floor(2*TW)-1;
        Ktapers = max(floor(Ktapers_max*0.9),1);
        Nosc_max=6;
        Nosc_kept=sum(init_params.NNharm);
        [f_init_tmp, a_init_tmp, sigma2_init_tmp,R_init, contrib_tot]=get_init_osc(y,1,length(y),Fs,TW,Ktapers,Nosc_max, doRedress,0);

        [~,OscOrder] = sort(contrib_tot,2, 'descend');
        f_tot_ordered      = zeros(size(f_init_tmp));
        a_tot_ordered      = zeros(size(f_init_tmp));
        sigma2_tot_ordered = zeros(size(f_init_tmp));

        f_tot_ordered      (1 , :) = f_init_tmp     (1, OscOrder(1, :));
        a_tot_ordered      (1 , :) = a_init_tmp     (1, OscOrder(1, :));
        sigma2_tot_ordered (1 , :) = sigma2_init_tmp(1, OscOrder(1, :));
 
        f_init_auto      = f_tot_ordered     (:, 1:Nosc_kept);
        a_init_auto      = a_tot_ordered     (:, 1:Nosc_kept);
        sigma2_init_auto = sigma2_tot_ordered(:, 1:Nosc_kept);
        
        init_params_auto=struct();
        init_params_auto.f_init     =  f_init_auto;
        init_params_auto.a_init     =  a_init_auto;
        init_params_auto.sigma2_init=  sigma2_init_auto;
        init_params_auto.R          =  R_init;
        init_params_auto.NNharm     =  ones(1,Nosc_kept);
        [Phi_init, Q_init, R_init, mu_init, sigma_init,NNharm,f_init] = get_init_params(init_params_auto,Fs,windID);
        
    else % Use maual initialisation
        [Phi_init, Q_init, R_init, mu_init, sigma_init,NNharm,f_init] = get_init_params(init_params,Fs,windID);
        
    end
    
    Nosc_ind = length(NNharm); % Number of independant oscillations
    Nosc     = sum(NNharm);    % Total number of oscillations
    
    Q_cur    =Q_init;          % State noise
    Phi_cur  =Phi_init;        % State transition matrix
    R_cur    =R_init;          % observation noise
    mu_cur   =mu_init;         % Initial State means
    sigma_cur=sigma_init;      % Initial State covariances
    
    % Initialize plots
    if doPlot>0;figure;
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
    end
    
    % Disp EM iteration
    deltaCur=1; deltaDisp=floor(em_its/10);
    
    % Take init param as estimate if no iteration
    if em_its==0
        Phi_opt=Phi_cur;Q_opt=Q_cur;R_opt=R_cur;mu_opt=mu_cur;
        [x_t_n, P_t_n,P_t_tmin1_n, ll] = kalman_filter(y,Phi_opt,Q_opt,R_opt,mu_opt,sigma_init);
    else
        
        % Initialize em_parameters innovations
        ll_tot= zeros(1     ,em_its);
        f_tot = zeros(Nosc  ,em_its);
        a_tot = zeros(Nosc  ,em_its);
        Q_tot = zeros(2*Nosc,em_its);
        R_tot = zeros(1     ,em_its);
        
        % Boundaries for the harmonic fundamental for each independant
        % oscillation (used only if NNharm(1,i)>1)
        w_percHarm=1; 
        w_fond_range      = zeros(2, length(f_init));
        w_fond_range(1,:) = max(0,(1-w_percHarm)* f_init * 2 * pi /Fs);
        w_fond_range(2,:) = (1+w_percHarm)* f_init * 2 * pi /Fs;
        
        % Launch em_its EM iteration
        for ite=1:em_its
            
            % E-Step 
            [x_t_n_opt, P_t_n,P_t_tmin1_n, ll] = kalman_filter(y,Phi_cur,Q_cur,R_cur,mu_cur,sigma_cur);

            % M-Step
            [Phi_opt, Q_opt,R_opt, mu_opt ]=m_step(y, NNharm, x_t_n_opt, P_t_n,P_t_tmin1_n,w_fond_range);
            
            % Innovations
            ll_tot(1,ite)=ll;
            for nosc=1:Nosc
                PhiCurr=Phi_opt(2*(nosc-1)+1:2*nosc,2*(nosc-1)+1:2*nosc);
                [a_cur,w_cur]=get_rotmat_pam(PhiCurr);
                f_tot(nosc,ite)=w_cur*Fs/(2*pi);
                a_tot(nosc,ite)=a_cur;
            end
            R_tot(1,ite)=R_opt;
            Q_tot(:, ite)= diag(Q_opt);
            
            % Update Parameter
            Phi_cur=Phi_opt;
            Q_cur  =Q_opt;
            R_cur  =R_opt;
            mu_cur =mu_opt;
            
            % Break EM loop if convergence tolerence is reached
            if ite>2 && abs(ll_tot(1,ite)-ll_tot(1,ite-1))<convergenceTolerance; break
            end
            
            % Plots
            if deltaCur>=deltaDisp || (ite==em_its)
                deltaCur=0;
                
                if doPlot ==3
                    disp(['ite : ' ,num2str(ite), '/', num2str(em_its)])
                end
                
                if (doPlot ==1 ||doPlot ==2)  || (doPlot==3 && ite==em_its)
                    if doPlot>1
                        f_y=0.1:0.01:100;
                        [H_tot,H_i]=get_theoretical_psd(f_y,Fs,f_tot(:,ite)',a_tot(:,ite)',Q_tot(1:2:end,ite)');
                        subplot(4,2,3);cla; hold on
                        plot(f_y,10*log10(H_tot+R_tot(1,ite)/Fs) ,'k')
                        for nosc=1:Nosc_ind
                            for nn_osc_r=1:NNharm(1,nosc)
                                curId= (sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r;
                                plot(f_y, 10*log10(H_i(curId,:)) , 'color', colorCur(nosc,:));
                            end
                        end
                        xlabel('[Hz]'); ylabel('dB')
                        xlim([0 50])
                        
                        subplot(4,2,1);cla; hold on
                        
                    else
                        subplot(2,2,1);cla; hold on
                    end
                    
                        % Plot raw
                        plot(ta,y, 'color', 'k', 'linewidth', 1.2)
                        % Plot fit
                        plot(ta,sum(x_t_n_opt(1:2:end,2:end),1), 'color', 'g', 'linewidth', 1)
                        % Plot decomposition
                        for nosc=1:Nosc_ind
                            for nn_osc_r=1:NNharm(1,nosc)
                                curId= 2*((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r) -1;
                                plot(ta, x_t_n_opt(curId,2:end) , 'color', colorCur(nosc,:));
                            end
                        end
                        xlabel('[sec]'); ylabel('a.u')
                        
                    
                    
                    subplot(2,2,2);cla; hold on
                    scatter(1:ite, ll_tot(1,1:ite))
                    xlabel(['iteration']); ylabel('ll')
 
                    subplot(2,2,3);cla; hold on
                    for nosc=1:Nosc_ind
                        for nn_osc_r=1:NNharm(1,nosc)
                            curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r) ;
                            scatter(1:ite, abs(f_tot(curId,1:ite)) , 30, colorCur(nosc,:));
                        end
                        
                    end
                    xlabel(['iteration']);ylabel('freq [Hz]')
                    
                    subplot(4,2,6);cla; hold on
                    for nosc=1:1:Nosc_ind
                        for nn_osc_r=1:NNharm(1,nosc)
                            curId= ((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r) ;
                            scatter(1:ite,  a_tot(curId,1:ite) , 30, colorCur(nosc,:));
                        end
                    end
                    xlabel(['iteration']);ylabel('amp')
                    
                    subplot(4,2,8);cla; hold on
                    scatter(1:ite, R_tot(1,1:ite), 30, colorCur(nosc,:),  'm')
                    
                    for nosc=1:1:Nosc_ind
                        for nn_osc_r=1:NNharm(1,nosc)
                            curId= 2*((sum(NNharm(1,1:nosc)) - NNharm(1,nosc)) + nn_osc_r) -1;
                            scatter(1:ite,  Q_tot(curId,1:ite) , 30, colorCur(nosc,:));
                        end
                    end
                    xlabel(['iteration']);ylabel('noise')
                    drawnow()
                end
                
            else;deltaCur=deltaCur+1;
            end
            
        end
        
        % Last KF with final em_parameters
        [x_t_n, P_t_n,P_t_tmin1_n, ll] = kalman_filter(y,Phi_opt,Q_opt,R_opt,mu_opt,sigma_init);
        
        
    end
    
    % Save Results for this window
    model_params=struct();
    model_params.Phi  =Phi_opt;
    model_params.Q    =Q_opt;
    model_params.R    =R_opt;
    model_params.mu   =mu_opt;
    model_params.sigma=sigma_cur;
    
    res{1,windID}.ta=ta;
    res{1,windID}.y=y;
    res{1,windID}.x_t_n=x_t_n;
    res{1,windID}.P_t_tmin1_n=P_t_tmin1_n;
    res{1,windID}.P_t_n=P_t_n;
    res{1,windID}.ll=ll;
    res{1,windID}.model_prams=model_params;
    
    
end

% Gather Result in modelOsc
modelOsc.res=res;

end




