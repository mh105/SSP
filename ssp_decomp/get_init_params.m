function [Phi_init, Q_init, R_init, mu_init, S_init,NNharm,f_tot_curWind] = get_init_params(init_params,Fs,windId, Nosc_max)
%% GET_INIT_PARAMS build parameters to initialize EM Algorithm for subwindow windId
%   Nosc_indp number of independant oscilaltions
%   init_params is a structure containing :
%                        - f_init       : (Nwind_totxNosc_indp) or (1xNosc_indp)  initialisation freq centers
%                        - a_init       : (Nwind_totxNosc_indp) or (1xNosc_indp)   initialisation freq amplitudes
%                        - sigma2_init  : (Nwind_totxNosc_indp) or (1xNosc_indp)   initialisation process noise
%                        - R_init       : (Nwind_totxNosc_indp) or (1xNosc_indp)   initialisation observation noise
%                        - NNharm       : (Nwind_totxNosc_indp) or (1xNosc_indp)   harmoonics for each osc
%
%   Fs     : sampling frequency
%   windId : sub-Window which is going to be processed

% If windId is not given or if we only have 1 initialization value for each
% oscillatiom, we initialize all subwindows with the same parameters
if nargin <3 || size(init_params.f_init,1)==1
    windId=1;
end

% If Nosc_max not given, we take all the oscillations
if nargin <4 
    Nosc_max = size(init_params.f_init,2);
    orderOsc = 1:Nosc_max;
else
    [~,orderOsc_tmp] = sort(init_params.contrib(windId,:),'descend');
    orderOsc = orderOsc_tmp(1:Nosc_max);

end

a_tot_curWind      = init_params.a_init       (windId,orderOsc);
f_tot_curWind      = init_params.f_init       (windId,orderOsc);
sigma2_tot_curWind = init_params.sigma2_init  (windId,orderOsc);

NNharm = init_params.NNharm (windId,orderOsc); % Array containing the number of harmonics for each independant oscillation
Nosc_ind=length(f_tot_curWind);                % Number of independant oscillation
Nosc    = sum(NNharm );                        % Total number of oscillation
 
assert(Nosc_ind==length(a_tot_curWind));
assert(Nosc_ind==length(sigma2_tot_curWind));

R_init  = init_params.R(windId,1);     % Observation noise
mu_init = zeros(2*Nosc,1);   % Initial state means
S_init  = 3*eye(2*Nosc);     % Initial State covariances

Phi_init= zeros(2*Nosc);     % State Transition matrix
Q_init  = zeros(2*Nosc);     % State/process noise

for n_osc_ind= 1:Nosc_ind
    w_n = (1:NNharm (1,n_osc_ind)) * f_tot_curWind(1,n_osc_ind) * 2 * pi / Fs;
    a_n      = a_tot_curWind(1, n_osc_ind)+ zeros(1, NNharm (1,n_osc_ind));      % We initialize all harmonics with same amplitude
    sigma2_n = sigma2_tot_curWind(1, n_osc_ind); % Idem
    
    for nn_osc= 1 : NNharm (1,n_osc_ind)
        curId= 2*((sum(NNharm (1,1:n_osc_ind)) - NNharm (1,n_osc_ind)) + nn_osc) -1;
        Phi_tmp=get_rot_mat(a_n(1,nn_osc),w_n(1,nn_osc)); 
        Phi_init(curId:curId+1, curId:curId+1)= Phi_tmp;
        Q_init (curId:curId+1, curId:curId+1)= eye(2)*sigma2_n;
        
    end
 
end

end