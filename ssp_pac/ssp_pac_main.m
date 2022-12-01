function [processed_ssp_pac,regressed_pac, SS_tot, X_t_tot] = ssp_pac_main(modelOsc, slowID, fastID,Kappa_tot, Nlag,Nconst,doArFiltering)
%% SSP_PAC_MAIN is a helper which calls:
%   1) SSP_PAC: Use a constrained linear regression to estimate coupling
%   2) SSP_PAC_PROCESS: Gather Results and use a second state space to filter results if necesserary  
%   
%  INPUTS    -modelOsc  : structure obtained from SSP_DECOMP
%            -Kappa     : ~1x2 number of resampled time series for CI estimation 
%               Kappa(1): Resample Decomposition, Kappa(2) ressample regression 
%            -slowID    : corresponding to slower oscillation in x_t_n
%            -fastID    : corresponding to faster oscillation in x_t_n
%            -Nlag      : max lag used to resample oscillations time series
%            -Nconst (optionnal) : 1: constraint on K_mod, >1 linear, constraint, 0: none
%            -doArFiltering: , binary variable, use or not a 2nd state space to filter regressed_pac
%
%  OUTPUTS   -processed_ssp_pac: structure containing estimated parameters K_mod, Phi_mod,R2, filtered parameters and credible intervals.
%                For 1 window  : A_fast_t = A_0 (1 + K_mod cos(Phi_mod + Phi_slow_t))
%            -regressed_pac    : ~ 3 x Kappa+1 x Nwind the (resampled)
%                linear regression coeeficients beta (:,1,:) correspond to the
%                actual estimate.
%            -SS_tot           : ~3 x (Kappa+1) x Nwin sum of squares : [ESS  TSS  RSS]'
%            -X_t_tot          : 2 x Nosc x nT x (Kappa+1) Resampled oscillations time series (might be large)
%%         

if ~exist('Kappa_tot','var')
   Kappa_tot = [0,0];
end

if ~exist('slowID','var')
   slowID = 1;
end

if ~exist('fastID','var')
   fastID = 2;
end

if ~exist('Nlag','var')
   Nlag = Inf;
end

if ~exist('NlinConst','var')
   Nconst = 1;
end

if ~exist('doArFiltering','var')
   doArFiltering = 1;
end

if nargout>2
    [regressed_pac, SS_tot, X_t_tot] = ssp_pac (modelOsc, slowID, fastID, Kappa_tot, Nlag,Nconst);
else
    [regressed_pac, SS_tot] = ssp_pac (modelOsc, slowID, fastID, Kappa_tot, Nlag,Nconst);
end

tVect = 0.5*(modelOsc.startPointId+modelOsc.endPointId)/modelOsc.Fs;
processed_ssp_pac = ssp_pac_process(regressed_pac, SS_tot,doArFiltering,tVect,Nconst);


end







