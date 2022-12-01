function [f_tot, a_tot, sigma2_tot, R_tot, s_tot_param, s_param , Nosc_ind,NNharm, s_mt,s_mt_tot, fmt] = gather_params(modelOsc,freq)
%% GATHER_PARAMS grasps each oscillation parameters
%   INPUTS  : modelOsc, returned by SSP_DECOMP
%             freq,     frequency at which we estimate the parametric power spectral densities
%
%   OUTPUTS : f_tot       : peak frequency of each oscillations (e.o)
%             a_tot       : scaling parameter of e.o
%             sigma2_tot  : process noise of e.o
%             s_tot_param : total parametric PSD
%             s_param     : parametric psd of e.o
%             Nosc_ind    : numer of independatn oscillations
%             NNharm      : number of harmonic component for e.o
%             s_mt        : multitaper spectrogram of e.o
%             s_mt_tot    : total mt psd
%             fmt         : frequencies at which mt_psd is estimated

if nargin <2
    freq=0.01:0.01:120;
end

staPointId = modelOsc.startPointId;
endPointId = modelOsc.endPointId;
Nwind  = length(staPointId);
Fs     = modelOsc.Fs;

NNharm=modelOsc.init_params.NNharm(1,:);
Nosc     = sum(NNharm);    % Total number of oscillation
Nosc_ind = length(NNharm); % Number of indpt oscillations

specparams.Fs     = Fs;
specparams.tapers = [5 6]; %specparams.tapers = [1 1];
tpoints=mean(staPointId-endPointId);
nfft   =max(2^(nextpow2(tpoints)),tpoints)/2+1;

f_tot=zeros(Nosc,Nwind);
a_tot=zeros(Nosc,Nwind);
sigma2_tot=zeros(Nosc,Nwind);
R_tot=zeros(1,Nwind);

s_tot_param=zeros(length(freq),Nwind);
s_param    =zeros(length(freq),Nwind,Nosc);
s_mt       =zeros(nfft,Nwind);
s_mt_tot   =zeros(nfft,Nwind,Nosc);



for tt=1:Nwind
    Phi_tmp=modelOsc.res{1,tt}.model_prams.Phi;
    Q_tmp=modelOsc.res{1,tt}.model_prams.Q;
    R_tmp=modelOsc.res{1,tt}.model_prams.R;
    
    for nosc=1:Nosc
        Phi_n= Phi_tmp((nosc-1)*2+1:2*nosc,(nosc-1)*2+1:2*nosc);
        [a_tmp,w_tmp]=get_rotmat_pam(Phi_n);
        a_tot(nosc,tt)=a_tmp;
        f_tot(nosc,tt)=w_tmp*Fs/(2*pi);
        sigma2_tot(nosc,tt)=Q_tmp((nosc-1)*2+1,(nosc-1)*2+1);
        
        s_i_tmp=mtspectrumc(modelOsc.res{1,tt}.x_t_n((nosc-1)*2+1,2:end),specparams); 
        s_mt_tot(:,tt,nosc) =s_i_tmp;
    end
    
    R_tot(1,tt)=R_tmp;

    [H_tot, H_i]=get_theoretical_psd(freq,Fs,f_tot(:,tt)',a_tot(:,tt)',sigma2_tot(:,tt)');
    s_tot_param(:,tt)=H_tot';
    s_param(:,tt,:)=H_i';

    [s_mt_tmp,fmt]=mtspectrumc(modelOsc.res{1,tt}.y,specparams);
    s_mt(:,tt)=s_mt_tmp;
end



end