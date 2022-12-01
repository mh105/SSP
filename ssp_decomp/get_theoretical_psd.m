function [H_tot, H_i]=get_theoretical_psd(f_y,Fs,freq_tot,ampl_tot_k,nois_tot_k)
%% Returns theoretical/parametric power spectral density (PSD) for a signal generated with Matsuda's and al. (2017)
%       f_y : frequencies at which PSD is calculated
%       Fs  : sampling frequency
%       freq_tot   : representant frequencies of the of the simulated signal
%       ampl_tot_k : associated amplitude
%       nois_tot_k : associated process noise covariance
%
%       H_i is the PSD of a given oscillation at frequency freq_tot(1,ii)
%       H_tot is the total PSD

Nfreq=length(freq_tot);

% Solve small issue when f_i = Fs /4 -> w = pi/2
for ww_i=1:length(freq_tot)
    if freq_tot(1,ww_i)==Fs/4
        freq_tot(1,ww_i)=freq_tot(1,ww_i)+freq_tot(1,ww_i)*0.002;
    end
end



z=exp(2*1i*pi*f_y/Fs);
w_tot= freq_tot *2*pi / Fs;


%nois2_tmp=diag(nois_tot_k);%nois2    =nois2_tmp(1:2:end)';
nois2    =nois_tot_k;


A_i = (1- 2.*ampl_tot_k.^2 .* cos(w_tot).^2+ampl_tot_k.^4.*cos(2*w_tot) )./ (ampl_tot_k.*(ampl_tot_k.^2-1).*cos(w_tot));
B_i = 0.5* (A_i- 2*ampl_tot_k.* cos(w_tot) +sqrt((A_i- 2*ampl_tot_k.*cos(w_tot)).^2-4 ));
V_i= -(nois2.*ampl_tot_k.*cos(w_tot))./B_i;


H_i= zeros(Nfreq, length(z));
for ii=1:Nfreq
    H_i(ii, :)=  (V_i(1,ii)/(Fs)).* abs(1+ B_i(1,ii).*z).^2./abs(1-2*ampl_tot_k(1,ii)*cos(w_tot(1,ii)).*z+ ampl_tot_k(1,ii).^2.*z.^2).^2;
end


H_tot= sum(H_i,1);

end