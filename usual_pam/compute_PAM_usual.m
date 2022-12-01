function [pa_usual]=compute_PAM_usual(data,Fs,pbins)
phasefreq = [0.5 2]; ampfreq = [14 22]; order=500;
phasefreq = [0.1 1]; ampfreq = [8 12]; order=500;

Npb=length(pbins)-1;

[x_slow, tail1] =quickbandpass(data',Fs,phasefreq, order);   %Bandpass at the phase frequency
[x_alpha, tail2]=quickbandpass  (data',Fs,ampfreq, order);     %Bandpass at the amplitude frequency

x_al_hilb= hilbert(x_alpha(1:end-tail1));
env_Al   =abs(x_al_hilb);
x_so_hil=hilbert(x_slow(1:end-tail2));
phase_so= angle(x_so_hil);
pa_usual=phaseamp(env_Al,phase_so,pbins);
%pa_usual=(Npb/(2*pi))*pa_usual/sum(pa_usual);




end