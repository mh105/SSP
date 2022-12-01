%QUICKBANDPASS   A quick bandpass returns a the data filtered at a given frequency center or band using an equiripple bandpass filter
%
%   Usage:
%   [filtdata] = quickbandpass(data, Fs, bandpass_frequencies, bandwidth)
%
%   Input:
%   data:  EEG data
%   Fs:  sampling frequency
%   bandpass_frequencies:  frequency band or one frequency as center of
%   band
%   bandwidth: bandwidth around center frequency
%   
%   Output:
%   filtdata:  bandpass filtered EEG data
%   
%   Example:
% 
%   See also firls
%
%   Copyright 2011 Michael J. Prerau, Ph.D.
%   
%   Last modified 01/07/2011
%********************************************************************

function [filtdata tail] = quickbandpass(data, Fs, bandpass_frequencies,N ,bandwidth)
%A quick bandpass returns a the data filtered at a given frequency center or
%band using an equiripple bandpass filter
%
% filtdata = quickbandpass(data, Fs, freq_band)
% filtdata = quickbandpass(data, Fs, freq_center)
% filtdata = quickbandpass(data, Fs, freq_center, bandwidth)

%If only a frequency center is selected
if length(bandpass_frequencies)==1
    %Default bandwidth is 1Hz
    if nargin<5
        bandwidth=1;
    end
    
    %Specify frequency parameters
    Fpass1 = max(bandpass_frequencies-(bandwidth/2),.01);    % First Passband Frequency
    Fstop1 = max(Fpass1-.01,.01);                             % First Stopband Frequency
    
    Fpass2 = bandpass_frequencies+(bandwidth/2);    % Second Passband Frequency
    Fstop2 = Fpass2+.01;                             % Second Stopband Frequency
%If a band is selected
else
    Fpass1 = max(bandpass_frequencies(1),.01);   % First Passband Frequency
    Fstop1 = max(Fpass1-.01,.01);                 % First Stopband Frequency
    
    Fpass2 = bandpass_frequencies(2);   % Second Passband Frequency
    Fstop2 = Fpass2+.01;                 % Second Stopband Frequency
end



% Dstop1 = 0.001;           % First Stopband Attenuation
% Dpass  = 0.057501127785;  % Passband Ripple
% Dstop2 = 0.001;           % Second Stopband Attenuation
% dens   = 20;              % Density Factor

% % Calculate the order from the parameters using FIRPMORD.
% [N, Fo, Ao, W] = firpmord(max([Fstop1 Fpass1 Fpass2 Fstop2],1e-4)/(Fs/2), [0 1 ...
%     0], [Dstop1 Dpass Dstop2]);
% 
% % Calculate the coefficients using the FIRPM function.
% b  = firpm(N, Fo, Ao, W, {dens});
% Hd = dfilt.dffir(b);

%N      = 2000;   % Order

Wstop1 = 1;      % First Stopband Weight
Wpass  = 1;      % Passband Weight
Wstop2 = 1;      % Second Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2]);
Hd = dfilt.dffir(b);
delay=ceil(N/2);

%Filter the data and return result
filtfull=filter(Hd,[data(1)*ones(1,N) data]);
filtdata=[filtfull((delay+N+1):end) zeros(1,ceil(N/2))];
tail=ceil(N/2);