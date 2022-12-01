addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
%% Generate 3 oscillations during n windows of length N.
Fs=250; % Hz
N=3*Fs; % Window length
n=4;    % Window Number
tt= (1:n*N);

StaPointId = ((1:n)-1)*N +1;
EndPointId = StaPointId+N-1;

f_so = 1;  w_so = f_so *2*pi /Fs; 
f_al = 10; w_al = f_al *2*pi /Fs;
f_ga = 30; w_ga = f_ga *2*pi /Fs;

x_so = 10 * cos(w_so*tt);
x_al = 5 * cos(w_al*tt);
x_ga = 3 * cos(w_ga*tt);

sigma_noise= 1;
y_t = x_so +  x_al + x_ga + normrnd(0,sigma_noise,[1,n*N]);
y_t=y_t';

figure
subplot(2,1,1); hold on
plot(tt/Fs, x_ga)
plot(tt/Fs, x_al)
plot(tt/Fs, x_so)

subplot(2,1,2); hold on
plot(tt/Fs, y_t, 'k')


%% Initialize ssp and fit
convergenceTolerance=eps;
em_its=100;
doPlot=0;
namecur='test_fit';

init_params=struct();
init_params.f_init     =  [2, 6, 15];
init_params.a_init     =  [0.98, 0.98, 0.98];
init_params.sigma2_init=  [1,1,1];
init_params.R          =  1;
init_params.NNharm     =  [1, 1, 1];

modelOsc=ssp_decomp(y_t,Fs,StaPointId,EndPointId,em_its,convergenceTolerance, init_params,namecur,doPlot);

%% Plot summary
crange = [-20 20];
plot_summary(modelOsc, crange)




