addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))

%% Generate a slow oscillation with harmonics and a faster one during n windows of length N.
Fs=250; % Hz
N=3*Fs; % Window length
n=1;    % Window Number
tt= (1:n*N);

StaPointId = ((1:n)-1)*N +1;
EndPointId = StaPointId+N-1;

f_so = 1;  w_so = f_so *2*pi /Fs; 
f_al = 10; w_al = f_al *2*pi /Fs;

x_so = 10 * cos(w_so*tt)+5 * cos(2*w_so*tt);
x_al = 2 * cos(w_al*tt);

sigma_noise= 1;
y_t = x_so +  x_al + normrnd(0,sigma_noise,[1,n*N]);
y_t=y_t';

figure
subplot(2,1,1); hold on
plot(tt/Fs, x_so )
plot(tt/Fs, x_al)

subplot(2,1,2); hold on
plot(tt/Fs, y_t, 'k')

%% Initialize and fit
convergenceTolerance=eps;
em_its=100;
doPlot=0;
namecur='test_fit';

init_params=struct();
init_params.f_init     =  [1.5, 6];
init_params.a_init     =  [0.98, 0.98];
init_params.sigma2_init=  [1,1];
init_params.R          =  1;
init_params.NNharm     =  [2, 1];
modelOsc=ssp_decomp(y_t,Fs,StaPointId,EndPointId,em_its,convergenceTolerance, init_params,namecur,doPlot);


%% Plot summary
plot_summary(modelOsc)


