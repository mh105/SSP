addpath(genpath('../ssp_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../ssp_init'))
%% Generate 3 oscillations during n windows of length N and add an AR(1)
Fs=250; % Hz
N=3*Fs; % Window length
n=1;    % Window Number
ta= (1:n*N);

StaPointId = ((1:n)-1)*N +1;
EndPointId = StaPointId+N-1;

f_so = 1;  w_so = f_so *2*pi /Fs; 
f_al = 10; w_al = f_al *2*pi /Fs;
f_ga = 30; w_ga = f_ga *2*pi /Fs;

x_so = 10 * cos(w_so*ta);
x_al = 5 * cos(w_al*ta);
x_ga = 3 * cos(w_ga*ta);

ar_1 = 0.8;
x_ar1 = zeros(1,n*N);
eps_t = normrnd(0,1, 1,n*N);
for tt=2:n*N
    x_ar1(1,tt) = x_ar1(1,tt-1)*ar_1 + eps_t(1,tt);
end

sigma_noise= 1;
y_t = x_so +  x_al + x_ga + x_ar1 + normrnd(0,sigma_noise,[1,n*N]);
y_t=y_t';

figure
subplot(2,1,1); hold on
plot(ta/Fs, x_ga)
plot(ta/Fs, x_al)
plot(ta/Fs, x_so)
plot(ta/Fs, x_ar1, 'linewidth',2, 'color','m')

subplot(2,1,2); hold on
plot(ta/Fs, y_t, 'k')




%% Estimate Initialization Parameters
TW = 1;
Ktapers = 4;
Nosc_max = 6;
doRedress= 1;
doPlot = 1;

[f_tot, a_tot, sigma2_tot,R_hat, contrib_tot]=get_init_osc(y_t,StaPointId,EndPointId,Fs,TW,Ktapers,Nosc_max, doRedress,doPlot);
[~,osc_id]=sort(contrib_tot,'descend');


%% Oscillation Decomposition and IC to pull out the best number of oscillation
convergenceTolerance=eps;
em_its=500;
doPlot=0;
namecur='test_fit';

Nosc_max = 6;
ll_tot = zeros(1,Nosc_max);

selectionCrit ='BIC';

ICmin = Inf;
for Nosc=1:Nosc_max
    disp(['Nosc tested: ' , num2str(Nosc), '/', num2str(Nosc_max)])
    
    init_params=struct();
    init_params.f_init     =  f_tot(1,osc_id(1:Nosc));
    init_params.a_init     =  a_tot(1,osc_id(1:Nosc));
    init_params.sigma2_init=  sigma2_tot(1,osc_id(1:Nosc));
    init_params.R          =  R_hat;
    init_params.NNharm     =  ones(1,Nosc);
    
    modelOsc=ssp_decomp(y_t,Fs,StaPointId,EndPointId,em_its,convergenceTolerance, init_params,namecur,doPlot);
    ll_tot(1,Nosc) =modelOsc.res{1,1}.ll;
    
    AIC_cur = -2*ll_tot(1,Nosc)+2*(3*(Nosc)+1);
    BIC_cur = -2*ll_tot(1,Nosc)+(3*(Nosc)+1) * log(N);
    
    if strcmp(selectionCrit,'AIC')
        IC = AIC_cur;
    elseif strcmp(selectionCrit,'BIC')
        IC = BIC_cur;
    end
    
    if IC<ICmin
        ICmin = IC;
        modelOsc_best = modelOsc;
        Nosc_opt = Nosc;
    end
 
end

%% Plot AIC and BIC
AIC = -2*ll_tot+2*(3*(1:Nosc_max)+1);
BIC = -2*ll_tot+(3*(1:Nosc_max)+1) * log(N);

figure;
yyaxis left
plot(1:Nosc_max,AIC)
xlabel('Number of Oscillation')
ylabel('AIC')
yyaxis right
plot(1:Nosc_max,BIC)
xlabel('Number of Oscillation')
ylabel('AIC')

%% Plot summary of the best model
plot_summary(modelOsc_best)