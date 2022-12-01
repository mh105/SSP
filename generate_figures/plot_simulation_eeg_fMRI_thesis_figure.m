addpath(genpath('../dss_decomp'))
addpath(genpath('../spectralanalysis'))

dpath='../../data/eeg_fmri_simulation/';


load([dpath,'master_thesis_generatedComa.mat']);
load([dpath,'master_theis_simulation_2018_6_15_50ite.mat']);

Fs=modelOsc.Fs;
starPointId=modelOsc.startPointId;
endPointId=modelOsc.endPointId;
T=length(starPointId);
Nosc=floor(size(modelOsc.res{1,1}.x_t_n,1)/2);

x_t_n_tot=zeros(2*Nosc, endPointId(1,end)-starPointId(1,1) );

for tt=1:T
    x_t_n_tot(:,starPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.x_t_n(:,2:end);    
end

eeg_noisy=generatedComa.eeg_bcg;
bcg_noise=generatedComa.bcg_background;
eeg_theor=generatedComa.eeg;
eeg_clean=x_t_n_tot(1,:)+x_t_n_tot(3,:);%+x_t_n_tot(5,:);
bcg_fit  =sum(x_t_n_tot(7:2:end,:),1);
%bcg_fit  =eeg_noisy-eeg_clean;

%% Mt Spectrum)-startPointId(1)));
ta=(1:endPointId(end)-starPointId(1))/Fs;

obj.Fs = Fs;
obj.tapers = [3 5];
obj.trialave = 1;
obj.fpass=[0 120];

[S_noisy, stimes_noisy, sfreqs_noisy]= mtspecgramc(eeg_noisy,[6 5],obj);
[S_theor, stimes_theor, sfreqs_theor]= mtspecgramc(eeg_theor,[6 5],obj);
[S_bcg  , stimes_bcg  , sfreqs_bcg  ]= mtspecgramc(bcg_noise,[6 5],obj);
[S_clean, stimes_clean, sfreqs_clean]= mtspecgramc(eeg_clean,[6 5],obj);


%% Parametric Spectrum

Fs=modelOsc.Fs;
startPointId=modelOsc.startPointId;
endPointId=modelOsc.endPointId;
T=length(startPointId);
Nharm=modelOsc.init_params.Nharm;
Nosc=floor(size(modelOsc.res{1,1}.x_t_n,1)/2);

f_tot=zeros(Nosc,T);
a_tot=zeros(Nosc,T);
q_tot=zeros(Nosc,T);
R_tot=zeros(1,T);
f_y=0.1:0.1:120;8
s_tot      =zeros(length(f_y),T);
s_tot_clean=zeros(length(f_y),T);

for tt=1:T
    Phi_tmp=modelOsc.res{1,tt}.model_prams.Phi;
    Q_tmp=modelOsc.res{1,tt}.model_prams.Q;
    R_tmp=modelOsc.res{1,tt}.model_prams.R;
    
    for nosc=1:Nosc
       Phi_n= Phi_tmp((nosc-1)*2+1:2*nosc,(nosc-1)*2+1:2*nosc);
       [a_tmp,w_tmp]=get_rotPam(Phi_n);
       a_tot(nosc,tt)=a_tmp;
       f_tot(nosc,tt)=w_tmp*Fs/(2*pi);
       q_tot(nosc,tt)=Q_tmp((nosc-1)*2+1,(nosc-1)*2+1);
    end
    
    R_tot(1,tt)=R_tmp;
    [H_tot, H_i]=getTheoreticalPSD(f_y,Fs,f_tot(:,tt)',a_tot(:,tt)',q_tot(:,tt));
    s_tot(:,tt)=H_tot';
    s_tot_clean(:,tt)=H_i(1,:)+H_i(2,:);
    
end



%%
ylmax=20;
ylmin=0;
clmin=-20;
clmax=40;


figure
p1=subplot(3,3,1);
imagesc(stimes_noisy/60, sfreqs_noisy, pow2db(S_noisy)');
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (min)');
ta=(1:endPointId(end)-startPointId(end));
title('mt PSD raw');
%c=colorbar;
caxis([clmin clmax])
%ylabel(c,'Power (dB)');
ylim([ylmin ylmax])

p2=subplot(3,3,2);
imagesc(stimes_theor/60, sfreqs_theor, pow2db(S_theor)');
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (min)');
title('mt PSD theoretical');
%c=colorbar;
caxis([clmin clmax])
%ylabel(c,'Power (dB)');
ylim([ylmin ylmax])


p3=subplot(3,3,5);
imagesc(stimes_clean/60, sfreqs_clean, pow2db(S_clean)');
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (min)');
title('mt PSD cleaned EEG');
%c=colorbar;
caxis([clmin clmax])
%ylabel(c,'Power (dB)');
axis tight
ylim([ylmin ylmax])

p4=subplot(3,3,4);
imagesc(startPointId/(60*Fs),f_y, pow2db(s_tot))
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (min)');
title('Parametric PSD EEG');
%c=colorbar;
caxis([clmin clmax])
%ylabel(c,'Power (dB)');
ylim([ylmin ylmax])

colormap(jet)

subplot(2,3,3);
c=colorbar('WestOutside');
ylabel(c,'Power [dB]');
%c.Label.String= 'Power '
caxis([clmin clmax])



linkaxes([p1,p2,p3,p4])
colormap(jet)
axis off

ta=(startPointId(1): startPointId(end))/(Fs);
subplot(3,3,7);hold on
TTp=13;
taP=(startPointId(TTp)+1*Fs: startPointId(TTp)+Fs*4);
plot(ta(taP), detrend(eeg_noisy(taP)), 'color', 'k')
plot(ta(taP), detrend(eeg_theor(taP)), 'color', [0  0.4470    0.7410])
plot(ta(taP), detrend(eeg_clean(taP)), 'color', [0.8500    0.3250    0.0980],'linewidth',1)
axis tight
xlabel('Time (sec)')
ylabel(['au'])

subplot(3,3,8);hold on
TTp=25;
taP=(startPointId(TTp)+1*Fs: startPointId(TTp)+Fs*4);
plot(ta(taP), detrend(eeg_noisy(taP)), 'color', 'k')
plot(ta(taP), detrend(eeg_theor(taP)), 'color', [0  0.4470    0.7410])
plot(ta(taP), detrend(eeg_clean(taP)), 'color', [0.8500    0.3250    0.0980],'linewidth',1)
axis tight
xlabel('Time (sec)')
ylabel(['au'])
legend('Raw','True', 'Cleaned')

   
%set(gcf,'PaperOrientation','landscape');print(gcf, '-dpdf', 'test1.pdf')

