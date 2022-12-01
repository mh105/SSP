%% Load a structure modelOsc beforehand (load eeg, modelOSC AND PAM analysis)
addpath(genpath('../dss_decomp'))
addpath(genpath('../spectral_analysis'))
addpath(genpath('../osc_pam_bootstrap'))

Fs=modelOsc.Fs;
startPointId=modelOsc.startPointId;
%startPointId=modelOsc.startPointMin;
endPointId=modelOsc.endPointId;
T=length(startPointId);
Nharm=modelOsc.init_params.Nharm;
%Nharm=0;
Nosc=floor(size(modelOsc.res{1,1}.x_t_n,1)/2);

f_tot=zeros(Nosc,T);
a_tot=zeros(Nosc,T);
q_tot=zeros(Nosc,T);
R_tot=zeros(1,T);
f_y=0.01:0.01:120;
s_tot=zeros(length(f_y),T);
s_tot_clean=zeros(length(f_y),T);
taTot=1:length(eeg);
params=struct();
params.Fs=Fs;
params.tapers=[4 5];


%% Mt spectrum
obj.Fs = Fs;
obj.tapers = [3 5];
obj.trialave = 1;
obj.fpass=[0 120];
[S_mt, stimes_mt, sfreqs_mt]= mtspecgramc(eeg,[5 5],obj);
ylmax=60;
ylmin=0;
clmin=-15;
clmax=20;




%% Old Way PAM

Npb=18;               % Phase Bin Number
pbins=linspace(-pi,pi,Npb+1); % Phase bins
lengthWind_PAM=120;
lengthWind=lengthWind_PAM*Fs;
N=floor(length(eeg)/lengthWind);
cutMeas=eeg(1: N*lengthWind);
resMeas=reshape(cutMeas, [lengthWind, N]);
PAMusual=zeros(Npb, N);
sofreq=[0.1 1]; alfreq=[8 12]; orderfilt=100;

for tt=1:N
    disp([num2str(tt),'/', num2str(N) ])
    [x_so_tmp, tail1] = quickbandpass(resMeas(:,tt)', Fs, sofreq,orderfilt);
    [x_al_tmp, tail2] = quickbandpass(resMeas(:,tt)', Fs, alfreq,orderfilt);
    Phi_filt_so_tmp=angle(hilbert(x_so_tmp));
    A_filt_al_tmp  =abs(hilbert(x_al_tmp));
    pam_filt_tmp=phaseamp(A_filt_al_tmp ,Phi_filt_so_tmp,pbins);
    PAMusual(:,tt)=pam_filt_tmp;
end

%% Parametric PAM
DecimateFactor=20;
K_mod=K_mod_t_k(1,:);
Phi_mod=Phi_mod_t_k(1,:);
starindDec=floor(startPointId/DecimateFactor)+1;
endindDec=floor(endPointId/DecimateFactor);

%% Plots


figure(1); hold on
set(gcf, 'position', [675   644   1.9*570   2.2*450]);
subplot(3,2,1); hold on
plot(taTot/(60*Fs),eeg, 'color', 'k')
title('EEG trace')
ylabel('\mu V')
xlabel('(min)');

subplot(3,2,3); hold on
imagesc(stimes_mt/60,  sfreqs_mt, pow2db(S_mt)');
axis xy;
ylabel('Frequency (Hz)');
xlabel('(min)');
title('Multitaper PSD');
caxis([clmin clmax])
colormap(jet)
ylim([0 30])


subplot(6,2,2)
FsMAX=DRUG.Fs;
plot(DRUG.spump.T/(60*FsMAX),DRUG.spump.prpfol,'color','k','linewidth',2)
title('Propofol Concentration')
ylabel('\mu g / mL')
xlabel('(min)');

subplot(6,2,4); hold on;
plot(BHVR.T'/(FsMAX*60),BHVR.prob_verbal.p500 ,'color','k','linewidth',2)
[fillhandle,msg]=jbfill(BHVR.T'/(FsMAX*60),BHVR.prob_verbal.p975',BHVR.prob_verbal.p025','r');
title('Verbal Response Probability')
ylabel('\mu g / mL')
xlabel('(min)');

pPAM=subplot(3,2,4)
imagesc((1:N)*lengthWind/(Fs*60),pbins,PAMusual);
col = redbluemap; colormap(pPAM,col);
set(gca, 'clim', [0.2 1.8]); %0.2-1.8
title('Usual PAM (120s windows)')
xlabel('(min)');
ylabel('(rad)')




%
for tt=1:T
    Phi_tmp=modelOsc.res{1,tt}.model_prams.Phi;
    Q_tmp=modelOsc.res{1,tt}.model_prams.Q;
    R_tmp=modelOsc.res{1,tt}.model_prams.R;
    y_t =modelOsc.res{1,tt}.y;
    x_so=modelOsc.res{1,tt}.x_t_n(1,2:end);
    x_al=modelOsc.res{1,tt}.x_t_n(3,2:end);
    
    [S_al,fal]=mtspectrumc(x_al,params);
    [S_so,fso]=mtspectrumc(x_so,params);
    [S_y,fy]  =mtspectrumc(y_t,params);
    
    TimeCur=(startPointId(tt)+endPointId(tt))*0.5/(Fs*60);
    
    for nosc=1:Nosc
       Phi_n= Phi_tmp((nosc-1)*2+1:2*nosc,(nosc-1)*2+1:2*nosc);
       [a_tmp,w_tmp]=get_rotPam(Phi_n);
       a_tot(nosc,tt)=a_tmp;
       f_tot(nosc,tt)=abs(w_tmp*Fs/(2*pi));
       q_tot(nosc,tt)=Q_tmp((nosc-1)*2+1,(nosc-1)*2+1);
    end
    R_tot(1,tt)=R_tmp;
    [H_tot, H_i]=getTheoreticalPSD(f_y,Fs,f_tot(:,tt)',a_tot(:,tt)',100+0*q_tot(:,tt));

  
    
    Kcur= mean(K_mod(1,starindDec(tt):endindDec(tt)),2);
    Phicur= mean(Phi_mod(1,starindDec(tt):endindDec(tt)),2);
    Npb2=600;theta=linspace(-pi, pi,Npb2);
    PAMequiv =  (1+2*sin(pi/Npb2)*Kcur.*cos(Phicur-theta')*Npb2/(2*pi));
    
    figure(1); hold on; 
    
    subplot(3,2,1); hold on;
    ylimTep=ylim;
    if exist('pLine'); delete(pLine);end;
    pLine=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',2, 'color', 'r');
    
    subplot(6,2,2); hold on;
    ylimTep=ylim;
    if exist('pLine2'); delete(pLine2);end;
    pLine2=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',2, 'color', 'r');
    
     subplot(3,2,3); hold on;
    ylimTep=ylim;
    if exist('pLine3'); delete(pLine3);end;
    pLine3=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',2, 'color', 'r');
    
    subplot(6,2,4); hold on;
    ylimTep=ylim;
    if exist('pLine4'); delete(pLine4);end;
    pLine4=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',2, 'color', 'r');
    
    subplot(3,2,4); hold on;
    ylimTep=ylim;
    if exist('pLine5'); delete(pLine5);delete(pLine55);end;
    pLine5=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',4, 'color', 'k');
    pLine55=plot([TimeCur TimeCur], [ylimTep(1) ylimTep(2)],'linewidth',2, 'color', 'r');
    
    
    
    p1=subplot(3,2,5);hold on;cla;
    title(['Parametric PSD. Time : ', num2str(TimeCur), ' min'])
    plot(f_y, 10*log10(H_tot/Fs),'linewidth',2, 'color', 'k')
    plot(f_y, 10*log10(H_i(1,:)/Fs),'linewidth',2, 'color', 'b');
    plot(f_y, 10*log10(H_i(2,:)/Fs),'linewidth',2, 'color', 'r');
    %plot(fy, 10*log10(S_so),':','linewidth',2,'color', 'b');
    %plot(fy, 10*log10(S_al),':','linewidth',2, 'color', 'r');
    ylim([-30 20])
    xlim([ylmin ylmax])
    
   
    subplot(3,2,6);hold on;cla;
    title(['Parametric PAM. Time : ', num2str(TimeCur), ' min'])
    plot(theta,PAMequiv,'linewidth',3, 'color', 'k');
    
    xlim([-pi pi]);
    ylim([0 2]);
    set(gca,'xtick',[0 1 2],'xtick',[-pi,0, pi],'XTickLabel',{'-\pi','0','\pi'}) 
    

    F(tt) = getframe(gcf) ;
    drawnow()
    
end

%%
  % create the video writer with 1 fps
  writerObj = VideoWriter([modelOsc.name,'.avi']);
  writerObj.FrameRate = 15;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for tt=1:T
    % convert the image to a frame
    frame = F(tt) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);