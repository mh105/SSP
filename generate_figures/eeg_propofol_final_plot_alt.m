addpath(genpath('../dss_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../osc_pam_bootstrap'))
addpath(genpath('../usual_pam'))

%% Load Dataset to be plotted
% dpath='/media/nsrlhs/f0e92df7-6f0d-4f06-bbeb-b6e03ea748fa/data/eeg_propofol_data/';
% load([dpath,'pam_eeg04_2018_6_12_18h53_6s_50ite.mat'])
% 0
% load([dpath,'pam_eeg04.mat'])
% 1
% load([dpath,'pam_eeg04_bootstrap_500.mat'])
% 2
% load([dpath,'pam_eeg04_bootstrap_500_final.mat'])
% 3

%% Retrieve params from modelOsc
startPointId=modelOsc.startPointId;
endPointId=modelOsc.endPointId;
timeWindNum=length(startPointId);
linRegOsc=struct();
windowLength=floor(mean(endPointId-startPointId+1));
Fs=modelOsc.Fs;
Kappa=size(Phi_al_final,2)-1;
T=length(startPointId);
Nharm=modelOsc.init_params.Nharm;
Nosc=floor(size(modelOsc.res{1,1}.x_t_n,1)/2);


%% Parametric Spectrum

f_tot=zeros(Nosc,T);
a_tot=zeros(Nosc,T);
q_tot=zeros(Nosc,T);
R_tot=zeros(1,T);
f_y=0.1:0.1:120;8
s_tot      =zeros(length(f_y),T);
s_tot_clean=zeros(length(f_y),T);


x_t_n_tot=zeros(2*Nosc, endPointId(1,end)-startPointId(1,1) );
y_tot=zeros(1, endPointId(1,end)-startPointId(1,1) );

for tt=1:T
    x_t_n_tot(:,startPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.x_t_n(:,2:end);    
    y_tot    (1,startPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.y;    
end

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

ta=(1:endPointId(end)-startPointId(1))/Fs;


%% Mt Spectrum


obj.Fs = Fs;
obj.tapers = [3 5];
obj.trialave = 1;
obj.fpass=[0 120];
[S_mt, stimes_mt, sfreqs_mt]= mtspecgramc(y_tot,[5 5],obj);


%% Old Way PAM

Npb=18;               % Phase Bin Number
pbins=linspace(-pi,pi,Npb+1); % Phase bins
lengthWind_PAM=120;
lengthWind=lengthWind_PAM*Fs;
N=floor(length(y_tot)/lengthWind);
cutMeas=y_tot(1: N*lengthWind);
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





%% Figure 2 : comparison 
DecimateFact=20; % Might be changed
ylmax=20;
ylmin=0;
clmin=-15;
clmax=20;
plotBOOT=0;


figure
p1=subplot(2,2,1); hold on
imagesc(stimes_mt/60,  sfreqs_mt, pow2db(S_mt)');
axis xy;
ylabel('Frequency (Hz)');
xlabel('(min)');
title('Multitaper PSD');
caxis([clmin clmax])
ylim([ylmin ylmax])

p2=subplot(2,2,2); hold on
imagesc((1:length(s_tot_clean))*windowLength/(60*Fs),  f_y, pow2db(s_tot_clean));
axis xy;
ylabel('Frequency (Hz)');
xlabel('(min)');
title('Parametric PSD');
caxis([clmin clmax])
ylim([ylmin ylmax])
linkaxes([p1,p2])
xlim([1*windowLength/(60*Fs), length(s_tot_clean)*windowLength/(60*Fs)])


colormap(p1, jet)
colormap(p2, jet)



pos1=get(subplot(2,2,1), 'Position');
pos2=get(subplot(2,2,2), 'Position');

highBar=pos1(4);
widtBar=highBar/20;
middleFig=(pos1(1)+pos2(1)+pos2(3))/2-widtBar/2;

cax = colorbar('Position', [middleFig pos1(2)  widtBar/1.2 highBar]); %caxis([0.5 1.5]);
ylabel(cax,'Power (dB)');




p3=subplot(4,2,6); hold on
if plotBOOT
plot(DecimateFact*(1:length(K_mod_t_k(1,:)))/(60*Fs),K_mod_t_k(:,:), 'color', [0.5 0.5 0.5])
end
plot(DecimateFact*(1:length(K_mod_t_k(1,:)))/(60*Fs),K_mod_t_k(1,:), 'color', [0    0.4470    0.7410])
title(' K^{mod}_t')
axis tight

p4=subplot(4,2,8); hold on
if plotBOOT
    for kk=1:Kappa+1
        scatter(DecimateFact*(1:length(K_mod_t_k(1,:)))/(60*Fs), Phi_mod_t_k(kk,:), 0.001, [0.5 0.5 0.5])
    end
end
scatter(DecimateFact*(1:length(K_mod_t_k(1,:)))/(60*Fs), Phi_mod_t_k(1,:), 1, [0    0.4470    0.7410])
axis tight

ylim([-pi pi])
title('\phi^{mod}_t')
xlabel('(min)');
ylabel('(rad)')

p5=subplot(4,2,5);
imagesc((1:N)*lengthWind/(Fs*60),pbins,PAMusual);
col = redbluemap; colormap(p5,col);
set(gca, 'clim', [0.2 1.8]); %0.2-1.8
title('Usual PAM (120s windows)')
ylabel('(rad)')

Npb2=18;
theta=linspace(-pi, pi,Npb2);
PAMequiv =  (1+2*sin(pi/Npb)*mean(K_mod_t_k(:,:),1).*cos(mean(Phi_mod_t_k(1,:),1)+theta')*Npb/(2*pi));

p6=subplot(4,2,7);
imagesc((1:N)*lengthWind/(Fs*60),theta,PAMequiv);
col = redbluemap; cax=colormap(p6,col);
set(gca, 'clim', [0 2]); %0.2-1.8
set(gca, 'clim', [0.2 1.8]); %0.2-1.8
title('Parametric PAM (6s windows)')
xlabel('(min)');
ylabel('(rad)')
%axis xy

pos1=get(subplot(4,2,7), 'Position');
pos2=get(subplot(4,2,8), 'Position');


middleFig=(pos1(1)+pos2(1)+pos2(3))/2-widtBar/2;
cax = colorbar('Position', [middleFig pos1(2)  widtBar/1.2 highBar]); %caxis([0.5 1.5]);
set(gca, 'clim', [0 2]); %0.2-1.8
ylabel(cax,'Alpha power Distribution');
uu=colormap(cax,col);
linkaxes([p1,p2,p3,p4,p5,p6], 'x')

p6=subplot(4,2,7);
xlim([20 120])

h=gcf;
set(h,'Position',[50 50 1200 800]);
set(h,'PaperOrientation','landscape');
%print(gcf, '-dpdf', 'figure2subject03_300ite2.pdf')



%% Thesis Figure 1



scaleTime=1;
scalemuVolt=20;

WOI=468;
WOI=735;
WOI=905;
%WOI=800;
WOI=904;
tspan=startPointId(WOI):endPointId(WOI);
taDecimated=ta(1:DecimateFact:end);
ta=tspan/(Fs*60);
y=eeg(1,tspan);

figure(5); cla;
subplot(4,1,1);cla; hold on;
y=eeg(1,tspan);
plot((1:length(eeg))/(Fs*60),eeg, 'color', 'k')
xlabel('(min)')
ylabel('\muV')

pos1=get(subplot(4,2,3), 'Position');
pos2=get(subplot(4,2,4), 'Position');
width1=pos1(3);heigh1=pos1(4);
ypos =pos1(2);xpos =(pos2(3)+pos2(1)+pos1(1)-width1)/2;
posiRaw=[xpos, ypos, width1, heigh1];
%posiRaw=get(subplot(4,2,3), 'Position'); subplot(4,2,4);set(gca, 'Visible', 'off')
f1=subplot('Position', posiRaw); hold on;
plot(ta,y, 'color', 'k')
set(gca,'xtick',[],'ytick',[])

x_so=X_t_so_final(tspan, :);
x_al=X_t_al_final(tspan, :);
[~, ind95]=sort( sum(abs(x_so-x_so(:,1))+abs(x_al-x_al(:,1)),1) );
max95=floor(size(x_so,2)*0.90);

x_so=X_t_so_final(tspan, ind95);
x_al=X_t_al_final(tspan, ind95);
phi_so=Phi_so_final(tspan,ind95);
phi_al=Phi_al_final(tspan,ind95);
a_al=x_al./cos(phi_al);



f3=subplot(4,2,6);cla; hold on;
%plot(ta, x_al(:,1:max95), 'Color', [0.5 0.5 0.5])
plot(ta, x_al(:,1), 'Color', 'r')
%title('x^{al}_t')

f2=subplot(4,2,5);cla; hold on;
%plot(ta, x_so(:,:), 'Color', [0.5 0.5 0.5])
plot(ta, x_so(:,1), 'Color', 'b')
title('x^{so}_t')

linkaxes([f1 f2 f3])
axis tight


subplot('Position', posiRaw);
axis tight
box on
ylimTmp=ylim;
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])


subplot(4,2,5)
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10 ,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')
subplot(4,2,6)
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')


subplot(4,2,8);cla; hold on;
plot(ta, a_al(:,1:max95),':', 'Color', [0.5 0.5 0.5], 'linewidth', 0.05)
plot(ta, a_al(:,1), 'Color', 'r')
title('A_t^{\alpha}')
ylimTmp=ylim;
plot([ta(1); ta(1)], [0; 0+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [0;0], '-k', 'LineWidth', 2)
text(ta(1),0-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')

subplot(4,2,7);cla; hold on;
for kk=1:max95+1
scatter(ta, phi_so(:,kk),1, [0.5 0.5 0.5])
end
scatter(ta, phi_so(:,kk),1,'b')
axis tight
ylim([-pi pi])
set(gca,'xtick',[],'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'}, 'FontSize', 16)
box on
ylimTmp=ylim;
plot([ta(1); ta(scaleTime*Fs)], [-pi;-pi], '-k', 'LineWidth', 2)
text(ta(1),-pi-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s'])


set(gcf,'Position',[961           1         640        730]);
pathDrop='/home/nsrlhs/Dropbox (MIT)/producedFig/';
print(gcf, '-dpdf', [pathDrop,'decomposotionCenter300ite.pdf'])

%%

