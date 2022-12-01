addpath(genpath('../dss_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../dss_pac'))
addpath(genpath('../usual_pam'))

Fs=250;
Nwindow=5000;
lengthWindSec=6;
lengthTot=Nwindow*lengthWindSec*Fs;
t_tot=1:lengthTot;
starind=t_tot(1:lengthWindSec*Fs:end-lengthWindSec*Fs+1)';
stopind=t_tot(lengthWindSec*Fs:lengthWindSec*Fs:end)';

gferr

Nstep=6;
Tstep=linspace(1,lengthTot,Nstep+1)';
StepFunction=sum(heaviside(t_tot-Tstep(1:end-1)),1)/Nstep;


K_th  =0.9*ones(1,lengthTot);
Phi_th=(pi/3)*ones(1,lengthTot);

w_so_th= (2*pi/Fs)*linspace(0.5,0.5 ,lengthTot);
w_al_th= (2*pi/Fs)*floor(10*linspace(12,12 ,lengthTot))/10;

%RNS=linspace(0,40 ,lengthTot);
RNS=40*StepFunction;


phi_so_th=   w_so_th.*t_tot;
phi_so_th_n= w_so_th.*t_tot+w_so_th.*normrnd(0,RNS);%+noise

phi_al_th= w_al_th.*t_tot;
phi_al_th_n= w_al_th.*t_tot+w_so_th.*normrnd(0,RNS);%+noise

x_so_th=10*cos(phi_so_th_n);
x_al_th=5*cos(phi_al_th_n);
x_al_th_mod=x_al_th.*(1+K_th.*cos(phi_so_th_n+Phi_th));

y_tot_th= x_so_th+ x_al_th_mod +normrnd(0,RNS);%+noise obs




%figure; plot(t_tot/Fs, y_tot_th)

StepIndexes=find(diff(StepFunction));
StepSize=floor(mean(diff(StepIndexes)));
timeWind=4*Fs;
SampleStart=StepIndexes+timeWind;
SampleEnd=SampleStart+timeWind;
tickVal=unique(RNS)
figure
ylimTmp=[-120 120];

for jj=1:length(SampleStart)
    subplot(3,length(SampleStart),jj ); hold on
    ta=t_tot(SampleStart(jj):SampleEnd(jj))/Fs; ta=ta-ta(1);
    Ycurr=y_tot_th(SampleStart(jj):SampleEnd(jj));
    x_so_cur=x_so_th(SampleStart(jj):SampleEnd(jj));
    x_al_cur=x_al_th_mod(SampleStart(jj):SampleEnd(jj));
    plot(ta,  Ycurr, 'color', 'k')
    scalemuVolt=20; scaleTime=1;
    axis off;title(['\sigma^{2}_N = ',num2str(floor(tickVal(jj+1)))]); ylim([ylimTmp]);
    if jj==1
    plot([ta(1); ta(1)], [min(Ycurr); min(Ycurr)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(Ycurr); min(Ycurr)], '-k', 'LineWidth', 3)
    text(ta(1),min(Ycurr)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end
   
    
    
    subplot(3,length(SampleStart),jj+ length(SampleStart)*2); hold on
    plot(ta,  x_so_cur, 'color', 'b')
    plot(ta,  x_al_cur, 'color', 'R')
    ylimTmp2=ylimTmp/4;
    ylim([ylimTmp2]);
    if jj==1
    plot([ta(1); ta(1)], [min(x_so_cur); min(x_so_cur)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so_cur); min(x_so_cur)], '-k', 'LineWidth', 2)
    text(ta(1),min(x_so_cur)-(ylimTmp2(2)-ylimTmp2(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end
    axis off;
    

    
    subplot(3,length(SampleStart),jj+ length(SampleStart)); hold on
    plot(ta,  Ycurr-x_al_cur-x_so_cur, 'color', [0.5 0.5 0.5])
    ylim([ylimTmp]);

    axis off;
    if jj==1
    plot([ta(1); ta(1)], [min(Ycurr); min(Ycurr)+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(Ycurr); min(Ycurr)], '-k', 'LineWidth', 3)
    text(ta(1),min(Ycurr)-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), 'a.u'])
    end
    
    
    
end
set(gcf, 'Position', [ 0.6667    0.6667  853.3333  469.3333])
set(gcf,'PaperOrientation','landscape');
%print(gcf, '-dpdf', ['/Users/apple/Desktop/figure4supalt'])



%% Osc decomposition
em_its=300;
init_params=struct();
init_params.f_init     =  [2, 8];
init_params.a_init     =  [0.98, 0.98];
init_params.sigma2_init=  [1, 1];
init_params.R          =  10;
init_params.Nharm      =  0;
modelOsc=dss_decomp(y_tot_th',Fs,starind(:),stopind(:),em_its,eps, init_params,'',1);
save('./modelOscIncreasingNoise','modelOsc','-v7.3')

%% Gather Osc
% main_pam_from_oscDecomp

%% Filter Windows500
sofreq=[0.1 1.2]; alfreq=[8 14]; orderfilt=50;
lengthWindScan= [6  120];
orderfiltScan = [50 300];
Npb=18*2;
PAM_MI=cell(1,length(lengthWindScan));

pbins=linspace(-pi,pi,Npb+1);
pbins2=linspace(-pi,pi,Npb);

for ll=1:length(lengthWindScan)
    lengthWindCur=lengthWindScan(ll);
    PAM_MI{1,ll}.lengthWindScan=lengthWindScan(ll);
    

    lengthWind=lengthWindCur*Fs;
    Ncur=floor(length(y_tot_th)/lengthWind);
    cutMeas=y_tot_th(1: Ncur*lengthWind);
    resMeas=reshape(cutMeas, [lengthWind, Ncur]);
    PAMusual=zeros(Npb, Ncur);
    
    for tt=1:Ncur
        disp(['Window length :' ,num2str(ll),'/', num2str(length(lengthWindScan)),' || ' ,num2str(tt),'/', num2str(Ncur) ])
        orderfilt_cur=orderfiltScan(ll);
        
        [x_so_tmp, tail1] = quickbandpass(resMeas(:,tt)', Fs, sofreq,orderfilt_cur);
        [x_al_tmp, tail2] = quickbandpass(resMeas(:,tt)', Fs, alfreq,orderfilt_cur);
        
        Phi_filt_so_tmp=angle(hilbert(x_so_tmp(1,1:end-tail1)));
        A_filt_al_tmp  =abs(hilbert(x_al_tmp(1,1:end-tail2)));
        
        pam_filt_tmp=phaseamp(A_filt_al_tmp ,Phi_filt_so_tmp,pbins);
        PAMusual(:,tt)=pam_filt_tmp;
    end
    
    uni=ones(Npb,1);
    MIusual_cur=sum(abs(PAMusual-uni),1)/Npb;  %L1 norm %MIusual_cur=sum(PAMusual.*log(PAMusual),1);%KL divergence
    PAM_MI{1,ll}.N=Ncur;
    PAM_MI{1,ll}.MI= MIusual_cur;
    PAM_MI{1,ll}.PAMusual=PAMusual;
    
end



%% Equiv MI pAM
uni=ones(Npb,1);
K_mod=K_mod_t_k(1,:);
Phi_mod=Phi_mod_t_k(1,:);
dPhi=2*pi/Npb;
MtPhi=    (1 + 2*K_mod.*cos(pbins2'+Phi_mod)*sin(dPhi/2)/(2*pi/Npb));
MtPhi_th= (1 + 2*K_th .*cos(pbins2'+Phi_th) *sin(dPhi/2)/(2*pi/Npb));
MI_th  = sum(abs(MtPhi_th-uni),1)/Npb;

MtPhi_t_k=zeros(Npb, size(K_mod_t_k,2), size(K_mod_t_k,1));
MI_osc_t_k=zeros(size(K_mod_t_k,2), size(K_mod_t_k,1));
for kk=1:size(K_mod_t_k,1)
    
    MtPhi_t_k(:,:,kk)=(  1 + 2*K_mod_t_k(kk,:).*cos(pbins2'+K_mod_t_k(kk,:))*sin(dPhi/2)/(2*pi/Npb));
    MI_osc_t_k(:,kk)= sum(abs(MtPhi_t_k(:,:,kk)-uni),1)/Npb;
end



%% Plots
DecimateFactor=60;
colorFilt=[[0.8 0.8 0.8];[0.4 0.4 0.4];[0.8 0 0]]; lineW= [1 , 1.3, 1.3];


figure;

% Plot the Evolution of MI estimate
subplot(2,2,1); hold on;
for ll=1:length(lengthWindScan)
    N=PAM_MI{1,ll}.N;
    MI=PAM_MI{1,ll}.MI;
    PAMusual=PAM_MI{1,ll}.PAMusual;
    lengthWindCur=lengthWindScan(ll);
    plot((1:N)*lengthWindCur/60,MI, 'color', colorFilt(ll,:),'linewidth',lineW(ll));
end
plot(t_tot(1:DecimateFactor:end)/(60*Fs),MI_osc_t_k(:,1), 'color', 'b','linewidth',1.3);
plot(t_tot/(60*Fs),MI_th,'-.' , 'color', 'k','linewidth',1.1)
ylim([0 0.7])
box on

   
% Plot the Evolution of Phi estimate
subplot(2,2,2); hold on;
for ll=1:length(lengthWindScan)
    N=PAM_MI{1,ll}.N;
    MI=PAM_MI{1,ll}.MI;
    PAMusual=PAM_MI{1,ll}.PAMusual;
    lengthWindCur=lengthWindScan(ll);
    [~,locM]=max(PAMusual);
    phi_est_filt=pbins2(locM);
    plot((1:N)*lengthWindCur/60,-phi_est_filt, 'color', colorFilt(ll,:),'linewidth',lineW(ll));
end
plot(t_tot(1:DecimateFactor:end)/(Fs*60), Phi_mod_t_k(1,:), 'color', 'b', 'linewidth',1.2)
xlabel('Noise Factor')
ylabel('rad')
plot(t_tot/(60*Fs),mean(Phi_th)*ones(1,length(t_tot)),'-.' , 'color', 'k','linewidth',1.2)
set(gca,'YTick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'})
ylim([-pi pi])
legend({'Standard 6s','Standard 120s','Param 6s','True'}, 'Location', 'southwest')



StepIndexes=find(diff(StepFunction));
StepSize=floor(mean(diff(StepIndexes)));
ddelay=1.5*60*Fs;
DecimatedDelay=floor(ddelay/DecimateFactor);
StepStart=StepIndexes+ddelay;
StepEnd=StepStart+StepSize-1-ddelay;

DecimatedStart=floor(StepStart/DecimateFactor)+1;
DecimatedStop =floor(StepEnd/DecimateFactor);

% Plot black rectangles
subplot(2,2,1)
for indexS=1:length(StepIndexes)+1
   if indexS==1
   rectangle('Position', [eps 0 StepStart(1)/(Fs*60) 1],'FaceColor','k')
   elseif indexS==length(StepIndexes)+1
      rectangle('Position', [StepEnd(indexS-1)/(Fs*60) 0 ddelay/(60*Fs) 1],'FaceColor','k') 
   else
       rectangle('Position', [StepEnd(indexS-1)/(Fs*60) 0 2*ddelay/(60*Fs) 1],'FaceColor','k')
   end
end
tickPos=0.5*(StepEnd+StepStart)/(Fs*60);
tickval=floor(unique(RNS));
set(gca,'xtick',tickPos,'xticklabels', tickval(2:end) );
xlim([t_tot(1)/(60*Fs) t_tot(end)/(60*Fs)])

subplot(2,2,2)
for indexS=1:length(StepIndexes)+1
   if indexS==1
   rectangle('Position', [eps -pi StepStart(1)/(Fs*60) 2*pi],'FaceColor','k')
   elseif indexS==length(StepIndexes)+1
      rectangle('Position', [StepEnd(indexS-1)/(Fs*60) -pi ddelay/(60*Fs) 2*pi],'FaceColor','k') 
   else
       rectangle('Position', [StepEnd(indexS-1)/(Fs*60) -pi 2*ddelay/(60*Fs) 2*pi],'FaceColor','k')
   end
end
tickPos=0.5*(StepEnd+StepStart)/(Fs*60);
tickval=floor(unique(RNS));
set(gca,'xtick',tickPos,'xticklabels', tickval(2:end) );
xlim([t_tot(1)/(60*Fs) t_tot(end)/(60*Fs)])
box on



% BoxPlot MIi estimate
subplot(2,2,3); hold on
x=[];g=[];colorsBox=[];comp=1;
stdMI=[];

for indexS=1:length(StepIndexes)

    estMI=MI_osc_t_k(DecimatedStart(indexS):DecimatedStop(indexS),1);
    x=[x;estMI];
    g=[g; comp*ones(size(estMI))];
    stdMI=[stdMI; std(estMI)];
    
    comp=comp+1;
    colorsBox=[colorsBox; [0 0 1]];

    
    for ll=1:length(lengthWindScan)
    N=PAM_MI{1,ll}.N;
    MI=PAM_MI{1,ll}.MI;
    
    NDecimateFactor=floor(length(t_tot)/N);
    NDecimatedStart=floor(StepStart/NDecimateFactor)+1;% sp1=StepStart(1)/(Fs*60)

    NDecimatedStop =floor(StepEnd/NDecimateFactor);
    
    x=[x;MI(NDecimatedStart(indexS):NDecimatedStop (indexS))'];
    g=[g ; comp*ones(length(NDecimatedStart(indexS):NDecimatedStop (indexS)),1)];
    stdMI=[stdMI; std(MI(NDecimatedStart(indexS):NDecimatedStop (indexS))')]
    comp=comp+1;
    colorsBox=[colorsBox; colorFilt(ll,:)];
    length(lengthWindScan)*(indexS+2)+ll
    end
end


boxplot(x,g, 'color', colorsBox,'symbol', '')
plot(1:comp,mean(MI_th)*ones(1,comp),'-.' , 'color', 'k','linewidth',1.3)
ylim([0 0.7])
boxTickPos=(1:length(StepIndexes))*3-1;
set(gca,'xtick',boxTickPos,'xticklabels',tickval(2:end)   );

% BoxPlot Phi estimate
subplot(2,2,4); hold on
x=[];g=[];colorsBox=[];comp=1;
stdPhi=[];
for indexS=1:length(StepIndexes)

    estPhi=Phi_mod_t_k(1,(DecimatedStart(indexS):DecimatedStop(indexS)));
    x=[x;estPhi'];
    g=[g; comp*ones(size(estPhi'))];
    stdPhi=[stdPhi; std(estPhi)]
    comp=comp+1;
    colorsBox=[colorsBox; [0 0 1]];
    
    
    for ll=1:length(lengthWindScan)
        
    N=PAM_MI{1,ll}.N;
    MI=PAM_MI{1,ll}.MI;
    PAMusual=PAM_MI{1,ll}.PAMusual;
    lengthWindCur=lengthWindScan(ll);
    [~,locM]=max(PAMusual);
    phi_est_filt=-pbins2(locM);

    
    NDecimateFactor=floor(length(t_tot)/N);
    NDecimatedStart=floor(StepStart/NDecimateFactor)+1;
    NDecimatedStop =floor(StepEnd/NDecimateFactor);
    
    x=[x;phi_est_filt(NDecimatedStart(indexS):NDecimatedStop (indexS))'];
    g=[g ; comp*ones(length(NDecimatedStart(indexS):NDecimatedStop (indexS)),1)];
    stdPhi=[stdPhi; std(phi_est_filt(NDecimatedStart(indexS):NDecimatedStop (indexS)))]
    comp=comp+1;
    colorsBox=[colorsBox;colorFilt(ll,:)];
    length(lengthWindScan)*(indexS+2)+ll
    end
end


boxplot(x,g, 'color', colorsBox,'symbol', '')
plot(1:comp,mean(Phi_th)*ones(1,comp),'-.' , 'color', 'k','linewidth',1.3)
set(gca,'YTick',[-pi 0 pi/3 pi],'YTickLabel', {'-\pi','0','\pi/3', '\pi'})
set(gca,'xtick',boxTickPos,'xticklabels',tickval(2:end)   );
%xlabel('Noise Factor')



set(gcf,'Position',[680   462   750   350])
pathDrop='/home/nsrlhs/Dropbox (MIT)/producedFig/';
print([pathDrop,'VarianceEstimates'],'-dpdf', '-painters') 
