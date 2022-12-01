dpath = '/media/hdhs/data/snr_dss_vs_standard/';
load([dpath, 'snr_dss_vs_standard_simulated_signal'])
load([dpath, 'snr_dss_vs_standard_PAC.mat'])
load([dpath, 'snr_dss_vs_standard_modelOsc.mat'])


%% Filter Windows500
sofreq=[0.1 1.5]; alfreq=[8 12]; 
lengthWindScan= [6  120];
orderfiltScan = [50 200];
Npb=18;
PAC_MI_standard=cell(1,length(lengthWindScan));

pbins=linspace(-pi,pi,Npb+1);
pbins2=linspace(-pi,pi,Npb);

for ll=1:length(lengthWindScan)
    lengthWindCur=lengthWindScan(ll);
    PAC_MI_standard{1,ll}.lengthWindScan=lengthWindScan(ll);
    

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
    PAC_MI_standard{1,ll}.N=Ncur;
    PAC_MI_standard{1,ll}.MI= MIusual_cur;
    PAC_MI_standard{1,ll}.PAMusual=PAMusual;
    
end



%% Equiv MI pAM
uni=ones(Npb,1);
dPhi=2*pi/Npb;

MtPhi_th= (1 + 2*K_th .*cos(pbins2'+Phi_th) *sin(dPhi/2)/(2*pi/Npb));
MI_th  = sum(abs(MtPhi_th-uni),1)/Npb;

PAC_DSS=zeros(Npb, size(K_mod_t_0,2), size(K_mod_t_0,1));
MI_DSS=zeros(size(K_mod_t_0,2), size(K_mod_t_0,1));
for kk=1:size(K_mod_t_0,1)
    PAC_DSS(:,:,kk)=(  1 + 2*K_mod_t_0(kk,:).*cos(pbins2'+Phi_mod_t_0(kk,:))*sin(dPhi/2)/(2*pi/Npb));
    MI_DSS(:,kk)= sum(abs(PAC_DSS(:,:,kk)-uni),1)/Npb;
end




%% Plots
DecimateFactor=1;
colorFilt=[[0.8 0.8 0.8];[0.4 0.4 0.4];[0.8 0 0]]; lineW= [1 , 1.3, 1.3];


figure;

% Plot the Evolution of MI estimate
subplot(2,2,1); hold on;
for ll=1:length(lengthWindScan)
    N=PAC_MI_standard{1,ll}.N;
    MI=PAC_MI_standard{1,ll}.MI;
    PAMusual=PAC_MI_standard{1,ll}.PAMusual;
    lengthWindCur=lengthWindScan(ll);
    plot((1:N)*lengthWindCur/60,MI, 'color', colorFilt(ll,:),'linewidth',lineW(ll));
end
plot(t_tot(1:DecimateFactor:end)/(60*Fs),MI_DSS(:,1), 'color', 'b','linewidth',1.3);
plot(t_tot/(60*Fs),MI_th,'-.' , 'color', 'k','linewidth',1.1)
ylim([0 0.7])
box on

   
% Plot the Evolution of Phi estimate
subplot(2,2,2); hold on;
for ll=1:length(lengthWindScan)
    N=PAC_MI_standard{1,ll}.N;
    MI=PAC_MI_standard{1,ll}.MI;
    PAMusual=PAC_MI_standard{1,ll}.PAMusual;
    lengthWindCur=lengthWindScan(ll);
    [~,locM]=max(PAMusual);
    phi_est_filt=pbins2(locM);
    plot((1:N)*lengthWindCur/60,-phi_est_filt, 'color', colorFilt(ll,:),'linewidth',lineW(ll));
end
plot(t_tot(1:DecimateFactor:end)/(Fs*60), Phi_mod_t_0(1,:), 'color', 'b', 'linewidth',1.2)
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
tickval=floor(unique(SNR));
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
tickval=floor(unique(SNR));
set(gca,'xtick',tickPos,'xticklabels', tickval(2:end) );
xlim([t_tot(1)/(60*Fs) t_tot(end)/(60*Fs)])
box on





% BoxPlot MIi estimate
subplot(2,2,3); hold on
x=[];g=[];colorsBox=[];comp=1;
stdMI=[];

for indexS=1:length(StepIndexes)

    estMI=MI_DSS(DecimatedStart(indexS):DecimatedStop(indexS),1);
    x=[x;estMI];
    g=[g; comp*ones(size(estMI))];
    stdMI=[stdMI; std(estMI)];
    
    comp=comp+1;
    colorsBox=[colorsBox; [0 0 1]];

    
    for ll=1:length(lengthWindScan)
    N=PAC_MI_standard{1,ll}.N;
    MI=PAC_MI_standard{1,ll}.MI;
    
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

    estPhi=Phi_mod_t_0(1,(DecimatedStart(indexS):DecimatedStop(indexS)));
    x=[x;estPhi'];
    g=[g; comp*ones(size(estPhi'))];
    stdPhi=[stdPhi; std(estPhi)]
    comp=comp+1;
    colorsBox=[colorsBox; [0 0 1]];
    
    
    for ll=1:length(lengthWindScan)
        
    N=PAC_MI_standard{1,ll}.N;
    MI=PAC_MI_standard{1,ll}.MI;
    PAMusual=PAC_MI_standard{1,ll}.PAMusual;
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
%pathDrop='/home/nsrlhs/Dropbox (MIT)/producedFig/';
%print([pathDrop,'VarianceEstimates'],'-dpdf', '-painters') 