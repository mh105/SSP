addpath(genpath('../dss_decomp'))
addpath(genpath('../spectralanalysis'))
addpath(genpath('../dss_PAM_bootstrap'))
addpath(genpath('../usual_pam'))


%% Load Data
dpath='/media/nsrlhs/f0e92df7-6f0d-4f06-bbeb-b6e03ea748fa/data/eeg_propofol_data/';
dpath2='/media/nsrlhs/f0e92df7-6f0d-4f06-bbeb-b6e03ea748fa/data/rawEEg/';

DecimateFactor=20;
PAM_MI=cell(1,2);

% Subject A
load([dpath,'pam_eeg04_bootstrap2analysis300ite.mat'])
PAM_MI{1,1}.name='pam_eeg04';
PAM_MI{1,1}.K_mod  = K_mod_t_k(:,:);
PAM_MI{1,1}.Phi_mod= Phi_mod_t_k(:,:);

load([dpath2,'eeganes04_laplac250_ch36.mat'],'BHVR', 'DRUG')
PAM_MI{1,1}.DRUG=DRUG;
PAM_MI{1,1}.BHVR=BHVR;
clear('K_mod_t_k','Phi_mod_t_k','BHVR', 'DRUG' )


% Subject B
load([dpath,'pam_eeg03_bootstrap2analysis300ite.mat'])
PAM_MI{1,2}.name='pam_eeg03';
PAM_MI{1,2}.K_mod  =  K_mod_t_k(:,:);
PAM_MI{1,2}.Phi_mod=Phi_mod_t_k(:,:);
clear('K_mod_t_k','Phi_mod_t_k')

load([dpath2,'eeganes03_laplac250_ch36.mat'],'BHVR', 'DRUG')
PAM_MI{1,2}.DRUG=DRUG;
PAM_MI{1,2}.BHVR=BHVR;
clear('K_mod_t_k','Phi_mod_t_k','BHVR', 'DRUG' )

%% Parametric MI
Npb=18; theta=linspace(-pi, pi,Npb);pbins=linspace(-pi,pi,Npb+1);

for subj=1:length(PAM_MI)
    
    % Parametric PAM
    K_mod  = PAM_MI{1,subj}.K_mod;
    Phi_mod= PAM_MI{1,subj}.Phi_mod;
    TT     = size( K_mod,2);
    Kappa  = size(K_mod,1)-1;
    MI_tot = zeros(Kappa+1, TT);
    
    
    for kk=1:Kappa+1
        disp(['Subject : ',num2str(subj), '/', num2str(length(PAM_MI)),' Kappa : ',num2str(kk), '/', num2str(Kappa+1)])
        MI_tmp = sum(abs( 2*sin(pi/Npb)*K_mod(kk,:).*cos(Phi_mod(kk,:)-theta')*Npb/(2*pi)),1);
        MI_tot(kk,:)=MI_tmp /Npb;
    end
    
    PAMequiv =  (1+2*sin(pi/Npb)*K_mod(1,:).*cos(Phi_mod(1,:)+theta')*Npb/(2*pi));
    PAM_MI{1,subj}.MI = MI_tot(1,:);
    PAM_MI{1,subj}.PAM= PAMequiv;
    PAM_MI{1,subj}.TT=TT;

    % Retrieve params from modelOsc
    load([dpath,PAM_MI{1,subj}.name,'_6s_300ite.mat']);
    startPointId=modelOsc.startPointId;
    endPointId=modelOsc.endPointId;
    timeWindNum=length(startPointId);
    Fs=modelOsc.Fs; T=length(startPointId);
    y_tot=zeros(1, endPointId(1,end)-startPointId(1,1) );
    for tt=1:T
        disp(['Subject : ',num2str(subj), '/', num2str(length(PAM_MI)),' Window :',num2str(tt), '/', num2str(T)])
        y_tot    (1,startPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.y;
    end
    clear('modelOsc');
    
    
    
    % Usual PAM for different dt and different bandpass cutoff
    lengthWindScan= [6  120];
    orderfiltScan = [50 300];
    sofreqTot=[[0.1 1];[0.1 1]] ;
    alfreqTot=[[8 12] ; [8 14]] ;
    
    sofreqTot=[[0.1 1]] ;
    alfreqTot=[[8 14]] ;
    
    sofreq=[0.1 1]; alfreq=[13];
    
    
    
    PAM_MI{1,subj}.PAM_MI_usual=cell(1,length(lengthWindScan));
   PAM_MI{1,subj}.lengthWindScan=lengthWindScan;
    
    
    
    for ll=1:length(lengthWindScan)
        lengthWindCur=lengthWindScan(ll);
        PAM_MI{1,subj}.PAM_MI_usual{1,ll}.lengthWindScan=lengthWindScan(ll);
        PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt=cell(1,size(sofreqTot,1));
    
        for CutOffID=1:size(sofreqTot,1)

            sofreq=sofreqTot(CutOffID,:);
            alfreq=alfreqTot(CutOffID,:);
            PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt{1,CutOffID}.sofreq=sofreq;
            PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt{1,CutOffID}.alfreq=alfreq;
            
            
            
            lengthWind=lengthWindCur*Fs;
            Ncur=floor(length(y_tot)/lengthWind);
            cutMeas=y_tot(1: Ncur*lengthWind);
            resMeas=reshape(cutMeas, [lengthWind, Ncur]);
            PAMusual=zeros(Npb, Ncur);
            
            for tt=1:Ncur
                disp(['Subject : ',num2str(subj), '/', num2str(length(PAM_MI)),' PAM, frq: ' ,num2str(CutOffID),'/', num2str(size(sofreqTot,1)), ' length :' ,num2str(ll),'/', num2str(length(lengthWindScan)),' || ' ,num2str(tt),'/', num2str(Ncur) ])
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
            PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt{1,CutOffID}.N=Ncur;
            PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt{1,CutOffID}.MI= MIusual_cur;
            PAM_MI{1,subj}.PAM_MI_usual{1,ll}.PAM_MI_usual_filt{1,CutOffID}.PAMusual=PAMusual;
            
        end
    end
  
end

%save('./DataForfigureMIfinal')

%% Plots : PAM


Ynum=size(PAM_MI,2);
WindLengthN=size(PAM_MI{1,1}.PAM_MI_usual,2);
Xnum=WindLengthN+2;


xlimTot=[[20 150]; [40 140]];
climTot=[[0.2 1.8]; [0.2 1.8]];


figure;

for subj=1:Ynum
    
    
    IdPlot=(1:Ynum:Ynum*Xnum)+(subj-1);
    
    % Response Probability
    pProb=subplot(Xnum,Ynum,IdPlot(1,1)); hold on; cla;
    DRUG=PAM_MI{1,subj}.DRUG; FsMAX=DRUG.Fs;
    BHVR=PAM_MI{1,subj}.BHVR;
    pr=plot(BHVR.T'/(FsMAX*60),BHVR.prob_verbal.p500 ,'color','r','linewidth',2);
    jbfill(BHVR.T'/(FsMAX*60),BHVR.prob_verbal.p975',BHVR.prob_verbal.p025','r');
    hold on;
    pb=plot(BHVR.T'/(FsMAX*60),BHVR.prob_burst.p500 ,'color','b','linewidth',2);
    jbfill(BHVR.T'/(FsMAX*60),BHVR.prob_burst.p975',BHVR.prob_burst.p025','b');
    if subj==Ynum;legend([pr,pb],{'Verbal','Click'},'Location','North');end;
    title('Response Probability')
    xlim([xlimTot(subj,1),xlimTot(subj,2)])
    box on
    
    % Parameset(gca, 'clim', [0.2 1.8]); %0.2-1.8tric PAM
    pMI=subplot(Xnum,Ynum,IdPlot(1,end)); hold on; cla;
    PAMequiv = PAM_MI{1,subj}.PAM;
    TT=PAM_MI{1,subj}.TT;
    imagesc((1:TT)*DecimateFactor/(Fs*60),-theta,PAMequiv);
    col = redbluemap; cax=colormap(col);
    axis tight
    xlabel('(min)')
    set(gca, 'clim', [climTot(subj,1),climTot(subj,2)]);
    set(gca,'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'})
    box on
    title('Parametric : 6s')
    xlim([xlimTot(subj,1),xlimTot(subj,2)])

    % Usual PAM
    for ll= 1:size(PAM_MI{1,subj}.PAM_MI_usual,2)
        subplot(Xnum,Ynum,IdPlot(1,ll+1)); hold on; cla;
        
        
        PAM_MI_cur_tmp=PAM_MI{1,subj}.PAM_MI_usual{1,ll};
        lengthWindScan=PAM_MI{1,subj}.lengthWindScan(ll);
        title(['Standard : ', num2str(lengthWindScan), 's'])

        
        for CutOffID=size(PAM_MI_cur_tmp.PAM_MI_usual_filt,2):size(PAM_MI_cur_tmp.PAM_MI_usual_filt,2)
            PAM_MI_cur=PAM_MI_cur_tmp.PAM_MI_usual_filt{1,CutOffID};
            N=PAM_MI_cur.N;
            PAMMIusual=PAM_MI_cur.PAMusual;

            imagesc((1:N)*lengthWindScan/60,-pbins,PAMMIusual);
            col = redbluemap; colormap(col);
            set(gca, 'clim', [climTot(subj,1),climTot(subj,2)]);
            axis tight
            set(gca,'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'})
            set(gca,'xtick',[])
            box on
            xlim([xlimTot(subj,1),xlimTot(subj,2)])
            
        end
    end
    
    
    
    % PAM colorbar
    pos1=get(subplot(Xnum,Ynum,1+Ynum), 'Position');
    pos2=get(subplot(Xnum,Ynum,(Xnum-1)*Ynum), 'Position');
    widtBar=pos2(1)/50;
    middleFig=(pos1(1)+pos2(1)+pos2(3))/2;%-widtBar/1000;
    highBar= abs(pos1(2)+pos1(4)-(pos2(2)));
    cax = colorbar('Position', [middleFig pos2(2)  widtBar highBar]); %caxis([0.5 1.5]);
    set(gca, 'clim', [climTot(subj,1),climTot(subj,2)]); %0.2-1.8
    set(cax  , 'Ticks', [0.5 1 1.5]);
    %ylabel(cax,'Alpha power Distribution');
    uu=colormap(cax,col);
    
    
end

%Save 
h=gcf;
set(h,'Position',[656   479   678   612]);
%print(gcf, '-dpdf', 'PAMsubjectAB.pdf')


%% Plots MI
Ynum=size(PAM_MI,2);
WindLengthN=size(PAM_MI{1,1}.PAM_MI_usual,2);
Xnum=WindLengthN+1;

xlimTot=[[20 150]; [40 140]];
climTot=[[0.2 1.8]; [0.2 1.8]];
xTickTot=[20 40 60 80 100 120 140 160];

FilterNum=size(PAM_MI{1,1}.PAM_MI_usual{1,1}.PAM_MI_usual_filt,2);
colorUsualFilt=[1 0 0]+ (((1:FilterNum)'-1)/FilterNum)*[-0.4 0.4 0.4];

figure;

for subj=1:Ynum
    
    
    IdPlot=(1:Ynum:Ynum*Xnum)+(subj-1);
    
   
    
    % Parameset(gca, 'clim', [0.2 1.8]); %0.2-1.8tric PAM
    pMI=subplot(Xnum,Ynum,IdPlot(1,end)); hold on; cla;
    MIequiv = PAM_MI{1,subj}.MI;
    TT=PAM_MI{1,subj}.TT; 
    plot((1:TT)*DecimateFactor/(Fs*60),MIequiv, 'color', 'b')
    axis tight
    xlabel('(min)')
    %set(gca,'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'})
    box on
    title('Parametric : 6s')
    xlim([xlimTot(subj,1),xlimTot(subj,2)])
    set(gca, 'XTick', xTickTot)

    
    %Usual MI
    for ll= 1:size(PAM_MI{1,subj}.PAM_MI_usual,2)
        subplot(Xnum,Ynum,IdPlot(1,ll)); hold on; cla;
        
        
        PAM_MI_cur_tmp=PAM_MI{1,subj}.PAM_MI_usual{1,ll};
        lengthWindScan=PAM_MI{1,subj}.lengthWindScan(ll);
        title(['Standard : ', num2str(lengthWindScan), 's'])

        
        for CutOffID=1:size(PAM_MI_cur_tmp.PAM_MI_usual_filt,2)
            PAM_MI_cur=PAM_MI_cur_tmp.PAM_MI_usual_filt{1,CutOffID};
            N=PAM_MI_cur.N;
            MIusual=PAM_MI_cur.MI;
            pb=plot((1:TT)*DecimateFactor/(Fs*60),MIequiv, '-', 'color', 'b', 'linewidth',1)
            pr=plot((1:N)*lengthWindScan/60, MIusual, 'color',colorUsualFilt(CutOffID,:), 'linewidth',1.5);
            
            if ll==1 && subj==Ynum ;legend([pr,pb],{'Usual','Parametric 6s'},'Location','North');end;
            axis tight
            %set(gca,'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'})
            box on
            xlim([xlimTot(subj,1),xlimTot(subj,2)])
            set(gca, 'XTick', xTickTot);
        end
    end
    
    
  
  
    
    
end

%Save 
h=gcf;

set(h,'Position',[421   259   678   449]);
print(gcf, '-dpdf', 'MIsubjectAB.pdf')

