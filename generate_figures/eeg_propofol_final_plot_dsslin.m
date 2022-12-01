%% Define Paths
addpath(genpath('../'))
bhvr_path   = '/media/hdhs/data/eeg_propofol_data/behavioral_data/';
decomp_path ='/media/hdhs/data/eeg_propofol_data/eeg_init/';
linreg_path ='/media/hdhs/data/eeg_propofol_data/dss_lin_pac_3D/';
subjectId='04';


%% Load and Grasp behavioral data
name_bhvr = [bhvr_path, 'eeganes_bhvr_summary_sub',subjectId];
load([name_bhvr,'.mat'])

tprop_min = bhvr_summary.prop_time/(bhvr_summary.Fs_behav*60);
tburs_min = bhvr_summary.t_burst/(bhvr_summary.Fs_behav*60);
tverb_min = bhvr_summary.t_verb/(bhvr_summary.Fs_behav*60);

%% Load and Grasp Decomposed Data
em_its=400;
lengthWind_sec=6;
name_decomp=['pam_eeg', subjectId];
load([decomp_path,name_decomp,'_',num2str(lengthWind_sec),'s_', num2str(em_its), 'ite_autoinit','.mat']);


%% Load and Grasp dsslin processed data
name_linreg = [linreg_path, 'dss_lin_pac_new_ss_filtering_sub',subjectId];
load([name_linreg,'.mat'])








%% Compute Parametric Spectrum

startPointId=modelOsc.startPointId;
endPointId=modelOsc.endPointId;

tVect_sec = 0.5*(endPointId+endPointId)/Fs;
f_y=0.1:0.1:120;
[f_tot, a_tot, q_tot, R_tot, s_tot_param, s_param , Nosc_ind,NNharm, s_mt,s_mt_tot, fmt] = gather_params(modelOsc,f_y);



%% Compute Multitaper Spectrum and Standard PAM
Fs = modelOsc.Fs;
specparams=struct();
specparams.Fs = Fs;
specparams.tapers = [3 5];
specparams.trialave = 1;
specparams.fpass=[0 120];

y_tot=zeros(1, endPointId(1,end)-startPointId(1,1) );
for tt=1:length(startPointId)
    x_t_n_tot(:,startPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.x_t_n(:,2:end);
    y_tot    (1,startPointId(1,tt):endPointId(1,tt))=modelOsc.res{1,tt}.y;
end
[S_mt, stimes_mt, sfreqs_mt]= mtspecgramc(y_tot,[5 5],specparams);


Npb=18;               % Phase Bin Number
pbins=linspace(-pi,pi,Npb+1); % Phase bins
lengthWind_PAM=120;
lengthWind=lengthWind_PAM*Fs;
N=floor(length(y_tot)/lengthWind);
cutMeas=y_tot(1: N*lengthWind);
resMeas=reshape(cutMeas, [lengthWind, N]);
pac_standard=zeros(Npb, N);
sofreq=[0.1 1]; alfreq=[8 15]; orderfilt=100;

for tt=1:N
    disp([num2str(tt),'/', num2str(N) ])
    
    [x_so_tmp, tail1] = quickbandpass(resMeas(:,tt)', Fs, sofreq,orderfilt);
    [x_al_tmp, tail2] = quickbandpass(resMeas(:,tt)', Fs, alfreq,orderfilt);
    
    Phi_filt_so_tmp=angle(hilbert(x_so_tmp));
    A_filt_al_tmp  =abs(hilbert(x_al_tmp));
    
    pam_filt_tmp=phaseamp(A_filt_al_tmp ,Phi_filt_so_tmp,pbins);
    pac_standard(:,tt)=pam_filt_tmp;
end




%% Figure 1 : Select Window of Interest and illustrate Osc decomposition

scaleTime=1;
scalemuVolt=20;

WOI=904;
tspan=startPointId(WOI):endPointId(WOI);
ta=tspan/(Fs*60);
y_t0=y_tot(1,tspan);

%x_so=X_t_so_final(tspan, :);
%x_al=X_t_al_final(tspan, :);

modelOsc_1wind = modelOsc;
modelOsc_1wind.startPointId = 1;
modelOsc_1wind.endPointId   = length(tspan);
modelOsc_1wind.res = cell(1,1);
modelOsc_1wind.res{1,1} = modelOsc.res{1,WOI};
% plot_summary(modelOsc_1wind)

Kappa1=dsslin_pac.Kappa1;
Kappa2=dsslin_pac.Kappa2;
Kappa = [Kappa1;Kappa2];
slowID = 1; fastID = 2; Nlinconst = 100;Nlag =inf;

[latent_pac, SS_tot,X_t_tot] =ssp_pac(modelOsc_1wind, Kappa, slowID, fastID, Nlag,Nlinconst);

%%
x_so_tot = reshape(X_t_tot(1:2,1:end-1,1,:), 2,length(tspan), Kappa1+1);
x_al_tot = reshape(X_t_tot(3:4,1:end-1,1,:), 2,length(tspan), Kappa1+1);

x_so = reshape(x_so_tot(1,:,:),length(tspan), Kappa1+1);
x_al = reshape(x_al_tot(1,:,:),length(tspan), Kappa1+1);

phi_so_tot =  reshape(atan2(x_so_tot(2,:,:),x_so_tot(1,:,:)),length(tspan), Kappa1+1);
phi_al_tot =  reshape(atan2(x_al_tot(2,:,:),x_al_tot(1,:,:)),length(tspan), Kappa1+1);
amp_al_tot =  reshape(sqrt(x_al_tot(1,:,:).^2+x_al_tot(2,:,:).^2),length(tspan), Kappa1+1);

amp_al_sup = max(amp_al_tot,[],2);
amp_al_inf = min(amp_al_tot,[],2);
[phi_so_mean,phi_so_sup,phi_so_inf] = get_phase_msi(phi_so_tot');

K_mod = sqrt(latent_pac(3,:).^2+latent_pac(2,:).^2)./latent_pac(1,:);
P_mod = atan2(latent_pac(3,:),latent_pac(2,:));


%%

M_t = @(Phi) [ones(length(Phi),1) cos(Phi), sin(Phi)];

A_hat_tot = zeros(size(phi_so_tot,1),size(latent_pac,2));
compt=1;

for kk1 = 1:Kappa1+1
    X_t = M_t(phi_so_tot(:,kk1));
    
    for kk2 = 1:Kappa2+1
    curK = kk2 + (kk1-1)*(Kappa1+1);
    A_hat = X_t * latent_pac(:,curK);
    
    A_hat_tot(:,curK) = A_hat;
    end
    
end

A_hat_mean = mean(A_hat_tot,2);
A_hat_inf  = min(A_hat_tot,[],2);
A_hat_sup  = max(A_hat_tot,[],2);


%%

figmain= figure; cla;
% Plot total raw signal
subplot(4,1,1);cla; hold on;
plot((1:length(y_tot))/(Fs*60),y_tot, 'color', 'k')
xlabel('(min)')
ylabel('\muV')
box on

% Plot raw isolated window
pos1=get(subplot(4,3,4), 'Position');
pos2=get(subplot(4,3,5), 'Position');
width1=pos1(3);heigh1=pos1(4);
ypos =pos1(2);xpos =(pos2(3)+pos2(1)+pos1(1)-width1)/2;
posiRaw=[xpos, ypos, width1, heigh1];
posiRaw=get(subplot(4,3,4), 'Position'); subplot(4,2,4);set(gca, 'Visible', 'off')
f1=subplot('Position', posiRaw); hold on;
plot(ta,y_t0, 'color', 'k')
set(gca,'xtick',[],'ytick',[])



% Plot alpha
f3=subplot(4,3,6);cla; hold on;
plot(ta, x_so(:,1), 'Color', [0.7 0.7 0.7])
plot(ta, x_al(:,1), 'Color', 'r')

% Plot slow
f2=subplot(4,3,5);cla; hold on;
plot(ta, x_al(:,1), 'Color', [0.7 0.7 0.7])
plot(ta, x_so(:,1), 'Color', 'b')
linkaxes([f1 f2 f3])
axis tight

% Plot scales
subplot('Position', posiRaw);
axis tight
box on
ylimTmp=ylim;
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])



% Plot scales
subplot(4,3,5)
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10 ,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')
subplot(4,3,6)
plot([ta(1); ta(1)], [min(x_so(:,1)); min(x_so(:,1))+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [min(x_so(:,1)); min(x_so(:,1))], '-k', 'LineWidth', 2)
text(ta(1),min(x_so(:,1))-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')


% Plot Alpha amplitude 
subplot(4,3,9);cla; hold on;
ciplot(amp_al_inf,amp_al_sup,ta,[0.5 0.5 0.5])
plot(ta, amp_al_tot(:,1), 'Color', 'r')
ylimTmp=ylim;
plot([ta(1); ta(1)], [0; 0+scalemuVolt], '-k',  [ta(1); ta(scaleTime*Fs)], [0;0], '-k', 'LineWidth', 2)
text(ta(1),0-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt), '\muV'])
set(gca, 'Visible', 'off')

% Plot Slow Phase
subplot(4,3,8);cla; hold on;
plot_ci_phase_hs(phi_so_sup,phi_so_inf, ta,[0.5 0.5 0.5])
scatter(ta, phi_so_tot(:,1),1,'b')

% Plot Scales
axis tight
ylim([-pi pi])
set(gca,'xtick',[],'ytick',[-pi,0, pi],'YTickLabel',{'-\pi','0','\pi'}, 'FontSize', 16)
box on
ylimTmp=ylim;
plot([ta(1); ta(scaleTime*Fs)], [-pi;-pi], '-k', 'LineWidth', 2)
text(ta(1),-pi-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s'])






zoomP = 2;
zoomK = 1;

K_mean = mean(K_mod);
P_mean = get_phase_msi(P_mod');


subplot(4,3,10)
Hk = histogram(K_mod, 'FaceColor', [0.75 0.75 0.75],'edgecolor', 'none', 'normalization', 'probability');
xlim([0 1])


[~,edgedIdK] = min(abs(Hk.BinEdges-K_mean));
x = Hk.BinEdges(edgedIdK:edgedIdK+1); 
dx = diff(x)*zoomK;
x(1,1) = x(1,1)-dx; x(1,2) = x(1,2)+dx;
x = [x fliplr(x)]; y = [0 0 repmat(Hk.Values(edgedIdK), 1, 2)]; 
S=patch(x, y,'g');
S.FaceColor = 'g';
S.EdgeColor = 'none';

subplot(4,3,11)
Hp = histogram(P_mod, 'FaceColor', [0.75 0.75 0.75],'edgecolor', 'none', 'normalization', 'probability');
xlim([-pi pi])
[~,edgedIdP] = min(abs(Hp.BinEdges-P_mean));
x = Hp.BinEdges(edgedIdP:edgedIdP+1); 
dx = diff(x)*zoomP;
x(1,1) = x(1,1)-dx; x(1,2) = x(1,2)+dx;
x = [x fliplr(x)]; y = [0 0 repmat(Hp.Values(edgedIdP), 1, 2)]; 
S=patch(x, y,'g');
S.FaceColor = 'g';
S.EdgeColor = 'none';



subplot(4,3,12)
ciplot(A_hat_inf,A_hat_sup,ta,[0.5 0.5 0.5]); hold on
plot(ta,A_hat_mean, 'color', 'g', 'LineWidth', 2)
scalemuVolt2 = 10;

ylimTmp=ylim;
plot([ta(1); ta(1)], [0; 0+scalemuVolt2], '-k',  [ta(1); ta(scaleTime*Fs)], [0;0], '-k', 'LineWidth', 2)
text(ta(1),0-(ylimTmp(2)-ylimTmp(1))/10,[num2str(scaleTime),'s/', num2str(scalemuVolt2), '\muV'])
set(gca, 'Visible', 'off')



figmain.Renderer='Painters';
set(figmain,'Position',[356   210   816   832]);
pathDrop='';
%print(gcf, '-dpdf', [pathDrop,'fig1_aug'])

%% Figure 2

figmain=figure; hold on


% Propofol Concentration
subplot(4,3,4)
plot(tprop_min,bhvr_summary.prop_conc,'Color', 'k')
ylabel('ug/ mL')
box on
xlim([tprop_min(1) tprop_min(end)])

% Probability of response
subplot(4,3,5)

%ciplot(bhvr_summary.p05_burst,bhvr_summary.p95_burst,tburs_min,[0 0 1], 0.5); hold on
%p1=plot(tburs_min,bhvr_summary.p_burst, 'linewidth',2, 'color' , 'b' );

ciplot(bhvr_summary.p05_verb,bhvr_summary.p95_verb,tverb_min,[0.5 0.5 0.5],0.5); hold on
p2=plot(tverb_min,bhvr_summary.p_verb, 'linewidth',2, 'color' , 'k' );
ylabel('P(response)')
%legend([p1,p2], {'Burst', 'Verbal'})
box on ;xlim([tverb_min(1) tverb_min(end)])

subplot(4,3,6); hold on
%ciplot(R2_inf,R2_sup,tVect,[0.5 0.5 0.5])
plot(tVect_sec/(60),dsslin_pac.R2_mis(1,:), 'k')
box on
xlim([tVect_sec(1)/(60) tVect_sec(end)/(60)])
ylim([0 1])
ylabel('R^2')


% Spectrogram
cmax=20;
cmin=-15;
p1=subplot(2,3,4);
imagesc(tVect_sec/(60), f_y, 10*log10(s_tot_param))
colormap(p1,jet)
axis xy;
ylabel('Freq [Hz]')
ylim([0 30])
xlabel(['Time [min]'])
caxis([cmin cmax])
cc = colorbar( 'location', 'northoutside');
cc.Label.String = 'dB';

% Modulogram
c0m = pi/Npb;
cr = 0.8*pi/Npb;
cmmin = c0m-cr;
cmmax = c0m+0.9*cr;

p2 =subplot(2,3,5);hold on
 imagesc(tVect_sec/(60) ,dsslin_pac.pbins,dsslin_pac.PAC_filt);
col = redbluemap; colormap(p2,col);
cc = colorbar( 'location', 'northoutside');
cc.Label.String = 'Alpha Power Distribution';
xlim([tVect_sec(1)/(60) tVect_sec(end)/(60)])
set(gca, 'clim', [cmmin cmmax]);
axis tight
ylabel('DSS Lin Filt', 'FontWeight', 'bold');
box on ;
xlabel('[min]')



K_mean = dsslin_pac.K_mis_filt(1,:); 
if size(dsslin_pac.K_mis,1) == 3
    K_inf = dsslin_pac.K_mis_filt(2,:);
    K_sup = dsslin_pac.K_mis_filt(3,:);
end

Phi_mean = dsslin_pac.Phi_mis_filt(1,:); 
if size(dsslin_pac.R2_mis,1) == 3
    Phi_inf = dsslin_pac.Phi_mis_filt(2,:);
    Phi_sup = dsslin_pac.Phi_mis_filt(3,:);
end
tVect =dsslin_pac.tVect;

% K_mod
subplot(4,3,9);hold on
ciplot(K_inf,K_sup,tVect,[0.5 0.5 0.5])
l1=plot(tVect,K_sup,'Color',[0.5 0.5 0.5]);
l2=plot(tVect,K_mean,'g');
ylim([0 1])
xlim([tVect(1) tVect(end)])
box on
legend([l1,l2], {[num2str(floor(100*(1-dsslin_pac.alpha_ci))), '% CI',num2str(Kappa1),'\times', num2str(Kappa2)],'Mean'})
ylabel(['K^{mod}_t. '])

% Phi_mod
subplot(4,3,12);hold on
plot_ci_phase_hs(Phi_sup,Phi_inf, tVect,[0.5 0.5 0.5])
scatter(tVect,Phi_mean,10,'g','filled')
ylim([-pi, pi])
xlim([tVect(1) tVect(end)])
box on
xlabel('[min]')
ylabel(['\phi^{mod}_t. '])




figmain.Renderer='Painters';
set(figmain,'Position',[207         446        1033         500]);
pathDrop='';
orient(figmain,'landscape')
print(figmain, '-dpdf', [pathDrop,'fig2_aug'])





