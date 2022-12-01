
function [beta_hat] = sample_from_pdf(X1,X2,pdf,Nsample,res)
%% Draw samples from a 2D discretized probability density function 
% 
% Inputs :  X1, X2 the points used to estimate pdf
%           pdf, the discrete probability density function 
%           Nsample number of drawn samples 
%           res >=1. Parameters to increase the resolution of pdf estimates
%           through cubic interpolatation (if >1). Adapted by hs from Tristan
%           Ursell (2012).
%
% Outputs : beta_hat ~ 2 x Nsample drawn samples

% pdf = p(beta1,beta2).  pdf_1 = p(beta1) 
pdf_1 = sum(pdf,2);

% Augment resolution through interpolation if res>1
if res ==1
    X1_aug = X1;
    X2_aug = X2;
    pdf_1_aug = pdf_1';
elseif res>1
    X1_aug=linspace(min(X1),max(X1),round(res*length(X1)));
    X2_aug=linspace(min(X2),max(X2),round(res*length(X2)));
    % 'abs()' makes sure the interpolated values are positive
    pdf_1_aug =abs(interp1(X1,pdf_1,X1_aug,'pchip')); 
    % Renormalize interpolated pdf
    pdf_1_aug = pdf_1_aug./sum(pdf_1_aug);
else 
  error('Incorrect res value.')
end

% Sample beta_1 (~ Nsamples) from p(beta1)
beta_1_hat_id = gen_1d_pdf(pdf_1_aug,Nsample);
beta_1_hat    = X1_aug(beta_1_hat_id);

% Initialize pdf_2 = p(beta2 | beta_1_hat) ~ Nsamples x lenght(X2)
pdf_2 = zeros(Nsample,length(X2));

% Map back beta_1_hat to X1 and grasp p(beta_1_hat,b2) ~ Nsamples x lenght(X2)
[val_temp,ind_temp]=sort((beta_1_hat'-X1).^2,2);

% If beta_1_hat is 
%in X1, no interpolation is needed for p(beta2 | beta_1_hat)
noInterpolation = val_temp(:,1)<eps;
noInterpolationId = find(noInterpolation);
pdf_2(noInterpolationId,:) = pdf(ind_temp(noInterpolationId,1),:);

% If beta_1_hat is in X1_aug\X1 we linearly interpolate p(beta2 | beta_1_hat)
doInterpolation = 1-noInterpolation;
doInterpolationId = find(doInterpolation); 

% Closest lower and Upper value in X1 
low_val_id  = min(ind_temp(doInterpolationId,1:2), [],2);
high_val_id = max(ind_temp(doInterpolationId,1:2), [],2);

X1_low  = X1(low_val_id)';
X1_high = X1(high_val_id)';
X1_hat_tmp = beta_1_hat(doInterpolationId)';

% Weights for the linear interpolation
weight_low  = 1-(X1_hat_tmp -X1_low)./(X1_high-X1_low);
weight_high = 1-(X1_high-X1_hat_tmp)./(X1_high-X1_low);
pdf_2_low =  pdf(low_val_id,:);
pdf_2_high = pdf(high_val_id,:);

pdf_2(doInterpolationId,:) = weight_low.*pdf_2_low + weight_high.*pdf_2_high;

% Build p(beta2 | beta_1_hat) ~ Nsample x {res x length(X1)} 
% (Augment resolution through interpolation if res>1)
if res == 1
    normC = sum(pdf_2,2);
    pdf_2_aug = (pdf_2./normC);

elseif res>1
    % 'abs()' makes sure the interpolated values are positive
    pdf_2_aug =abs(interp1(X2,pdf_2',X2_aug,'pchip')); 
    
     % Renormalize interpolated pdf
    normC = sum(pdf_2_aug,1);
    pdf_2_aug = (pdf_2_aug./normC)'; 
end


% Sample beta_1 (~ Nsamplesnon_norm_mvt_pdf = non_norm_mvt_pdf;) from p(beta2 | beta_1_hat)
beta_2_hat_id = gen_1d_pdf(pdf_2_aug,1);
beta_2_hat    = X2_aug(beta_2_hat_id);

beta_hat=[beta_1_hat;beta_2_hat];

end

function ids = gen_1d_pdf(pdf,Nsamples)
% Generate Nsamples from a 1D probability density function pdf
% Outputs ids, the id of the samples.

% pdf ~ ii x jj are ii probability density function estimated at points jj  
ii = size(pdf,1);
jj = size(pdf,2);

% We either draw many sample from same pdf or draw 1 sample from multiple pdf
assert(ii==1 || Nsamples ==1)

if min(min(pdf))<0
    error('Error:  All elements of first argument, P, must be positive.')
end

%re-normalize pdf
pdf1d_norm=[zeros(ii,1) pdf]./sum(pdf,2);

%create cumlative density function
cdf=cumsum(pdf1d_norm,2);

% Sample from Uniform distribution 
U=rand(max(Nsamples,ii),1);

if ii ==1
    U = repmat(U, 1,jj+1);
end

% Find the smallest bin satisfying cdf > U.
binbelong = sum(U>=cdf,2);% assert size X x jj;
ids = binbelong;

end




