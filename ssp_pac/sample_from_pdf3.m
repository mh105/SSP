function [beta_hat] = sample_from_pdf3(X1,X2,X3,pdf,Nsample,res)
%% Draws samples from a 3D discretized probability density function 
% 
% INPUTS :  -X1, X2, X3 points used to estimate pdf
%           -pdf, the discrete probability density function 
%           -Nsample number of drawn samples 
%           -res >=1. Parameters used to increase the resolution of pdf estimates
%           through cubic interpolatation (if >1). Adapted by hs from Tristan
%           Ursell (2012).
%
% OUTPUTS : beta_hat ~ 3 x Nsample drawn samples


if ndims(pdf) ~= 3
   error('Sampling method written for 3 dimmensional distributions') 
end

% pdf = p(x1,x2,x3). pdf_1 = p(x1)
pdf_1 = sum(sum(pdf,3),2);

% Augment resolution through interpolation if res>1
if res ==1
    X1_aug = X1;
    pdf_1_aug = pdf_1';
elseif res>1
    X1_aug=linspace(min(X1),max(X1),round(res*length(X1)));
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

% Grasp pdf_23 = p(beta1, x2,x3) ~ Nsamples x lenght(X2) x lenght(X3)
pdf_23 = zeros(Nsample,length(X2),length(X3) );

% Map back beta_1_hat to X1 and grasp p(beta_1_hat,b2) ~ Nsamples x lenght(X2)
[val_temp,ind_temp]=sort((beta_1_hat'-X1).^2,2);

% If beta_1_hat is in X1, no interpolation is needed for p(beta1, x2,x3)
noInterpolation = val_temp(:,1)<eps;
noInterpolationId = find(noInterpolation);

% If beta_1_hat is in X1_aug\X1 we linearly interpolate p(beta2 | beta_1_hat)
doInterpolation = 1-noInterpolation;
doInterpolationId = find(doInterpolation); 

% pdf when no interpolation is needed
pdf_23(noInterpolationId,:,:) = pdf(ind_temp(noInterpolationId,1),:,:);

% Closest lower and Upper value in X1 
low_val_id  = min(ind_temp(doInterpolationId,1:2), [],2);
high_val_id = max(ind_temp(doInterpolationId,1:2), [],2);

X1_low  = X1(low_val_id)';
X1_high = X1(high_val_id)';
X1_hat_tmp = beta_1_hat(doInterpolationId)';

% Weights for the linear interpolation
weight_low  = 1-(X1_hat_tmp -X1_low)./(X1_high-X1_low);
weight_high = 1-(X1_high-X1_hat_tmp)./(X1_high-X1_low);

pdf_23_low =  pdf(low_val_id,:,:);
pdf_23_high = pdf(high_val_id,:,:);

% Do linear interpolation of the pdf
pdf_23(doInterpolationId,:,:) = weight_low.*pdf_23_low + weight_high.*pdf_23_high;

% Normalise from p(beta1, x2,x3) to p(x2,x3 | beta1)
normC_sample_n = sum(sum(pdf_23,3),2);
pdf_23l1 = pdf_23./normC_sample_n;

% Draw samples from from p(x2,x3 | beta1)
beta_hat = zeros(3,Nsample);
for nn = 1:Nsample
    pdf_23l1_nn = squeeze(pdf_23l1(nn,:,:));
    beta23 = sample_from_pdf(X2,X3,pdf_23l1_nn,2,res);
    beta23 = beta23(:,1)';
    
    beta_hat(:,nn) = [beta_1_hat(nn), beta23];
end


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



