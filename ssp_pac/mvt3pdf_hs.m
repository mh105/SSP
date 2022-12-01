function [pdf,beta_0_grid,beta_1_grid,beta_2_grid,beta_0,beta_1,beta_2] =mvt3pdf_hs(dbeta,betaLim ,a, nu, V, b)
%% Estimate 3D Multivariate t distribution with parameter : 
%   a : mean ~ 3 x 1
%   nu : degrees of freedom
%   sigma = b*inv(V) covariance parameter (! different from covariance)
%   dbeta, betalim : parameter used to define a 3D square mesh
%
% OUTPUTS : pdf the estimate probability density function
%           beta_0_grid, beta_1_grid, beta_2_grid mesh at which pdf is estimated
%           beta_0, beta_1, beta_2 points at which pdf is estimated

% Build the mesh
beta_0 = linspace(-betaLim,betaLim,2*betaLim/dbeta);
beta_1 = linspace(-betaLim,betaLim,2*betaLim/dbeta);
beta_2 = linspace(-betaLim,betaLim,2*betaLim/dbeta);

[beta_0_grid,beta_1_grid,beta_2_grid] = ndgrid(beta_0,beta_1,beta_2);
beta_test_tot = [beta_0_grid(:)';beta_1_grid(:)';beta_2_grid(:)'];

% Estimate pdf
sigma = b*inv(V);
pdf = mvt3pdf_tmp(beta_test_tot,a,sigma,nu);

% Reshape pdf
pdf = reshape(pdf, length(beta_0),length(beta_1),length(beta_2));

end


function pdf = mvt3pdf_tmp(betas,a,sigma,nu)
%% Estimate 3D Multivariate t distribution at points betas with parameter : 
%   a : mean ~ 3 x 1
%   nu : degrees of freedom
%   sigma covariance parameter (! different from covariance)


if size(betas,1)~=3
    error('The current version only supports trivariate multivariate t')
else
    p=3;
end

gamma_quot_appro = @(vu) (1- 1./(4*vu) +1./(32*vu.^2)+5./(128*vu.^3));

% If nu is too big we approximate gamma(nu+1/2)/gamma(nu)
if nu<100
   k1 = gamma((nu+p)/2)/(gamma(nu/2) * (2*pi)^(p/2) * (nu/2)^(p/2)  * sqrt(det(sigma)) );
else
   k1 = (1+1/nu)*gamma_quot_appro(nu) / ( sqrt(det(sigma)) * (2*pi)^(p/2));
end

assert(size(betas,1) == size(a,1));

k2 = (1+ (1/nu)*(sum( (betas-a)'*inv(sigma).*(betas-a)', 2) )).^(-(p+nu)/2);
pdf = k1*k2';

end