function beta = mvt3rnd_hs(a,sigma,nu, cases)
%% Draw -cases- samples from 3D Multivariate t distribution with parameter :
%   a : mean ~ 3 x 1
%   nu : degrees of freedom
%   sigma covariance parameter (! different from covariance)

p = size(a,1);
mu = zeros(p,1);
y = mvnrnd(mu,sigma,cases)';
u = chi2rnd(nu,1,cases);

beta_centered = y./sqrt(u/nu);
beta = beta_centered + a;
end