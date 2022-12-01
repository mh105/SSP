function [pdf,beta_0_grid,beta_1_grid,beta_2_grid,beta_0,beta_1,beta_2] = mvt3_const_pdf(dbeta,betamLim ,a, nu, V, b, W,f)
%% Estimate trivariate t distribution probability density function with multiple linear constraints.
% INPUTS :  dbeta the resolution at which we estimate pdf
%           a ~ 3 x 1 mean of the pdf
%           nu ~ degrees of freedom
%           b >=1
%           V ~ covariance parameter. Total covariance is b x inv(V) x (nu/ (nu - 2))
%           W ~ r x 3 and f ~ r x 1 defining linear constraints : W x beta > f
%
% OUTPUTS : pdf the estimate probability density function
%           beta_0_grid, beta_1_grid, beta_2_grid mesh at which pdf is estimated
%           beta_0, beta_1, beta_2 points at which pdf is estimated


[non_norm_mvt_pdf,beta_0_grid,beta_1_grid,beta_2_grid,beta_0,beta_1,beta_2] =mvt3pdf_hs(dbeta,betamLim ,a, nu, V, b);

beta_test_tot = [beta_0_grid(:)';beta_1_grid(:)';beta_2_grid(:)'];
non_norm_mvt_pdf = non_norm_mvt_pdf(:);

% Constraint
if length(f)>1
    troncate = all(W*beta_test_tot > f, 1);
elseif length(f)==1
    troncate =  sqrt(beta_test_tot(2,:).^2+beta_test_tot(3,:).^2) < 1*beta_test_tot(1,:);
end
non_norm_mvt_pdf = non_norm_mvt_pdf.*troncate';

% Normalise pdf
norm_mvt_pdf = non_norm_mvt_pdf./sum(non_norm_mvt_pdf);

% Reshape pdf
pdf = reshape(norm_mvt_pdf, length(beta_0),length(beta_1),length(beta_2));

end




