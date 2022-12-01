function [a_post,nu_post, V_post,b_post] = get_post_param(y_t,X_t,a_prior,nu_prior,V_prior,b_prior)
%% Get posterior parameters for a (constrained) linear model. Based on W.W.Davis (1978).

V_post = V_prior + X_t'*X_t;
a_post = V_post\(V_prior*a_prior  +  X_t'*y_t);
nu_post= nu_prior + length(y_t);

beta_ols_hat = (X_t'*X_t)\X_t'*y_t;

Q_post = (y_t - X_t*beta_ols_hat)'*(y_t - X_t*beta_ols_hat) ...
    + (a_prior-a_post)'*V_prior*(a_prior-a_post) ...
    + (beta_ols_hat - a_post)' *((X_t)'*(X_t))*(beta_ols_hat - a_post);

b_post = (nu_prior*b_prior+ Q_post)/nu_post;

end