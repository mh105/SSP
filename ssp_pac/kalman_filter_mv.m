function [x_t_n,P_t_n,P_t_tmin1_n,logL,x_t_t,P_t_t,K_t,x_t_tmin1,P_t_tmin1] ...
    = kalman_filter_mv(y,M,Phi,Q,R,mu,Sigma,doBackward)
%% Estimate forward and backward multivariate Kalman Filter Estimate x_t_n = E(x_t|{y_t}t=1..n , P_t_s_n = cov(x_t, x_s | {y_t}))
% Model :
%   Observation  : y_t ~ q x T
%   Latent state : x_t ~ q x T
%   Augmented latent state : X_t = [x[t], ..., x[t-ARp]]' ~ (ARpxq)x T
%
%   Process Matrix : Phi ~ (ARpxq)x(ARpxq)
%   Process Noise  : Q   ~ (ARpxq)x(ARpxq)
%   Observation Matrix : M ~ qx(ARpxq)
%   Observation Noise  : R ~ q x q
%
%   y_t = M X_t + v_t, v_t ~ N(0,R)
%         =   x_t + v_t, v_t ~ N(0,R)
%   X_t = Phi X_t-1 + w_t, w_t ~ N(0,Q)

if nargin < 8
    doBackward = 1;
end

% 1-dimensional observation y
T = length(y);
p = length(mu);
q = size(R,1);

%assert(p==4);
assert(all(size(Phi)==[p,p]))
assert(all(size(Q)==[p,p]))
assert(all(size(R)==[q,q]))
assert(all(size(mu)==[p,1]))
assert(all(size(Sigma)==[p,p]))

J_t = zeros(p,p,T);
x_t_n = zeros(p,T+1);
P_t_n = zeros(p,p,T+1);
P_t_tmin1_n = zeros(p,p,T+1);



% Forward recursion (index 1 corresponds to t=0, etc)
x_t_tmin1 = zeros(p,T+1);
P_t_tmin1 = zeros(p,p,T+1);
K_t = zeros(p,q,T+1);
x_t_t = zeros(p,T+1);
P_t_t = zeros(p,p,T+1);

% initialize
x_t_t(:,1) = mu; % x_0_0
P_t_t(:,:,1) = Sigma; % P_0_0
logL = 0;


% Forward filtering
for ii=2:T+1
    
    %Predic
    x_t_tmin1(:,ii) = Phi*x_t_t(:,ii-1);
    P_t_tmin1(:,:,ii) = Phi*P_t_t(:,:,ii-1)*Phi' + Q;
    
    %Update
    K_t(:,:,ii) = P_t_tmin1(:,:,ii)*M'* inv(M*P_t_tmin1(:,:,ii)*M' + R);
    x_t_t(:,ii) = x_t_tmin1(:,ii) + K_t(:,:,ii)*(y(:,ii-1) - M*x_t_tmin1(:,ii));
    
    P_t_t(:,:,ii) = P_t_tmin1(:,:,ii) - K_t(:,:,ii)*M*P_t_tmin1(:,:,ii);
    
    % Innonvation form of the log likelihood
    logL = logL - 1/2 * (log(det(M*P_t_tmin1(:,:,ii)*M' + R)) + (y(:,ii-1)-M*x_t_tmin1(:,ii))'/(M*P_t_tmin1(:,:,ii)*M' + R)*(y(:,ii-1)-M*x_t_tmin1(:,ii)));
    
    if ~isreal(logL)
        return
    end
end

% Backward recursion


if doBackward
    x_t_n(:,end) = x_t_t(:,end); % x_n_n
    P_t_n(:,:,end) = P_t_t(:,:,end); % P_n_n
    
    for ii=T+1:-1:2
        J_t(:,:,ii-1) = P_t_t(:,:,ii-1)*Phi'/(P_t_tmin1(:,:,ii));
        x_t_n(:,ii-1) = x_t_t(:,ii-1) + J_t(:,:,ii-1)*(x_t_n(:,ii) - Phi*x_t_t(:,ii-1));
        P_t_n(:,:,ii-1) = P_t_t(:,:,ii-1) + J_t(:,:,ii-1)*(P_t_n(:,:,ii) - P_t_tmin1(:,:,ii))*J_t(:,:,ii-1)';
        
        P_t_tmin1_n(:,:,ii) = P_t_n(:,:,ii)*J_t(:,:,ii-1)';
    end
end
