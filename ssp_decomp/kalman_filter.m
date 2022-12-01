function [x_t_n,P_t_n,P_t_tmin1_n,logL,x_t_t,P_t_t,K_t,x_t_tmin1,P_t_tmin1] = kalman_filter(y,Phi,Q,R,mu,Sigma)
%% Forward and backward Kalman Filter Estimate 
% INPUTS : y,Phi,Q,R,mu,Sigma such that:
%          y_t   ~   M x_t   + N(0,R)
%          x_t+1 ~ Phi x_t + N(0,Q)
%
% OUTPUTS: 
%         x_t_n = E(x_t|{y_t}t=1..n) , 
%         P_t_s_n = cov(x_t, x_s | {y_t}))
%         logL: log likelihood
%         K_t: optimal Kalman Gain


% 1-dimensional observation y
T = length(y);
p = length(mu);
q = 1;

Mo=[1 0];
Comp=floor(p/2)+1;
M=repmat(Mo, 1,Comp);
M=M(1,1:p);

%assert(p==4);
assert(all(size(Phi)==[p,p]))
assert(all(size(Q)==[p,p]))
assert(all(size(R)==[q,q]))
assert(all(size(mu)==[p,1]))
assert(all(size(Sigma)==[p,p]))

% Forward recursion (index 1 corresponds to t=0, etc)
x_t_tmin1 = zeros(p,T+1);
P_t_tmin1 = zeros(p,p,T+1);
K_t = zeros(p,T+1);
x_t_t = zeros(p,T+1);
P_t_t = zeros(p,p,T+1);

% initialize
x_t_t(:,1) = mu; % x_0_0
P_t_t(:,:,1) = Sigma; % P_0_0
logL = 0;


for ii=2:T+1
    
    %Predic
    x_t_tmin1(:,ii) = Phi*x_t_t(:,ii-1);
    P_t_tmin1(:,:,ii) = Phi*P_t_t(:,:,ii-1)*Phi' + Q;
 
    %Update
    K_t(:,ii) = P_t_tmin1(:,:,ii)*M'/(M*P_t_tmin1(:,:,ii)*M' + R);
    x_t_t(:,ii) = x_t_tmin1(:,ii) + K_t(:,ii)*(y(ii-1) - M*x_t_tmin1(:,ii));
    P_t_t(:,:,ii) = P_t_tmin1(:,:,ii) - K_t(:,ii)*M*P_t_tmin1(:,:,ii);
    
    % Innonvation form of the log likelihood
    logL = logL - 1/2 * (log(det(M*P_t_tmin1(:,:,ii)*M' + R)) + (y(ii-1)-M*x_t_tmin1(:,ii))'/(M*P_t_tmin1(:,:,ii)*M' + R)*(y(ii-1)-M*x_t_tmin1(:,ii)));

    if ~isreal(logL)
        return
    end
end

% Backward recursion
J_t = zeros(p,p,T);
x_t_n = zeros(p,T+1);
P_t_n = zeros(p,p,T+1);
P_t_tmin1_n = zeros(p,p,T+1);

x_t_n(:,end) = x_t_t(:,end); % x_n_n
P_t_n(:,:,end) = P_t_t(:,:,end); % P_n_n

for ii=T+1:-1:2
    J_t(:,:,ii-1) = P_t_t(:,:,ii-1)*Phi'/(P_t_tmin1(:,:,ii));
    x_t_n(:,ii-1) = x_t_t(:,ii-1) + J_t(:,:,ii-1)*(x_t_n(:,ii) - Phi*x_t_t(:,ii-1));
    P_t_n(:,:,ii-1) = P_t_t(:,:,ii-1) + J_t(:,:,ii-1)*(P_t_n(:,:,ii) - P_t_tmin1(:,:,ii))*J_t(:,:,ii-1)';
    
    P_t_tmin1_n(:,:,ii) = P_t_n(:,:,ii)*J_t(:,:,ii-1)';
end
