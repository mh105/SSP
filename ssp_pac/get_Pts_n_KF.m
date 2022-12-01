function [x_t_n,Pstruct] = get_Pts_n_KF(y,Phi,Q,R,mu,Sigma, Smax)

% Returns the smoothed estimates x_t_n = E(x(t) | y(1) .... y(T))
% And the smoothed Covariances in a structure :
%           P_tmink_t_n = E (x(t-k), x(t) | y(1) .... y(T))
%           P_tmink_t_n ~ (pxp) x (T+1)x(k_max+1);
%           P_tmink_t_n  (:,:, i,j)  =  P_[(i-1)-(j-1),(i-1)]_n; i=1..T+1 / j=1..kmax+1

% y is the observation
% Phi, Q, R, mu and Sigma are the model parameters
% k_max is the maximum temporal correlation


%% Usual forward-backward Kalman Filtering -> x_t_n, P_t_n
T = length(y);
p = length(mu);
q = 1;
M=  repmat([1,0], 1,floor(p/2));

if nargin<7
    Smax=inf;
end

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
C_t = zeros(T+1,p);
x_t_t = zeros(p,T+1);
P_t_t = zeros(p,p,T+1);

% initialize
x_t_t(:,1) = mu; % x_0_0
P_t_t(:,:,1) = Sigma; % P_0_0
logL = 0;


for ii=2:T+1
    %Prediction
    x_t_tmin1(:,ii) = Phi*x_t_t(:,ii-1);
    P_t_tmin1(:,:,ii) = Phi*P_t_t(:,:,ii-1)*Phi' + Q;
    x_t_tmin1(:,ii) = Phi*x_t_t(:,ii-1);
    P_t_tmin1(:,:,ii) = Phi*P_t_t(:,:,ii-1)*Phi' + Q;
    
    %Update
    K_t(:,ii) = P_t_tmin1(:,:,ii)*M'/(M*P_t_tmin1(:,:,ii)*M' + R);
    x_t_t(:,ii) = x_t_tmin1(:,ii) + K_t(:,ii)*(y(ii-1) - M*x_t_tmin1(:,ii));
    P_t_t(:,:,ii) = P_t_tmin1(:,:,ii) - K_t(:,ii)*M*P_t_tmin1(:,:,ii);
    % Innonvation form of the log likelihood
    logL = logL - 1/2 * (log(abs(det(M*P_t_tmin1(:,:,ii)*M' + R))) + ...
        (y(ii-1)-M*x_t_tmin1(:,ii))'/(M*P_t_tmin1(:,:,ii)*M' + R)*(y(ii-1)-M*x_t_tmin1(:,ii)));
    
    if ~isreal(logL)
        return
    end
end


% First Backward recursion
J_t = zeros(p,p,T);
x_t_n = zeros(p,T+1);
P_t_n = zeros(p,p,T+1);
x_t_n(:,end) = x_t_t(:,end); % x_n_n
P_t_n(:,:,end) = P_t_t(:,:,end); % P_n_n

for ii=T+1:-1:2
    J_t(:,:,ii-1) = P_t_t(:,:,ii-1)*Phi'/(P_t_tmin1(:,:,ii));
    x_t_n(:,ii-1) = x_t_t(:,ii-1) + J_t(:,:,ii-1)*(x_t_n(:,ii) - Phi*x_t_t(:,ii-1));
    P_t_n(:,:,ii-1) = P_t_t(:,:,ii-1) + J_t(:,:,ii-1)*(P_t_n(:,:,ii) - P_t_tmin1(:,:,ii))*J_t(:,:,ii-1)';    
end


%% Get all the covariances2
% we have J_t ~ (pxp)x T
% P_s_t^n= J_s+1 P_s+1,t^ n 
PP=zeros(p*(T+1),p*(T+1));

J_tf=zeros(p,p,T+1);
J_tf(:,:,2:T+1)=J_t;
J_tf(:,:,1)=P_t_t(:,:,1)*Phi'/(P_t_tmin1(:,:,2));

for tt=1:T+1
   % Case ss=tt 
   PP(p*(tt-1)+1 :p*tt, p*(tt-1)+1 :p*tt )= P_t_n(:,:,tt);
   Jcurr=p*(tt-1)+1 :p*tt;
   
   Smax_curr=max(tt-Smax,1);
   
   for ss=tt-1:-1:Smax_curr
       Iprev=p*(tt-1)+1 -(tt-ss-1)*p :p*tt -(tt-ss-1)*p;
       Icurr=p*(tt-1)+1 -(tt-ss)*p :p*tt -(tt-ss)*p;
       
       Psp1_t_n= PP(Iprev ,Jcurr);
       Psp_t_n=J_tf(:,:,ss+1)*Psp1_t_n;
       
       
       PP(Icurr ,Jcurr)=Psp_t_n;
       PP(Jcurr,Icurr) =Psp_t_n';
   end
end




Pstruct=struct();
Pstruct.P_t_n=P_t_n;
Pstruct.PP=PP;
end