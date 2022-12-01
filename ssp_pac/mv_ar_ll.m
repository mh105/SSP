function mll = mv_ar_ll(y_t,r,C,c)
%% Compute -logL of a multidimentionnal Autoregressive process of order ARp with noisy observation y_t
%
% INPUTS : For k = 0 .. ARp C_k = cov(y[t], y[t-k])
%         C = toeplitz matrix : C_tmp_ij = C_[|i-j|] i,j = 0..ARp-1Inputs
%         c = [C_1' ... C_ARp']' ~ (ARpxq) x q
%         r : candidate observation noise covariance
%
% OUTPUTS : mll = -logL where L is the likelihood estimated with Kalman Filters

q = size(y_t,1);
ARp = size(C,1)/q;

% A ~ q x (qxARp)
% A_j = a_j ~ q x q j th AR coefficient
A_tmp = (C-r*eye(q*ARp))\c;
A = transpose_coeff(A_tmp);

% Observation noise
R = r*eye(q);

% Process noise : Q = C_0 - sum a_k(r)C_k
Q_tmp = C(1:q,1:q)-R;

for k=1:ARp
    Q_tmp = Q_tmp-A(:,(k-1)*q+1:k*q)*c((k-1)*q+1:k*q,:);
end

Q = zeros(q*ARp,q*ARp);
Q(1:q,1:q) = Q_tmp;
Q(q+1:q*ARp,q+1:q*ARp) = eps^(1/2)* eye(q*(ARp-1));

mu = zeros(q*ARp,1);
%mu = repmat([1;0.01;0], ARp,1);
Sigma = 1*eye(q*ARp,q*ARp);

% Get Transition Matrix
Phi = zeros(q*ARp,q*ARp);
Phi(q+1:q*ARp,1:q*(ARp-1)) = eye(q*(ARp-1));
Phi(1:q,:) = A;

% Get Observation Matrix
M = zeros(q,q*ARp);
M(1:q,1:q) = eye(q);

% Use Kalman Filter to estimate likelihood
doBackward =1;


[~,~,~,ll]= kalman_filter_mv(y_t,M,Phi,Q,R,mu,Sigma,doBackward);

if isreal(ll)~= 1
    Q_tmp = Q(1:3,1:3);
    Q_tmp = nearest_spd(Q_tmp);
    Q(1:3,1:3) = Q_tmp;
    [~,~,~,ll]= kalman_filter_mv(y_t,M,Phi,Q,R,mu,Sigma,doBackward);
    
    if isreal(ll)~= 1
      ll = -Inf;
    end
    warning(['Invalid Model. AR order p = ',num2str(ARp) ,' discarded'])
end

mll = -ll;



end

