function [Phi_hat,M,Q_hat,R_hat,mll] = mvyule_ar(y_t,ARp)
%% Fit an Autoregressive (AR) process of order ARp to a multimmensionnal latent series x_t through noisy observation y_t using Yule's equations
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
%
% OUTPUTS : - Optimal parameter estimated through numerical optimization (adapted from Matsuda and al.2017)
%             Phi_hat,M,Q_hat,R_hat
%           - mll = -logL where L is the likelihood

q = size(y_t,1);
T = size(y_t,2);

% C_k = cov(y[t], y[t-k]) k = 0 .. ARp
C_k = zeros(q,q,ARp+1);
%y_t = y_t-mean(y_t,2)*ones(1,T);
for k=1:ARp+1
    C_k(:,:,k) = y_t(:,1:T-k+1)*y_t(:,k:T)'/T;
end

% C_tmp = toeplitz matrix : C_tmp_ij = C_[|i-j|] i,j = 0..ARp-1
C_tmp = zeros(q*ARp,q*ARp);
for k=1:ARp
    for j=1:k
        C_tmp((k-1)*q+1:k*q,(j-1)*q+1:j*q) = C_k(:,:,k-j+1);
    end
    for j=k+1:ARp
        C_tmp((k-1)*q+1:k*q,(j-1)*q+1:j*q) = C_k(:,:,j-k+1)';
    end
end

% c = [C_1' ... C_ARp']'
c = zeros(q*ARp,q);
for k=1:ARp
    c((k-1)*q+1:k*q,:) = C_k(:,:,k+1);
end

% C = full toeplitz matrix : C_tmp_ij = C[|i-j|] i,j = 0..ARp
C = zeros(q*(ARp+1),q*(ARp+1));
C(1:q*ARp,1:q*ARp) = C_tmp;
for k=1:ARp
    C((k-1)*q+1:k*q,ARp*q+1:(ARp+1)*q) = c((ARp-k)*q+1:(ARp-k+1)*q,:);
    C(ARp*q+1:(ARp+1)*q,(k-1)*q+1:k*q) = c((ARp-k)*q+1:(ARp-k+1)*q,:)';
end
C(ARp*q+1:(ARp+1)*q,ARp*q+1:(ARp+1)*q) = C_k(:,:,1);

% Get C's minimum eigen value -Rmin- and minimize -logL(R) for 0<R<Rmin
eigs = eig(C);
Rmin = min(eigs);
if Rmin <= 0
    r_hat = 0;
    mll = mv_ar_ll(y_t,r_hat,C_tmp,c);
else
    [r_hat,mll] = fminbnd(@(r)mv_ar_ll(y_t,r,C_tmp,c),0,Rmin);
end

% Get observation covariance matrix
R_hat = r_hat*eye(q);

% A_hat_tmp ~ (qxARp) x q
% A_hat_tmp_i = a_i ~ q x q i th AR coefficient
A_hat_tmp = (C_tmp-r_hat*eye(q*ARp))\c;

% A ~ q x (qxARp)
% A_j = a_j ~ q x q j th AR coefficient
A_hat = transpose_coeff(A_hat_tmp);

% Get process noise covariance matrix
Q_hat_tmp = C_k(:,:,1)-R_hat;
for k=1:ARp
    Q_hat_tmp = Q_hat_tmp-A_hat(:,(k-1)*q+1:k*q)*C_k(:,:,k+1);
end

Q_hat = zeros(q*ARp,q*ARp);
Q_hat(1:q,1:q) = Q_hat_tmp;
Q_hat(q+1:q*ARp,q+1:q*ARp) = eps^(1/2)* eye(q*(ARp-1));

% Get Transition Matrix
Phi_hat = zeros(q*ARp,q*ARp);
Phi_hat(q+1:q*ARp,1:q*(ARp-1)) = eye(q*(ARp-1));
Phi_hat(1:q,:) = A_hat;

% Get Observation Matrix
M = zeros(q,q*ARp);
M(1:q,1:q) = eye(q);

end

