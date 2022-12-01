function [Phi_new, Q_new,R_new, mu_new]=m_step(y,NNharm, x_t_n, P_t_n,P_t_tmin1_n,w_fond_range)
%% Maximization step for the extended State-Space oscillator decomposition (Matsuda and Komaki 2017).
% INPUTS:
%   y                 : signal for a given window
%   NNharm            : 1 x #indpt oscillations an array containing the number of harmonics
%                       for each independant oscillation
%   w_fond_range      : range to look for fundamental frequency
%                     in case of harmonic estimation
%   x_t_n, P_t_n, ... : state means and covariance estimates from E-step
%
% OUTPUTS:
%   Updated parameters

Nosc =  sum(NNharm);       % Total number of oscillation
Nosc_ind = length(NNharm); % Number of independant oscillation

T=length(y);
Phi_new=zeros(2*Nosc);
Q_new  =zeros(2*Nosc);

% See shumway and al. 1986
A = sum(P_t_n(:,:,1:end-1),3) + x_t_n(:,1:end-1)*x_t_n(:,1:end-1)';
B = sum(P_t_tmin1_n(:,:,2:end),3) + x_t_n(:,2:end)*x_t_n(:,1:end-1)';
C = sum(P_t_n(:,:,2:end),3) + x_t_n(:,2:end)*x_t_n(:,2:end)';

% Update initial means mu
mu_new = x_t_n(:,1);

% Define observation matrix M and update observation noise R
Mo=[1 0]; M=repmat(Mo, 1,Nosc);
R_new = 1/T * ((y'-M*x_t_n(:,2:end))*(y'-M*x_t_n(:,2:end))' + M*squeeze(sum(P_t_n(:,:,2:end),3))*M');


for n_osc_ind=1:Nosc_ind
    
    % Update Phi and Q for an oscillation with no harmonics
    if NNharm(1,n_osc_ind)==1
        n_harm=1;
        curId= 2*( sum(NNharm(1,1:n_osc_ind)) - NNharm(1,n_osc_ind) + n_harm) -1;
        
        A_tmp = A (curId:curId+1,curId:curId+1);
        B_tmp = B (curId:curId+1,curId:curId+1);
        C_tmp = C (curId:curId+1,curId:curId+1);
        
        a_tmp     = max(sqrt(rt(B_tmp)^2+trace(B_tmp)^2)/abs(trace(A_tmp)),0.01);
        %w_tmp     = atan(rt(B_tmp)/trace(B_tmp));  %if w_tmp<0; w_tmp=w_tmp+pi; end
        w_tmp     = atan2(rt(B_tmp),trace(B_tmp));
        
        sigma2_tmp= (1/(2*T))*(trace(C_tmp)- a_tmp^2*trace(A_tmp));
        
        Phi_tmp=get_rot_mat(a_tmp,w_tmp);
        
        Phi_new(curId:curId+1,curId:curId+1)= Phi_tmp;
        Q_new  (curId:curId+1,curId:curId+1)= eye(2)* sigma2_tmp;
        
        
        % Update Phi and Q for an oscillation containing harmonics
    else
        
        %update_a = 1; % Update the harmonic scaling parameter a_h
        Nharm_cur= NNharm(1, n_osc_ind);% Number of harmonics
        
        rHarm=1:1:Nharm_cur;
        rtBr=zeros(1,Nharm_cur);
        trBr=zeros(1,Nharm_cur);
        trAr=zeros(1,Nharm_cur);
        trCr=zeros(1,Nharm_cur);
        
        for n_harm=1:Nharm_cur
            curId= 2*( sum(NNharm(1,1:n_osc_ind)) - NNharm(1,n_osc_ind) + n_harm) -1;
            A_tmp=A(curId:curId+1,curId:curId+1);
            B_tmp=B(curId:curId+1,curId:curId+1);
            C_tmp=C(curId:curId+1,curId:curId+1);
            
            trAr(1,n_harm)=trace(A_tmp);
            trBr(1,n_harm)=trace(B_tmp);
            rtBr(1,n_harm)=   rt(B_tmp);
            trCr(1,n_harm)=trace(C_tmp);
        end
        
        
        %eif update_a==1
            
            % Define expected log likelihood G = E (logL | y1...yn)
            Gharm= @(w, sigma2r, ar) -T* sum(log(sigma2r),2) -(1/2)* sum(   sigma2r.^(-1).*(trCr-ar.^2.*trAr)   ,2);
            dGdw = @(w) sum( rHarm.* (  rtBr.*trBr.*cos(2.*rHarm.*w) + 0.5*(rtBr.^2-trBr.^2).*sin(2.*rHarm.*w)   ) ...
                ./  (  (trAr.*trCr-0.5*(rtBr.^2+trBr.^2))  -  ( 0.5*(trBr.^2-rtBr.^2).*cos(2.*rHarm.*w)+rtBr.*trBr.*sin(2.*rHarm.*w) )    ) );
            
            % set the interval to look for optimal w_h
            w_r_min=w_fond_range(1,n_osc_ind);
            w_r_max=w_fond_range(2,n_osc_ind);
            
            % Find all the zeros of dG/dw_heart
            w_r_opt_k=AllMin(dGdw,w_r_min,w_r_max,500);
           
            
            if isempty(w_r_opt_k)
                w_r_opt=0.5*(w_r_min+w_r_max);
                disp('no min')
            else
                Gmax=-Inf;
                % In case of multiple local minima, find the one which minimizes G=log(L) itself
                Gtot_tmp = zeros(1,length(w_r_opt_k));
                for k=1:length(w_r_opt_k)
                    % if k==2;disp('mult');end
                    w_tmp=w_r_opt_k(1,k);
                    a_r_tmp      = max((  trBr.*cos(rHarm.*w_tmp)+rtBr.*sin(rHarm*w_tmp)  ) ./trAr, 0.01);
                    sigma2_r_tmp = (1/(2*T)) .*(trCr-a_r_tmp.^2.*trAr);
                    Gtmp=Gharm(w_tmp , sigma2_r_tmp, a_r_tmp);
                    Gtot_tmp(1,k)=Gtmp;
                    if Gtmp>Gmax && isreal(Gtmp)
                        Gmax=Gtmp;
                        w_r_opt=w_tmp;
                    end
                end
                
            end
            
            
            a_r_new      = max((  trBr.*cos(rHarm.*w_r_opt)+rtBr.*sin(rHarm*w_r_opt)  ) ./trAr,0.01);
            sigma2_r_new = (1/(2*T)) .*(trCr-a_r_new .^2.*trAr);
            
            
            % Let us try to not update a
%         else
%             a_h_0 = 0.9999999; 
%             RT_h = zeros(1,Nharm_cur);
%             dd_h = zeros(1,Nharm_cur);
%             theta_h = zeros(1,Nharm_cur);
%             
%             for n_harm=1:Nharm_cur           
%                 dd_h(1,n_harm) = trCr(1,n_harm) + a_h_0^2 * trAr(1,n_harm);
%                 RT_h(1,n_harm) = (trBr(1,n_harm)^2 + rtBr(1,n_harm)^2)^(1/2);
%                 theta_h(1,n_harm) = atan2(rtBr(1,n_harm), trBr(1,n_harm));
%             end
%             
%             Gharm= @(w, sigma2r, ar) -T* sum(log(sigma2r),2) -(1/2)* sum(   sigma2r.^(-1).*(trCr-ar.^2.*trAr)   ,2);
%             dGdw = @(w)sum( -rHarm .* RT_h .* sin(rHarm * w -theta_h) ./ (dd_h - 2 * a_h_0 * RT_h .* cos(rHarm * w -theta_h))  );
% 
%             % set the interval to look for optimal w_heart
%             w_r_min=w_fond_range(1,n_osc_ind);
%             w_r_max=w_fond_range(2,n_osc_ind);
%             
%             % Find all the zeros of dG/dw_heart
%             w_r_opt_k=AllMin(dGdw,w_r_min,w_r_max,500);
%             
%             w_test = linspace(w_r_min,w_r_max,500);
%             dg0 = zeros(1,length(w_test));
%             dg1 = zeros(1,length(w_test));
%             dg2 = zeros(1,length(w_test));
%             
%             for ww= 1 :length(w_test)
%                 dg0(1,ww) = dGdw(w_test(ww));
%                 dg1(1,ww) =  sum( rHarm .* RT_h .* cos(rHarm * w_test(ww) -theta_h));
%                 dg2(1,ww) = sum((dd_h - 2 * a_h_0 * RT_h .* sin(rHarm * w_test(ww) -theta_h)));
%                
%             end
%             
%             if isempty(w_r_opt_k)
%                 w_r_opt=0.5*(w_r_min+w_r_max);
%                 disp('no min')
%                 a_r_new      = a_h_0 * ones(1,Nharm_cur);
%                 sigma2_r_new = (1/(2*T)) .*(trCr + a_r_new.^2.*trAr - 2 * a_r_new.* RT_h .* cos(rHarm * w_r_opt -theta_h));
%                 
%             else
%                 Gmax=-Inf;
%                 % In case of multiple local minima, find the one which minimizes G=log(L) itself
%                 Gtot_tmp = zeros(1,length(w_r_opt_k));
%                 for k=1:length(w_r_opt_k)
%                     % if k==2;disp('mult');end
%                     w_tmp=w_r_opt_k(1,k);
%                     a_r_tmp      = a_h_0 * ones(1,Nharm_cur);
%                     sigma2_r_tmp = (1/(2*T)) .*(trCr + a_r_tmp.^2.*trAr - 2 * a_r_tmp.* RT_h .* cos(rHarm * w_tmp -theta_h));
%                     
%                     Gtmp=Gharm(w_tmp , sigma2_r_tmp, a_r_tmp);
%                     Gtot_tmp(1,k)=Gtmp;
%                     if Gtmp>Gmax && isreal(Gtmp)
%                         Gmax=Gtmp;
%                         w_r_opt=w_tmp;
%                         sigma2_r_new = sigma2_r_tmp;
%                         a_r_new = a_r_tmp;
%                     end
%                 end        
%                 
%             end
%             
%             
%         end
        
        % Build optimal harmonic functions
        for n_harm=1:Nharm_cur
            
            curId= 2*( sum(NNharm(1,1:n_osc_ind)) - NNharm(1,n_osc_ind) + n_harm) -1;
            Phi_tmp=get_rot_mat(a_r_new(1,n_harm),n_harm*w_r_opt);
            
            Phi_new(curId:curId+1, curId:curId+1) = Phi_tmp;
            Q_new  (curId:curId+1, curId:curId+1) = sigma2_r_new(1,n_harm);
            Q_new  (curId:curId+1, curId:curId+1) = sigma2_r_new(1,n_harm);
        end
        
    end
    
end

end

function z=AllMin(f,xmin,xmax,N)
% Find all zero crossing from positive to negative
if (nargin<4)
    N=100;
end
dx=(xmax-xmin)/N;
x2=xmin;
y2=f(x2);
z=[];
for i=1:N
    x1=x2;
    y1=y2;
    x2=xmin+i*dx;
    y2=f(x2);
    if (y1*y2<=0) && (y1>0)
        options = optimset('Display','off');
        z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2),options)];
    end
end
end



