function [pac, pbins, mi_kl, mi_l1] = get_pac_from_kphi(K,Phi,Npb)
%% Computes parametric phase amplitude modulogram and modulation indices
% INPUTS  : modulation phase Phi and Amplitude K ~ 1 x Ntimepts
% OUTPUTS : phase amplitude coupling pac ~ Npb x Ntimepts
%           Modulation Index Defined as :
%           KL divergence with uniform distribution : mi_kl
%           L1 distance to uniform distributiom     : mi_l1
%

if size(K,1) > size(K,2)
   K = K'; 
end

if size(Phi,1) > size(Phi,2)
   Phi = Phi'; 
end

assert(size(K,1) == 1);
assert(size(Phi,1) == 1);
assert(size(K,2)==size(K,2));

pbins = linspace(-pi,pi, Npb);
dPhi = (2*pi)/Npb;

pac  = (1+  sin(dPhi/2)/(dPhi/2)*K.*cos(pbins' - Phi))/(2*pi);

% Make sure that for all t = 1 .. Ntimepts : 
% sum_t (pac(phi, t) dphi) = 1. 
renorm_pac = @(PAC) PAC./sum(PAC)/dPhi;
pac = renorm_pac(pac);

if nargout>2
    mi_kl = sum(pac.*log2(pac *2*pi),1)*dPhi;
end

if nargout>3
     mi_l1 = sum(abs(pac -1/(2*pi)),1)*dPhi;
end


end