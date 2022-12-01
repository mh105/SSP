function [phase_mean,phase_sup,phase_inf] = get_phase_msi(Phi_k_t)
% INPUTS  : Phase series Phi_k_t ~ Nsample x Ntimepts
% OUTPUTS : Mean phase, Inf phae and sup phase ~ 1 x Ntimepts

phase_mean = angle(mean(exp(1i*Phi_k_t),1));
Ntimepts = size(phase_mean,2);

% Find the maximum phase distance
phase_dist = sin((Phi_k_t-phase_mean)/2).^2;
[~,sortedPhi_id] = sort(phase_dist,1);
Phi_sup_tmp = zeros(1, Ntimepts);

for tt = 1:Ntimepts
    Phi_sup_tmp(1,tt) = Phi_k_t(sortedPhi_id(end,tt),tt);
end
dPhi = abs(wrapToPi(phase_mean-Phi_sup_tmp));

% Define upper and lower phase
phase_sup = wrapToPi(phase_mean+dPhi);
phase_inf = wrapToPi(phase_mean-dPhi);

end