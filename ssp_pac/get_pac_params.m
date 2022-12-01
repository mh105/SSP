function [K, Phi, A_0] = get_pac_params(regressed_pac)
%% GET_PAC_PARAMS deduce K_mod and Phi_mod from the linear regression vector
% INPUTS   : -Latent pac ~ 3 x Nsamp    x Nwind
%             Nsamp number of drawn samples and Nwind number of time windows
%
% OUTPUTS : K, Phi ~ Nsamp x Nwind, modulation strength and amplitude
%           such that: for t = 1..n, i = 1..Nwind :
%           A_fast_t = A_0(1 + K_i cos(-Phi_i + Phi_slow_t))

dim_pac = size(regressed_pac,1);
Nsamp   = size(regressed_pac,2);
Nwind   = size(regressed_pac,3);

if (nargout == 3) && (dim_pac == 2)
    error('A_0 is estimated for 3D linear regression only (size(latent_pac,1) = 3)')
end

if dim_pac == 2
    K = reshape(sqrt(regressed_pac(1,:,:).^2+regressed_pac(2,:,:).^2), Nsamp, Nwind);
    Phi = reshape(atan2(regressed_pac(2,:,:),regressed_pac(1,:,:)), Nsamp, Nwind);
    
elseif dim_pac==3
    K = reshape(sqrt(regressed_pac(2,:,:).^2+regressed_pac(3,:,:).^2)./regressed_pac(1,:,:), Nsamp, Nwind);
    Phi = reshape(atan2(regressed_pac(3,:,:),regressed_pac(2,:,:)), Nsamp, Nwind);
    
    if (nargout == 3)
        A_0 = reshape(regressed_pac(1,:,:), Nsamp, Nwind);
    end
    
else
    error('latent_pac first dimmension should be 2 or 3')
end

end