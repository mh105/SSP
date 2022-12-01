function [a,w]=get_rotmat_pam(Phi)
% Phi = a x R(w) where R is a rotation matrix 
assert(all(size(Phi)==[2,2]))
a= sqrt(Phi(1,1)^2+Phi(2,1)^2);
w=atan(Phi(2,1)/Phi(1,1));
w=atan2(Phi(2,1),Phi(1,1));
% if w<0
%     w=w+pi;
% end
end