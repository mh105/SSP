function Phir=get_rot_mat(a,w)
% Phi = a x R(w) where R is a rotation matrix 
Phir=zeros(2,2);
Phir(1,1)= a *cos(w);
Phir(1,2)=-a *sin(w);
Phir(2,1)= a *sin(w);
Phir(2,2)= a *cos(w);
end
