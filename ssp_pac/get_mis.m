function [X_mean, X_inf, X_sup] = get_mis(X_mis)

if (size(X_mis,1)~=3 && size(X_mis,2)~=3 )
    error('Invalid format')
elseif size(X_mis,1)~=3 
    X_mis=X_mis';
end

X_mean = X_mis(1,:);
X_inf  = X_mis(2,:);
X_sup  = X_mis(3,:);
    
end

