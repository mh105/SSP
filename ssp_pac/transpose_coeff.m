function A = transpose_coeff(B)
%% Transpose a block matrix B
    d = size(B,2);
    ARdeg = size(B,1)/d;
    
    A = zeros(d,d*ARdeg);
    for k=1:ARdeg
        A(:,(k-1)*d+1:k*d) = B((k-1)*d+1:k*d,:)';
    end

end