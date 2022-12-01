function plot_ci_phase_hs(phase_ub, phase_lb, ta,col)
%% Plots credible/confidence interval for a variable encoding phase
% Arguments : 
%   - phase_ub,phase_lb upper and lower bound of the variable
%   - ta, time vector 
%   - col, color used for the shaded interval 

if size(ta,1)>size(ta,2)
    ta=ta';
end
assert(size(ta,1)==1);

dta = diff(ta)/2;

ta1=[ ta(1), ta(2:end)-dta] ;
ta2=[ ta(1:end-1)+dta,ta(end)] ;

inToPi = find(phase_ub>phase_lb);
noToPi = find(phase_ub<=phase_lb);

for ii = inToPi   
    heighi = phase_ub(1,ii) - phase_lb(1,ii);
    widthi = ta2(1,ii)-ta1(1,ii);
    rectangle('Position', [ta1(1,ii) phase_lb(1,ii) widthi heighi], 'facecolor', col, 'linestyle', 'none' )
end

for ii = noToPi   
    heighi_1 =  pi - phase_lb(1,ii);
    heighi_2 =  phase_ub(1,ii) - (-pi);
    widthi = ta2(1,ii)-ta1(1,ii);
  
    rectangle('Position', [ta1(1,ii) phase_lb(1,ii) widthi heighi_1], 'facecolor', col, 'linestyle', 'none' )
    rectangle('Position', [ta1(1,ii) -pi widthi heighi_2], 'facecolor', col, 'linestyle', 'none' )  
end


end