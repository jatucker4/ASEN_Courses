function M_num = pm_eqn_m(v)
%PM_EGN_M Summary of this function goes here
%   Detailed explanation goes here
syms M
PM = v == sqrt(2.4/.4)*atan(sqrt((.4/2.4)*(M^2 - 1))) - atan(sqrt(M^2 - 1));
M_num = double(solve(PM,M));
end

