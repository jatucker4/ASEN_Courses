function v = pm_eqn_v(M)
%PM_EQN_V Summary of this function goes here
%   Detailed explanation goes here
v = sqrt(2.4/.4)*atan(sqrt((.4/2.4)*(M^2 - 1))) - atan(sqrt(M^2 - 1));
end

