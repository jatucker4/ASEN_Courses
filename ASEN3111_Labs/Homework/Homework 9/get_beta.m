function beta2 = get_beta(theta1,M_1)
%GET_BETA Summary of this function goes here
%   Detailed explanation goes here
syms beta1

eqn = tan(theta1) == 2*cot(beta1)*(((M_1^2)*(sin(beta1)^2) - 1)/((M_1^2) * (1.4 + cos(2*beta1)) + 2));

beta2 = double(solve(eqn,beta1));

if beta2 < 0
    beta2 = beta2 + pi;
end
end

