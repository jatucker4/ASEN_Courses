function M_out = eq10_23_solve_M(area_ratio,large)
%EQ10_23_SOLVE_M Summary of this function goes here
%   Detailed explanation goes here
syms M

equation = area_ratio^2 == (1/M^2)*((2/2.4)*(1 + (.4/2)*M^2))^(2.4/.4);

M1 = double(solve(equation,M));
if large == 1
    M_out = real(M1(5));
elseif large == 0
    M_out = real(M1(6));
end
end

