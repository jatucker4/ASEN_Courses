clc;clear all;close all;
%% Begin code block for part b
A_eoA_t = 1.53;
p_0 = 1;
p_e_a = 0.94;

p_0op_e_a = p_0/p_e_a;
M_e_a = pop0_solveM(p_0op_e_a);

A_eoA_star_a = eq10_23_solve_area(M_e_a);
A_toA_star_a = (1/A_eoA_t)*A_eoA_star_a;

%% Begin code block for part b
p_e_b = 0.886;

p_0op_e_b = p_0/p_e_b;
M_e_b = pop0_solveM(p_0op_e_b);

A_eoA_star_b = eq10_23_solve_area(M_e_b);
A_toA_star_b = (1/A_eoA_t)*A_eoA_star_b;

%% Begin code block for part c
p_e_c = 0.75;

p_0op_e_c = p_0/p_e_c;
M_e_c = pop0_solveM(p_0op_e_c);

A_eoA_star_c = eq10_23_solve_area(M_e_c);
A_toA_star_c = (1/A_eoA_t)*A_eoA_star_c;

% Create guess vector
A2oAt = 1.2:0.01:1.3;
for i = 1:length(A2oAt)
    [Me_out(i),Pe(i)] = guess_iteration(A2oAt(i));
end
%% Begin code bloc for part d
p_e_d = 0.154;

p_0op_e_d = p_0/p_e_d;
M_e_d = pop0_solveM(p_0op_e_d);

A_eoA_star_d = eq10_23_solve_area(M_e_d);
A_toA_star_d = (1/A_eoA_t)*A_eoA_star_d;

%% Functions included
% <include> eq10_23_solve_area.m eq10_23_solve_M.m get_delta_s.m guess_iteration.m pop0_solveM.m</include>
