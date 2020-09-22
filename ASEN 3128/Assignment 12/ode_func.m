function dydt = ode_func(t,y,A,B,K)

% Create a check for the closed loop time constants
check_mat = A + B*K;
check_evals = eig(check_mat);
omega_d_dutch_check = imag(check_evals(3,1));
nat_freq_dutch_check  = sqrt(omega_d_dutch_check^2 + real(check_evals(3,1))^2);
damping_ratio_dutch_check  = -real(check_evals(3,1))/nat_freq_dutch_check;

check_DR_time_const = -1/real(check_evals(3,:));
check_spiral_time_const = -1/check_evals(5);

if check_DR_time_const > 40
    fprintf("Dutch Roll time constant failed check and is: %f\n",check_DR_time_const)
end

if damping_ratio_dutch_check < 0.35
    fprintf("Dutch Roll damping ratio failed check and is: %f\n",damping_ratio_dutch_check)
end
if check_spiral_time_const > 20
    fprintf("Spiral time constant failed check and is: %f\n",check_spiral_time_const)
end

% Compute the dydt vector
dydt = (A + B*K)*y;


end

