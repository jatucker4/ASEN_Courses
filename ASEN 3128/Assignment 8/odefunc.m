function dydt = odefunc(t,y,full_mat)
%odefunc is the function meant to be used by ode45 for assignment 7 in ASEN
%3128

% The output vector is the product of A +B*K and the state vector
dydt = full_mat*y;
end

