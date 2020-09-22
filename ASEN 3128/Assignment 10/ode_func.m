function dydt = ode_func(t,y,A)
% This function is to be passed into ode45 to propogate lateral dynamics
%% Create output vector
dydt = A*y;
end

