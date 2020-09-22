function v_out = system_eqs_1(t,v, T)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
v_out = [0;0];
if(mod(t,T) >= 10 && mod(t,T) <=13)
    S = 0.25;
else
    S = 0;
end
v_out(1) = -8*v(1)*(v(1)-0.15)*(v(1)-1)-v(1)*v(2) + S;
v_out(2) = (0.002 + ((0.2*v(2))/(v(1) + 0.3)))*(-v(2) -8*v(1)*(v(1)-0.15-1));
end