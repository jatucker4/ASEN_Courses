function v_out = system_eqs(t,v)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
v_out = zeros(2,1);
v_out(1) = -8.*v(1).*(v(1)-0.15).*(v(1)-1)-v(1).*v(2);
v_out(2) = (0.002 + ((0.2.*v(2))./(v(1) + 0.3))).*(-v(2) -8.*v(1).*(v(1)-0.15-1));
end

