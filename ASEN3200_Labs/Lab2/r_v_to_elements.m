function elements_vec = r_v_to_elements(r,v,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns an orbital elements vector of the form:
%   
%   [h,a,e,i,ohm,cap_ohm,theta]
%
% Created on: 9/17/2019
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output vector
elements_vec = zeros(1,7);

k = [0, 0, 1];

% Get h, e, n
h = cross(r,v);
e = (cross(v,h)/mu) - (r/norm(r));
n = cross(k,h);

% Get the norms for each vector
e_norm = norm(e);
h_norm = norm(h);
n_norm = norm(n);
r_norm = norm(r);
v_norm = norm(v);

% Get the semi-major axis
a = ((h_norm^2)/mu)/(1-e_norm^2);

% Get the inclination
i = acos(h(3)/h_norm);

% Find the angles and perform quadrent checks.
if(n(2) <0)
    cap_ohm = (2*pi) - acos(n(1)/n_norm);
else
    cap_ohm = acos(n(1)/n_norm);
end

if(e(3) <0)
    ohm = (2*pi) - acos(dot(n,e)/(n_norm*e_norm));
else
    ohm = acos(dot(n,e)/(n_norm*e_norm));
end

if(dot(r,v) < 0)
    theta = (2*pi) - acos(dot(e,r)/(r_norm*e_norm));
else
    theta = acos(dot(e,r)/(r_norm*e_norm));
end
% Save the outputs to their corresponding index in the output vector.
elements_vec(1) = h_norm;
elements_vec(2) = a;
elements_vec(3) = e_norm;
elements_vec(4) = i*(180/pi);
elements_vec(5) = ohm*(180/pi);
elements_vec(6) = cap_ohm*(180/pi);
elements_vec(7) = theta*(180/pi);
end

