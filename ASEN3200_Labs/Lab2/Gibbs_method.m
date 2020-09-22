function [r_out,v_out] = Gibbs_method(r_1,r_2,r_3,mu)
% This function performs Gibbs method. Taking in 3 position vectors and the
% gravitational parameter of the focus of the orbit and outputs the
% velocity that correspons to the second position vector
%
% Date created: 9/17/2019
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the norms for each position vector
r_1_norm = norm(r_1);
r_2_norm = norm(r_2);
r_3_norm = norm(r_3);

% Get the cross products
C_1_2 = cross(r_1,r_2);
C_2_3 = cross(r_2,r_3);
C_3_1 = cross(r_3,r_1);

% Get the unit vector for the second position vector
u_r = r_2/r_2_norm;
% Create a flag to know if the method is applicable
flag_vec = C_2_3/norm(C_2_3);
flag = dot(u_r,flag_vec);
% Create a tolernace out to floating point double precision
tolerance = 1*(10)^(-16);

% If the flag is less than the tolerance calculate N,D,S
if  flag < tolerance
    N = r_1_norm*C_2_3 + r_2_norm*C_3_1 +...
        r_3_norm*C_1_2;
    
    D = cross(r_1,r_2) + cross(r_2,r_3) + cross(r_3,r_1);
    
    S = r_1*(r_2_norm - r_3_norm) + r_2*(r_3_norm - r_1_norm) + ...
        r_3*(r_1_norm - r_2_norm);
    

end
    % Use N,D,S to calcuate the velocity and output both the velocity and
    % position
    v_out = sqrt(mu/(norm(D)*norm(N)))*((cross(D,r_2)/r_2_norm)+S);
    r_out = r_2;
end

