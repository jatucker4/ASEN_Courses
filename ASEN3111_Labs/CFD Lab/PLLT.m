function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%PLLT Performs Prandtl lifting line theory calculations that take
% advantage of the Fourier sine series to obtain circulation. This 
% ultimately leads to the calculation and output of span efficiency factor,
% coefficient of lift, and coefficient of induced drag.
%
% Author: Johnathan Tucker
%
% Collaborators: N/A

% Inputs:
%           b: Span [ft]
%           a0_t: cross-sectional lift slope at the wingtips [-]
%           a0_r: cross-sectional lift slope at the wing roots [-]
%           c_t: chord at the tips [ft]
%           c_r: chord at the wing root [ft]
%           aero_t: zero-lift angle of attack at tips [deg]
%           aero_r: zero-lift angle of attack at wing root [de]
%           geo_t: geometric angle of attack at tips [deg]
%           geo_r: geometric angle of attack at wing root [deg]
%           N: number of odd terms in series
%
% Outputs:
%           e: span efficiency factor [-]
%           c_L: coefficient of lift [-]
%           c_Di: coefficient of induced drag [-]
% Last Revised: 3/26/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert all variables to compatible units
aero_t = aero_t*(pi/180);
aero_r = aero_r*(pi/180);
geo_t = geo_t*(pi/180);
geo_r = geo_r*(pi/180);
%% First create the theta vector
i = 1:1:N;
theta_vec = i.*(pi./(2*N));
%% Now create anonymous functions or vectors for wing properties
% Wing Taper
c = @(theta) c_r - (c_r - c_t).*cos(theta);
% Variable cross-sectional lift slope
a_0 = @(theta) a0_r - (a0_r - a0_t).*cos(theta);
% Aerodynamic Twist
alpha_lzero = aero_r - (aero_r - aero_t).*cos(theta_vec);
% Geometric Twist
alpha = geo_r - (geo_r - geo_t).*cos(theta_vec);

%% Now Create the matrix to solve for the Fourier Coefficients
% Preallocate to save time
A = zeros(N,N);
% Prior to this I'll create a vector that is N in length of odd values
odds = 1:2:2*N;
% First we need to iterate through every theta value
for j = 1:N
    % Then we'll need to iterate through every coefficient multiplier
    % corresponding to this theta value
    for k = 1:N
        A(j,k) = (4*b*sin(odds(k)*theta_vec(j)))/...
            (a_0(theta_vec(j))*c(theta_vec(j))) + odds(k)*sin(odds(k)*...
            theta_vec(j))/sin(theta_vec(j));
    end
end
%% With the A matrix we can now solve for the for the Fourier Coefficients
% First get the "b" vector using the geometric and aerodynamic twist
b_vec = alpha - alpha_lzero;
coeffs = A\b_vec';

%% Now get the required outputs
% Get the aspect ratio
AR = (b^2)/((c_r + c_t)*(b/2));

% Get the c_L
c_L = coeffs(1)*pi*AR;

% Get the delta vector
delta = odds(2:end).*(coeffs(2:end)'./coeffs(1)').^2;
delta = sum(delta);
% Use delta to solve for span efficiency
e = 1/(1+delta);

% Finally calculate c_Di
c_Di = (c_L^2)/(pi*e*AR);

end


