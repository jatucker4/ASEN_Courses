function aoa_zero_lift = get_zlaoa_tat(m,p,c,N)
%GET_ZLAOA_TAT Calculates the zero lift angle of attack using thin airfoil
%theory
%
% Author: Johnathan Tucker
% Collaborators: Michael Martinson
% This function takes in the max camber "m", location of max camper "p",
% the chord length "c", and the number of panels "N". It outputs the zero
% lift angle of attack calculated using thin airfoil theory.
%
% Last Revised: 3/27/2020

% Create vector for percentage of chord length
x = linspace(0,c,N);
% I need to find the index of x that is closest to the p*c value
[~,index] = min(abs(x-p*c));

% Now loop through the x values using a conditional statement that
% replicates the peicewise function
for i = 1:length(x)
    if i <= index
        y_c(i) = m.*(x(i)./p^2).*(2*p - x(i)./c);
    else
        y_c(i) = m.*((c-x(i))./(1-p)^2).*(1 + x(i)./c - 2*p);
    end
end

% Adding a check for NaN values
y_c(isnan(y_c)) = 0;

% Differentiate with respect to x
d_y_dx = diff(y_c)./diff(x);

% Create a theta vector from 0 to pi
theta = linspace(0,pi,length(d_y_dx));

% Calculate the zero lift AOA for outputting
aoa_zero_lift = (-1/pi)*trapz(d_y_dx.*(cos(theta)-1));
end

