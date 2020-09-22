function [x,y] = NACA_Airfoil(m,p,t,c,N)
%NACA_Airfoil Performs the calculations necessary to get the x and y
%vectors that describe the specified NACA airfoil
%
% Author: Johnathan Tucker
%
% Collaborators: N/A

% This function takes in the max chord value "m", the location of max chord
% "p", the thickness "t", the chord length "c", and the number of panels to
% use "N". This function outputs the x and y vectors that describe the 
% specified NACA airfoil.
%
% Last Revised: 3/26/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create vector for percentage of chord length
x = linspace(0,c,N/2);
% Create a vector of half thicknesses
y_t = (t/0.2)*c*(0.2969.*sqrt(x./c) - 0.1260.*(x./c) - 0.3516.*(x./c).^2 +...
    0.2843.*(x./c).^3 - 0.1036.*(x./c).^4);

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
%Adding a check for NaN values
y_c(isnan(y_c)) = 0;

% Create zeta to solve for the upper and lower x/y values
zeta = atan2(diff(y_c),diff(x));
zeta = [zeta,0];
% Solve for x_u and x_l
x_u = x - y_t.*sin(zeta);
x_l = x + y_t.*sin(zeta);

% Solve for y_u and y_l
y_u = y_c + y_t.*cos(zeta);
y_l = y_c - y_t.*cos(zeta);

% Combine upper and lower vectors to get final x and y
x = [flip(x_u),x_l(2:end)];
y = -[flip(y_u),y_l(2:end)];
end

