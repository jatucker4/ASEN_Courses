function [x,y] = NACA_Airfoils(m,p,t,c,N)
%NACA_Airfoils.m: This function computes the x and y coordinates for the
%boundary points required to implement the vortex panel method for a given
%airfoil geometry. 

%This function executes the procedure outlined in the last page of the CA3
%document to compute the x and y coordinates for the boundary points for a
%given NACA airfoil. The airfoil's m,p,t, and c values must be entered
%along with N, the number of panels used. 

% Author: Michael Martinson
% Collaborators: Isaac Goldner, Trace Valade
% Date: 3/27/2020

%% Define x values for yt and yc calculations: 

%Assume same number of panels on top and bottom of airfoil: 
halfN = N/2; 

%define the difference in x between each panel's boundary points: 
dx = c/halfN; 

%Define x values for yt and yc calculations:
xvals = [0:dx:c]; %m 
 
%Create a vector for x/c for ease of other calculations: 
xoverc = xvals./c; %unitless                                                                                                                                                                                                                                                                                                                                                              

%Compute the yt and yc vectors as defined in the CA3 handout: 
%yt
yt = ((t.*c)./(0.2)).*...
    (0.2969.*sqrt(xoverc) - 0.1260.*(xoverc) - 0.3516.*((xoverc).^2) ...
    + 0.2843.*((xoverc).^3) - 0.1036.*((xoverc).^4));                      

%yc
for i = 1:length(xvals)
    
    if xvals(i) <= p*c
       yc(i) = m.*(xvals(i)/(p^2)).*(2.*p - xoverc(i)); %m

    else
        yc(i) = m.*((c - xvals(i))./((1 - p).^2)).*...
            (1 + (xoverc(i)) - (2.*p)); %m  
    end 
    
end 

%Compute dyc/dx as defined in the CA3 handout: 
for i = 1:length(xvals)
    
    if xvals(i) <= p*c
      dycdx(i) = ((2*m)/(p)) - ((2*m)/(p^2)).*(xoverc(i)); 

    else
       dycdx(i) = ((-2*m)/((1-p)^2)).*(xoverc(i)) + ((2*p*m)/((1-p)^2)); 
    end 
end 

 %Assign xi vector: 
xi = atan(dycdx); %rad

%% Assign upper and lower surface coordinates: 

xu = xvals - yt.*sin(xi); 
xL = xvals + yt.*sin(xi); 

yu = yc + yt.*cos(xi); 
yL = yc - yt.*cos(xi); 

%Output the final x and y values, have the first and last points be the
%trailing edge: 
    x = [flip(xL(2:end))';xu']; %m
    y = [flip(yL(2:end))';yu']; %m
end
