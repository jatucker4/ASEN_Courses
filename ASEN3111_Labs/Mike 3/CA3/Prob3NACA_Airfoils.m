function [thetavals, dycdx] = Prob3NACA_Airfoils(m,p,t,c,N)
%Prob3NACA_Airfoils.m: computes the vectors of theta values and dyc/dx in
%terms of theta values required to implement Anderson eq. 4.61 for question
%3

%This function modifies the original NACA airfoils function to express the
%xvals vector in terms of theta and then express the derivative of the mean
%camber line function for the airfoil in question, dyc/dx, in terms of
%theta. This makes the integral in  eq. 4.61 more managable when computing
%the Thin airfoil theory estimate of the zero lift angle of attack for the
%airfoil being considered. 

% Author: Michael Martinson
% Collaborators: Isaac Goldner, Trace Valade
% Date: 3/28/2020

%%
%Assume same number of panels on top and bottom of airfoil: 
halfN = N/2; 

%define the difference in x between each panel's boundary points: 
dx = c/halfN; 

%Define x values for yt and yc calculations:
xvals = [0:dx:c]; %m 
 
%Create a vector for x/c for ease of other calculations: 
xoverc = xvals./c; %unitless                                                                                                                                                                                                                                                                                                                                                              

% %dyc/dx (Not used in new configuration)
% for i = 1:length(xvals)
%     
%     if xvals(i) <= p*c
%       dycdx(i) = ((2*m)/(p)) - ((2*m)/(p^2)).*(xoverc(i)); 
% 
%     else
%        dycdx(i) = ((-2*m)/((1-p)^2)).*(xoverc(i)) + ((2*p*m)/((1-p)^2)); 
%       
%     end 
%     
% end 

%% Problem 3 implementation: 
%convert x values to theta values: 
thetavals = acos(1 - 2.*xoverc)'; %rad

%Assign dyc/dx in terms of theta: 
for i = 1:length(thetavals)
    
    if thetavals(i) <= acos(1 - 2*p)
      dycdx(i,1) =((2*m)/(p)) -((2*m)/(p^2)).*(0.5*(1 - cos(thetavals(i)))); 

    else
       dycdx(i,1) =((-2*m)/((1-p)^2)).*(0.5*(1 - cos(thetavals(i)))) ...
           + ((2*p*m)/((1-p)^2)); 
      
    end 
end 

end
