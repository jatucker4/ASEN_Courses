%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to calculate the velocity of the collar using the
% geometric constants as well as the angular velocity of the disk and the
% angle theta
%
% Created by: Johnathan Tucker on 2/13/2019
%
% Inputs: 
%           radius = the distance between the center of the disk and the
%                    point A [m]
%           distance = the total distance from the slide to the center of
%                      disk [m]
%           length = length of the rod that connects points A and B [m]
%           theta = the angular position of the point A relative to the
%                   disk [radians]
%           w_disk = the angular velocity of the disk [rad/s]
%
% Outputs:
%           v_collar = the velocity of the collar [cm/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_collar] = LCSMODEL(radius,distance,length,theta,w_disk)
%% First calculate beta
beta = asin((distance - radius.*(sin(theta))) / length);
%% Calculate the velocity of the collar
v_collar = (-w_disk.*radius.*(sin(theta) + (cos(theta).*tan(beta))))*(100);

end

