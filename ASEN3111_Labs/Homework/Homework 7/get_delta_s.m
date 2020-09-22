function  delta_s = get_delta_s(M)
%GET_DELTA_S Summary of this function goes here
%   Detailed explanation goes here
cp = 1.4*287/(.4);
delta_s = cp*log( (1+ (2.8/2.4)*(M^2 - 1))*((2 + .4*M^2)/(2.4*M^2))) - ...
    287*log( (1+ (2.8/2.4)*(M^2 - 1)));
end

