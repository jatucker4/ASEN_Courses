function ang_accel = getAccel(omega,theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(omega)-1
    ang_accel(i) = (omega(i+1) - omega(i)) / (theta(i+1) - theta(i));
end

end

