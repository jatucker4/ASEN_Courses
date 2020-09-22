function output = mean(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = 0;
for i = 1:length(input)
    output = output + input(i);
end
output = output/(length(input));
end

