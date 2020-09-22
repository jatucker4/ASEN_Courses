function outputvector = eulersMethod(tinitial,yinitial,step, tfinal)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
outputvector(1,1) = yinitial;
% last_time = tfinal+step;
time = tinitial:step:tfinal;
k = 1;
j = 2;
i = time(k);
while i ~= tfinal
    outputvector(1,j) = outputvector(1,j-1) + step*ode_project1(time(k));
    k = k+1;
    i = time(k);
    j = j+1;
end
end

