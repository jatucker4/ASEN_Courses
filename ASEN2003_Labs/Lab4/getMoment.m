function moment = getMoment(omega,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct a monte carlo analysis to find the moment value that gives the
% minimum distance between the model 2 omega values and the the
% experimentally achieved omega values
%
% Created by: Johnathan Tucker
%
% Inputs:
%           omega: Experimentally achieved angular velocity values (rad/s)
%           theta: Experimentally acheived angular position values (rad)
%
% Outputs:
%           moment: The moment value that minimizes the distance between
%           the model and experimental values (N-m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seed the randomizer for repeatability
rng(1)
% Get  a random set of moment values
n = 10000;
moment_test = rand(1,n);
% Get the outputs from the model 2 function
monte_outputs = model_2(moment_test,theta(20));
% Get the index of the min distance between the actual and modeled output
distance = sqrt((monte_outputs - omega(20)).^2);
moment_index = find(distance==min(distance));
moment = moment_test(1,moment_index);
end

