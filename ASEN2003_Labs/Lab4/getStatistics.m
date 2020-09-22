function statistics = getStatistics(risiduals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the statistics for the risidual omega values
% between the experimental and modeled data
%
% Created by: Johnathan Tucker
%
% Inputs:
%           risiduals: Difference between experimental and modeled angular
%                      velocities (rad/s)
%
% Outputs:
%           statistics: A matrix of the standard deviation (rad/s), mean of
%                       the risiduals (rad), number of trials, standard 
%                       error of the mean (rad/s), and number of outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standard_dev = std(risiduals);
risidual_mean = mean(risiduals);
N = length(risiduals);
sem = standard_dev./sqrt(N);
sigma_3 = length(find(risiduals>3*standard_dev));
if isempty(sigma_3)
    sigma_3 = 0;
end
statistics = [standard_dev,risidual_mean,sem,N,sigma_3];
end

