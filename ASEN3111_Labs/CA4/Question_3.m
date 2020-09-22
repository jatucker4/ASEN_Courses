function Question_3()
%Question_3 Performs all calculations and outputs for question 3 in CA4
%
% Author: Johnathan Tucker
%
% Collaborators: N/A
%
% This function has no inputs or direct outputs. However, it does create a
% plot of span efficiency factor versus taper ratio for different aspect
% ratios.
%
% Last Revised: 4/8/2020
%% First Create any necessary 
% Use some odd term number greater than 20
N = 50;

% Create a vector of taper ratios
taper_ratio_vec = linspace(0,1,N);

% Using thin airfoil theory with non varying lift slope
a0_r = 2*pi; % [rad]
a0_t = 2*pi; % [rad]

% Use a constant geometric aoa 
geo_r = 5; % [deg]
geo_t = 5; % [deg]

% Vector of Aspect Ratios
AR_vec = [4,6,8,10];

% Assume a zero lift aoa of 0 deg
aero_r = 0; % [deg]
aero_t = 0; % [deg]

% The span should be constant reuse the span from Q2
b = 100; % [ft]

%% Solve for the c_t and c_r values at the constant span for each AR
% First iterate through each aspect ratio
for i = 1:length(AR_vec)
    % Solve for the c_r values
    c_r(i,:) = (2*b)./(AR_vec(i).*(1+taper_ratio_vec)');
end
% Solve for the c_t values
c_t = c_r.*taper_ratio_vec;

%% Solve for the span efficiency factor values
% Iterate through each aspect ratio value
for i = 1:length(AR_vec)
    % Iterate through the taper ratio vector
    for j = 1:length(taper_ratio_vec)
        % Calculate the span efficiency factor
        [e,~,~] = PLLT(b,a0_t,a0_r,c_t(i,j),c_r(i,j),aero_t,aero_r,geo_t,geo_r,N);
        e_vec(i,j) = e;
    end
end

%% Plot the results
figure
plot(taper_ratio_vec,e_vec(1,:))
hold on
plot(taper_ratio_vec,e_vec(2,:))
hold on
plot(taper_ratio_vec,e_vec(3,:))
hold on
plot(taper_ratio_vec,e_vec(4,:))
title("Span Efficiency Factor vs. Taper Ratio",'FontSize',18)
xlabel("Taper Ratio, $\frac{c_t}{c_r}$",'FontSize',12)
ylabel("Span Efficiency Factor, e",'FontSize',12)
legend("AR = 4", "AR = 6", "AR = 8", "AR = 10")
