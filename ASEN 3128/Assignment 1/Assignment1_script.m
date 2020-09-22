%% Housekeeping
clear all;
close all;
clc;
%% Create initial state vector and input to ode45
% This allows me to have a baseline trajectory to compare the sensitivity
% analysis to
y_init = [0 0 0 0 20 20];
wind = [0, 0, 0]; % [m/s]
mass = .03; % [kg]
tspan = [0 5]; % [s]

[t,y] = ode45(@(t,y) odefunc(t,y,wind,mass), tspan, y_init);

% Find the landing location: It'll be the index prior to the first negative
% index
land_loc = find(y(:,3) < 0, 1);
land_loc = land_loc - 1;

% Get baseline distances from the trajectory in each direction
N_distance = y(land_loc,1);
E_distance = y(land_loc,2);
U_distance = max(y(1:land_loc,3));

% Plot the baseline trajectory
figure(1)
plot3(y(land_loc,1),y(land_loc,2),y(land_loc,3),'.','Color','b')
hold on
plot3(y(1:land_loc,1),y(1:land_loc,2),y(1:land_loc,3))
title("Baseline Golfball Trajectory Plot",'Interpreter',...
    'latex','FontSize',26)
xlabel("North[m]",'Interpreter','latex','FontSize',26)
ylabel("East[m]",'Interpreter','latex','FontSize',26)
zlabel("Upward[m]",'Interpreter','latex','FontSize',26)

%% Perform sensitivity analysis on the wind vector
% Create a vector of velocity increments for the wind. ie: I'm going to
% analyze the balls sensitivity to winds between 1 and 5 m/s in increments
% of 0.5 m/s
vel_increments = 1:0.5:5;

% Loop through the wind velocity increments
for j = 1:length(vel_increments)
    
    % Create a normal distribution where negative wind just represents wind
    % in the negative N direction
    N = 100;
    N_wind = vel_increments(j) + 0.5*randn(1,N);
    
    % Loop through each random wind value
    for i = 1:N
        
        % Assign the wind the random value and run it through ode45
        wind = [N_wind(i), 0, 0];
        [t,y] = ode45(@(t,y) odefunc(t,y,wind,mass), tspan, y_init);

        % Find the landing location
        land_loc = find(y(:,3) < 0, 1);
        land_loc = land_loc - 1;
        
        % Only plot the first set of values to avoid unnecessary
        % computation time
        if i == 1
            figure(2)
            plot3(y(land_loc,1),y(land_loc,2),y(land_loc,3),'.','Color','b')
            hold on
            plot3(y(1:land_loc,1),y(1:land_loc,2),y(1:land_loc,3))
            title("Plot of Golfball Trajectories Subject to Wind Velocity Changes",...
            'Interpreter','latex','FontSize',26)
            xlabel("North[m]",'Interpreter','latex','FontSize',26)
            ylabel("East[m]",'Interpreter','latex','FontSize',26)
            zlabel("Upward[m]",'Interpreter','latex','FontSize',26)
        end
        
        % Create a vector of the landing locations in the North direction
        land_loc_vec(i,:) = y(land_loc,1);
    end
    
    % Use the above landing vector to calculate the average drift from the
    % baseline landing location
    avg_drift_vec(j,:) = mean(land_loc_vec-N_distance);
end
% Plot the effects of wind variation on the golfballs North direction drift
figure(3)
plot(vel_increments,avg_drift_vec)
xlabel("Wind Velocity [m/s]",'Interpreter','latex','FontSize',26)
ylabel("Average Crossrange Drift [m]",'Interpreter','latex','FontSize',26)
title("Average Crossrange Drift of a Golfball Subject to Various Wind Velocities",...
    'Interpreter','latex','FontSize',26)

%% Perform golfball mass sensitivity analysis
% Create a linearly spaced vector of mass increments for the golf ball. 
% ie: I'm going to analyze the balls sensitivity to mass changes between 
% 0.03 and 0.1 kg.
mass_increments = linspace(0.03,0.1,10);

% Reset the wind vector so that there's only one independent variable at a
% time
wind = [0,0,0];

% Loop through the mass increments
for j = 1:length(mass_increments)
    
    % Create a uniform distribution of masses with a standard deviation of
    % 10 grams. Uniform distribution is used so that the mass cannot be
    % negative
    mass_change = mass_increments(j) + 0.01*rand(1,N);
    
    % Loop through the number of mass possibilities
    for i=1:N
        
        % Get the random mass and input it to ode45
        mass = mass_change(i);
        [t,y] = ode45(@(t,y) odefunc(t,y,wind,mass), tspan, y_init);

        % Find the landing location
        land_loc = find(y(:,3) < 0, 1);
        land_loc = land_loc - 1;
        
        % Only plot the first set of values to avoid unnecessary
        % computation time
        if i == 1
            figure(4)
            plot3(y(land_loc,1),y(land_loc,2),y(land_loc,3),'.','Color','b')
            hold on
            plot3(y(1:land_loc,1),y(1:land_loc,2),y(1:land_loc,3))
            title("Plot of Golfball Trajectories Subject to Mass Changes",...
            'Interpreter','latex','FontSize',26)
            xlabel("North[m]",'Interpreter','latex','FontSize',26)
            ylabel("East[m]",'Interpreter','latex','FontSize',26)
            zlabel("Upward[m]",'Interpreter','latex','FontSize',26)
        end
        
        % Create a vector of the landing locations in the East direction
        land_loc_vec(i,:) = y(land_loc,2);
        
    end
    % Use the above landing vector to calculate the average drift from the
    % baseline landing location
    avg_drift_vec(j,:) = mean(land_loc_vec-E_distance);
    
end
% Plot the effects of mass variation on the golfballs East direction drift
figure(5)
plot(mass_increments,avg_drift_vec)
xlabel("Golf Ball Mass [kg]",'Interpreter','latex','FontSize',26)
ylabel("Average Downrange Drift [m]",'Interpreter','latex','FontSize',26)
title("Average Downrange Drift of a Golfball Subject to Various Masses ",...
    'Interpreter','latex','FontSize',26)

