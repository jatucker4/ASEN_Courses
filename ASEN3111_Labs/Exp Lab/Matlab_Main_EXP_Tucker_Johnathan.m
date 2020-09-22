%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 - Experimental Lab 1
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to directly create the plots for each of the
% first experimental lab bullet points
%
% Created Date: 3/7/2020
%
% Change Log: 
%           - 3/7/2020: Code up the plots for the first bullet 
%           - 3/10/2020: Code up the plots for the second bullet 
%           - 3/11/2020: Code up the plots for the third bullet 
%           - 3/11/2020: Code up the plots for the fourth bullet 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%% Begin code block for the first bullet of the analysis section
% First process the airfoil data
load('Parsed_Airfoil_Data_EXP_Tucker_Johnathan.mat');
% Iteratively calculate the wake velocity and deficit and plot them
figure
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1)
for i = 1:length([WTData.speed])
    if WTData(i).bad_data ~= 2
        if WTData(i).speed == 15
            vel_wake_vec_15_airfoil{i} = sqrt(2.*WTData(i).Pwake_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_tunnel_vec_15_airfoil{i} = sqrt(2.*WTData(i).Pwindtunnel_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_deficit_vec_15_airfoil{i} =  vel_tunnel_vec_15_airfoil{i} - vel_wake_vec_15_airfoil{i};
            y_loc_vec_15_airfoil{i} = WTData(i).y_mm_mean;
            rho_vec_15_airfoil{i} = WTData(i).rhoatm_kgm3_mean;
            x_loc_vec_15_airfoil{i} = WTData(i).x_location .* ones(1,length(WTData(i).y_mm_mean));
            quiver(x_loc_vec_15_airfoil{i},y_loc_vec_15_airfoil{i},vel_deficit_vec_15_airfoil{i},zeros(1,length(vel_deficit_vec_15_airfoil{i})),.5)
            hold on    
        end
    end
end
xlabel("$X-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Y-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Airfoil\:Velocity\:Deficit$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
for i = 1:length([WTData.speed])
    if WTData(i).bad_data ~= 2
        if WTData(i).speed == 25
            vel_wake_vec_25_airfoil{i-14} = sqrt(2.*WTData(i).Pwake_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_tunnel_vec_25_airfoil{i-14} = sqrt(2.*WTData(i).Pwindtunnel_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_deficit_vec_25_airfoil{i-14} =  vel_tunnel_vec_25_airfoil{i-14} - vel_wake_vec_25_airfoil{i-14};
            y_loc_vec_25_airfoil{i-14} = WTData(i).y_mm_mean;
            rho_vec_25_airfoil{i-14} = WTData(i).rhoatm_kgm3_mean;
            x_loc_vec_25_airfoil{i-14} = WTData(i).x_location .* ones(1,length(WTData(i).y_mm_mean));
            quiver(x_loc_vec_25_airfoil{i-14},y_loc_vec_25_airfoil{i-14},vel_deficit_vec_25_airfoil{i-14},zeros(1,length(vel_deficit_vec_25_airfoil{i-14})),.5)
            hold on    
        end
    end
end
xlabel("$X\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Y\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Airfoil\:Velocity\:Deficit$",'Interpreter','latex','FontSize',18)
sgtitle("$Airfoil\:Velocity\:Deficit\:vs\:Y-Location\:Plots$",'Interpreter','latex','FontSize',18)
% Repeat for the cylinder
load('Parsed_Cylinder_Data_EXP_Tucker_Johnathan.mat');
figure
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1)
for i = 1:length([WTData.speed])
    if WTData(i).bad_data ~= 2
        if WTData(i).speed == 15
            vel_wake_vec_15_cylinder{i} = sqrt(2.*WTData(i).Pwake_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_tunnel_vec_15_cylinder{i} = sqrt(2.*WTData(i).Pwindtunnel_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_deficit_vec_15_cylinder{i} =  vel_tunnel_vec_15_cylinder{i} - vel_wake_vec_15_cylinder{i};
            y_loc_vec_15_cylinder{i} = WTData(i).y_mm_mean;
            rho_vec_15_cylinder{i} = WTData(i).rhoatm_kgm3_mean;
            x_loc_vec_15_cylinder{i} = WTData(i).x_location .* ones(1,length(WTData(i).y_mm_mean));
            quiver(x_loc_vec_15_cylinder{i},y_loc_vec_15_cylinder{i},vel_deficit_vec_15_cylinder{i},zeros(1,length(vel_deficit_vec_15_cylinder{i})))
            hold on    
        end
    end
end
xlabel("$X\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Y\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Cylinder\:Velocity\:Deficit$",'Interpreter','latex','FontSize',18)
sgtitle("$Cylinder\:Velocity\:Deficit\:vs\:Y-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
for i = 1:length([WTData.speed])
    if WTData(i).bad_data ~= 2
        if WTData(i).speed == 25
            vel_wake_vec_25_cylinder{i-15} = sqrt(2.*WTData(i).Pwake_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_tunnel_vec_25_cylinder{i-15} = sqrt(2.*WTData(i).Pwindtunnel_Pa_mean./WTData(i).rhoatm_kgm3_mean);
            vel_deficit_vec_25_cylinder{i-15} =  vel_tunnel_vec_25_cylinder{i-15} - vel_wake_vec_25_cylinder{i-15};
            y_loc_vec_25_cylinder{i-15} = WTData(i).y_mm_mean;
            rho_vec_25_cylinder{i-15} = WTData(i).rhoatm_kgm3_mean;
            x_loc_vec_25_cylinder{i-15} = WTData(i).x_location .* ones(1,length(WTData(i).y_mm_mean));
            quiver(x_loc_vec_25_cylinder{i-15},y_loc_vec_25_cylinder{i-15},vel_deficit_vec_25_cylinder{i-15},zeros(1,length(vel_deficit_vec_25_cylinder{i-15})))
            hold on    
        end
    end
end
xlabel("$X\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Y\:Location\:[mm]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Cylinder\:Velocity\:Deficit$",'Interpreter','latex','FontSize',18)
%% Begin code for the second bullet point
% Iterate and get the half width at each x location for the airfoil
for i = 1:length(vel_tunnel_vec_15_airfoil)
    if ~isempty(x_loc_vec_15_airfoil{i})
        x_vec_15_airfoil(i) = max(x_loc_vec_15_airfoil{i}); 
        vel_deficit_vec_15_airfoil_max(i) = max(vel_deficit_vec_15_airfoil{i});
        half_width_vec_15_airfoil(i) = get_half_width(vel_deficit_vec_15_airfoil{i},y_loc_vec_15_airfoil{i});
    end
end
figure
set(gcf, 'Position', get(0, 'Screensize'));
x_vec_15_airfoil(find(x_vec_15_airfoil == 0)) = [];
vel_deficit_vec_15_airfoil_max(find(vel_deficit_vec_15_airfoil_max == 0)) = [];
half_width_vec_15_airfoil(find(half_width_vec_15_airfoil == 0)) = [];

% average every two indices
x_vec_15_increment = [1,3,6,8];
for i = 1:length(x_vec_15_airfoil)
    if max(ismember(x_vec_15_increment,i)) == 1
        half_width_mean_15_airfoil(i) = mean([half_width_vec_15_airfoil(i),half_width_vec_15_airfoil(i+1)]);
        airfoil_15_vmax_mean(i) = mean([vel_deficit_vec_15_airfoil_max(i),vel_deficit_vec_15_airfoil_max(i+1)]);
    elseif i == 5 || i == 10
        half_width_mean_15_airfoil(i) = half_width_vec_15_airfoil(i);
        airfoil_15_vmax_mean(i) = vel_deficit_vec_15_airfoil_max(i);
    end  
end
half_width_mean_15_airfoil(find(half_width_mean_15_airfoil == 0)) = [];
airfoil_15_vmax_mean(find(airfoil_15_vmax_mean == 0)) = [];
airfoil_15_halfwidth_fit = polyfit(unique(x_vec_15_airfoil),half_width_mean_15_airfoil,1);
airfoil_15_vmax_fit = polyfit(unique(x_vec_15_airfoil),airfoil_15_vmax_mean,2);

subplot(1,2,1)
scatter(x_vec_15_airfoil,half_width_vec_15_airfoil)
hold on
plot(unique(x_vec_15_airfoil), polyval(airfoil_15_halfwidth_fit,unique(x_vec_15_airfoil)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Half\:Width\:[mm]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Airfoil\:Half\:Width\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
sgtitle("$Airfoil\:15\:m/s\:Half\:Width\:and\:Max\:Velocity\:Deficit\:vs\:X-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
scatter(x_vec_15_airfoil,vel_deficit_vec_15_airfoil_max)
hold on 
plot(unique(x_vec_15_airfoil), polyval(airfoil_15_vmax_fit,unique(x_vec_15_airfoil)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Max\:Velocity\:Deficit\:[m/s]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Airfoil\:Max\:Deficit\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(vel_tunnel_vec_25_airfoil)
    if ~isempty(x_loc_vec_25_airfoil{i})
        x_vec_25_airfoil(i) = max(x_loc_vec_25_airfoil{i}); 
        vel_deficit_vec_25_airfoil_max(i) = max(vel_deficit_vec_25_airfoil{i});
        half_width_vec_25_airfoil(i) = get_half_width(vel_deficit_vec_25_airfoil{i},y_loc_vec_25_airfoil{i});
    end
end
figure
set(gcf, 'Position', get(0, 'Screensize'));
x_vec_25_airfoil(find(x_vec_25_airfoil == 0)) = [];
vel_deficit_vec_25_airfoil_max(find(vel_deficit_vec_25_airfoil_max == 0)) = [];
half_width_vec_25_airfoil(find(half_width_vec_25_airfoil == 0)) = [];

% average every two indices
x_vec_25_increment = [1,3,6,8];
for i = 1:length(x_vec_25_airfoil)
    if max(ismember(x_vec_25_increment,i)) == 1
        half_width_mean_25_airfoil(i) = mean([half_width_vec_25_airfoil(i),half_width_vec_25_airfoil(i+1)]);
        airfoil_25_vmax_mean(i) = mean([vel_deficit_vec_25_airfoil_max(i),vel_deficit_vec_25_airfoil_max(i+1)]);
    elseif i == 5 || i == 10
        half_width_mean_25_airfoil(i) = half_width_vec_25_airfoil(i);
        airfoil_25_vmax_mean(i) = vel_deficit_vec_25_airfoil_max(i);
    end  
end
half_width_mean_25_airfoil(find(half_width_mean_25_airfoil == 0)) = [];
airfoil_25_vmax_mean(find(airfoil_25_vmax_mean == 0)) = [];
airfoil_25_halfwidth_fit = polyfit(unique(x_vec_25_airfoil),half_width_mean_25_airfoil,1);
airfoil_25_vmax_fit = polyfit(unique(x_vec_25_airfoil),airfoil_25_vmax_mean,2);

subplot(1,2,1)
scatter(x_vec_25_airfoil,half_width_vec_25_airfoil)
hold on
plot(unique(x_vec_25_airfoil), polyval(airfoil_25_halfwidth_fit,unique(x_vec_25_airfoil)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Half\:Width\:[mm]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Airfoil\:Half\:Width\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
sgtitle("$Airfoil\:25\:m/s\:Half\:Width\:and\:Max\:Velocity\:Deficit\:vs\:X-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
scatter(x_vec_25_airfoil,vel_deficit_vec_25_airfoil_max)
hold on 
plot(unique(x_vec_25_airfoil), polyval(airfoil_25_vmax_fit,unique(x_vec_25_airfoil)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Max\:Velocity\:Deficit\:[m/s]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Airfoil\:Max\:Deficit\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)

% Start the cylinder section
% Iterate and get the half width at each x location for the airfoil
for i = 1:length(vel_tunnel_vec_15_cylinder)
    if ~isempty(x_loc_vec_15_cylinder{i})
        x_vec_15_cylinder(i) = max(x_loc_vec_15_cylinder{i}); 
        vel_deficit_vec_15_cylinder_max(i) = max(vel_deficit_vec_15_cylinder{i});
        half_width_vec_15_cylinder(i) = get_half_width(vel_deficit_vec_15_cylinder{i},y_loc_vec_15_cylinder{i});
    end
end
figure
set(gcf, 'Position', get(0, 'Screensize'));
x_vec_15_cylinder(find(x_vec_15_cylinder == 0)) = [];
vel_deficit_vec_15_cylinder_max(find(vel_deficit_vec_15_cylinder_max == 0)) = [];
half_width_vec_15_cylinder(find(half_width_vec_15_cylinder == 0)) = [];

% average every two indices
x_vec_15_increment = [1,3,6];
for i = 1:length(x_vec_15_cylinder)
    if max(ismember(x_vec_15_increment,i)) == 1
        half_width_mean_15_cylinder(i) = mean([half_width_vec_15_cylinder(i),half_width_vec_15_cylinder(i+1)]);
        cylinder_15_vmax_mean(i) = mean([vel_deficit_vec_15_cylinder_max(i),vel_deficit_vec_15_cylinder_max(i+1)]);
    elseif i == 5 || i == 8
        half_width_mean_15_cylinder(i) = half_width_vec_15_cylinder(i);
        cylinder_15_vmax_mean(i) = vel_deficit_vec_15_cylinder_max(i);
    end  
end
half_width_mean_15_cylinder(find(half_width_mean_15_cylinder == 0)) = [];
cylinder_15_vmax_mean(find(cylinder_15_vmax_mean == 0)) = [];
cylinder_15_halfwidth_fit = polyfit(unique(x_vec_15_cylinder),sort(half_width_mean_15_cylinder),1);
cylinder_15_vmax_fit = polyfit(unique(x_vec_15_cylinder),sort(cylinder_15_vmax_mean),2);

subplot(1,2,1)
scatter(x_vec_15_cylinder,half_width_vec_15_cylinder)
hold on
plot(unique(x_vec_15_cylinder), polyval(cylinder_15_halfwidth_fit,unique(x_vec_15_cylinder)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Half\:Width\:[mm]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Cylinder\:Half\:Width\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
sgtitle("$Cylinder\:15\:m/s\:Half\:Width\:and\:Max\:Velocity\:Deficit\:vs\:X-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
scatter(x_vec_15_cylinder,vel_deficit_vec_15_cylinder_max)
hold on 
plot(unique(x_vec_15_cylinder), flip(polyval(cylinder_15_vmax_fit,unique(x_vec_15_cylinder))))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Max\:Velocity\:Deficit\:[m/s]$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Cylinder\:Max\:Deficit\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(vel_tunnel_vec_25_cylinder)
    if ~isempty(x_loc_vec_25_cylinder{i})
        x_vec_25_cylinder(i) = max(x_loc_vec_25_cylinder{i}); 
        vel_deficit_vec_25_cylinder_max(i) = max(vel_deficit_vec_25_cylinder{i});
        half_width_vec_25_cylinder(i) = get_half_width(vel_deficit_vec_25_cylinder{i},y_loc_vec_25_cylinder{i});
    end
end
figure
set(gcf, 'Position', get(0, 'Screensize'));
x_vec_25_cylinder(find(x_vec_25_cylinder == 0)) = [];
vel_deficit_vec_25_cylinder_max(find(vel_deficit_vec_25_cylinder_max == 0)) = [];
half_width_vec_25_cylinder(find(half_width_vec_25_cylinder == 0)) = [];

% average every two indices
x_vec_25_increment = [1,3,6];
for i = 1:length(x_vec_25_cylinder)
    if max(ismember(x_vec_25_increment,i)) == 1
        half_width_mean_25_cylinder(i) = mean([half_width_vec_25_cylinder(i),half_width_vec_25_cylinder(i+1)]);
        cylinder_25_vmax_mean(i) = mean([vel_deficit_vec_25_cylinder_max(i),vel_deficit_vec_25_cylinder_max(i+1)]);
    elseif i == 5 || i == 8
        half_width_mean_25_cylinder(i) = half_width_vec_25_cylinder(i);
        cylinder_25_vmax_mean(i) = vel_deficit_vec_25_cylinder_max(i);
    end  
end
half_width_mean_25_cylinder(find(half_width_mean_25_cylinder == 0)) = [];
cylinder_25_vmax_mean(find(cylinder_25_vmax_mean == 0)) = [];
cylinder_25_halfwidth_fit = polyfit(unique(x_vec_25_cylinder),sort(half_width_mean_25_cylinder),1);
cylinder_25_vmax_fit = polyfit(unique(x_vec_25_cylinder),sort(cylinder_25_vmax_mean),2);

subplot(1,2,1)
scatter(x_vec_25_cylinder,half_width_vec_25_cylinder)
hold on
plot(unique(x_vec_25_cylinder), polyval(cylinder_25_halfwidth_fit,unique(x_vec_25_cylinder)))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Half\:Width\:[mm]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Cylinder\:Half\:Width\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)
sgtitle("$Cylinder\:25\:m/s\:Half\:Width\:and\:Max\:Velocity\:Deficit\:vs\:X-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
scatter(x_vec_25_cylinder,vel_deficit_vec_25_cylinder_max)
hold on 
plot(unique(x_vec_25_cylinder), flip(polyval(cylinder_25_vmax_fit,unique(x_vec_25_cylinder))))
xlabel("$X-Location\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Max\:Velocity\:Deficit\:[m/s]$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Cylinder\:Max\:Deficit\:vs\:X\:Location$",'Interpreter','latex','FontSize',18)

%% Begin the code for the third bullet point
% Need to loop through the vel_deficit_vec_15_airfoil(cylinder) cell arrays
% dividing each vector in the array by the corresponding maximum deficit
% for that vector. Will need to do a check to ensure the vector in the
% deficit cell array is not empty

% First do the 15 m/s airfoil non dimensional velocity deficit and y
% location
counter = 1;
counter_1 = 1;
figure
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1)
for i = 1:length(vel_deficit_vec_15_airfoil)
    if ~isempty(vel_deficit_vec_15_airfoil{i})
        non_dim_vel_deficit_vec_15_airfoil{counter} = vel_deficit_vec_15_airfoil{i}./vel_deficit_vec_15_airfoil_max(counter);
        counter = counter+1;
    end
    if ~isempty(y_loc_vec_15_airfoil{i})
        non_dim_y_loc_vec_15_airfoil{counter_1} = y_loc_vec_15_airfoil{i}./half_width_vec_15_airfoil(counter_1);
        x_loc_vec_15_airfoil_nondimplot{counter_1} = x_loc_vec_15_airfoil{i};
        quiver(x_loc_vec_15_airfoil_nondimplot{counter_1},non_dim_y_loc_vec_15_airfoil{counter_1},non_dim_vel_deficit_vec_15_airfoil{counter_1},zeros(1,length(non_dim_vel_deficit_vec_15_airfoil{counter_1})))
        hold on
        counter_1 = counter_1 + 1;
    end 
end 
xlabel("$X-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Non-Dimensional\:Y-Distance$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Airfoil\:Velocity\:Non-Dimensional\:Deficit$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
sgtitle("$Airfoil\:Non-Dimensional\:Deficit\:vs\:Non-Dimensional\:Y-Location\:Plots$",'Interpreter','latex','FontSize',18)
counter = 1;
counter_1 = 1;
for i = 1:length(vel_deficit_vec_25_airfoil)
    if ~isempty(vel_deficit_vec_25_airfoil{i})
        non_dim_vel_deficit_vec_25_airfoil{counter} = vel_deficit_vec_25_airfoil{i}./vel_deficit_vec_25_airfoil_max(counter);
        counter = counter+1;
    end
    if ~isempty(y_loc_vec_25_airfoil{i})
        non_dim_y_loc_vec_25_airfoil{counter_1} = y_loc_vec_25_airfoil{i}./half_width_vec_25_airfoil(counter_1);
        x_loc_vec_25_airfoil_nondimplot{counter_1} = x_loc_vec_25_airfoil{i};
        quiver(x_loc_vec_25_airfoil_nondimplot{counter_1},non_dim_y_loc_vec_25_airfoil{counter_1},non_dim_vel_deficit_vec_25_airfoil{counter_1},zeros(1,length(non_dim_vel_deficit_vec_25_airfoil{counter_1})))
        hold on
        counter_1 = counter_1 + 1;
    end 
end 
xlabel("$X-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Non-Dimensional\:Y-Distance$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Airfoil\:Velocity\:Non-Dimensional\:Deficit$",'Interpreter','latex','FontSize',18)
%%%%%%%%%%% REPEAT FOR CYLINDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 1;
counter_1 = 1;
figure
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1)
for i = 1:length(vel_deficit_vec_15_cylinder)
    if ~isempty(vel_deficit_vec_15_cylinder{i})
        non_dim_vel_deficit_vec_15_cylinder{counter} = vel_deficit_vec_15_cylinder{i}./vel_deficit_vec_15_cylinder_max(counter);
        counter = counter+1;
    end
    if ~isempty(y_loc_vec_15_cylinder{i})
        non_dim_y_loc_vec_15_cylinder{counter_1} = y_loc_vec_15_cylinder{i}./half_width_vec_15_cylinder(counter_1);
        x_loc_vec_15_cylinder_nondimplot{counter_1} = x_loc_vec_15_cylinder{i};
        quiver(x_loc_vec_15_cylinder_nondimplot{counter_1},non_dim_y_loc_vec_15_cylinder{counter_1},non_dim_vel_deficit_vec_15_cylinder{counter_1},zeros(1,length(non_dim_vel_deficit_vec_15_cylinder{counter_1})),4)
        hold on
        counter_1 = counter_1 + 1;
    end 
end 
xlabel("$X-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Non-Dimensional\:Y-Distance$",'Interpreter','latex','FontSize',18)
title("$15\:m/s\:Cylinder\:Velocity\:Non-Dimensional\:Deficit$",'Interpreter','latex','FontSize',18)
sgtitle("$Cylinder\:Non-Dimensional\:Deficit\:vs\:Non-Dimensional\:Y-Location\:Plots$",'Interpreter','latex','FontSize',18)
subplot(1,2,2)
counter = 1;
counter_1 = 1;
for i = 1:length(vel_deficit_vec_25_cylinder)
    if ~isempty(vel_deficit_vec_25_cylinder{i})
        non_dim_vel_deficit_vec_25_cylinder{counter} = vel_deficit_vec_25_cylinder{i}./vel_deficit_vec_25_cylinder_max(counter);
        counter = counter+1;
    end
    if ~isempty(y_loc_vec_25_cylinder{i})
        non_dim_y_loc_vec_25_cylinder{counter_1} = y_loc_vec_25_cylinder{i}./half_width_vec_25_cylinder(counter_1);
        x_loc_vec_25_cylinder_nondimplot{counter_1} = x_loc_vec_25_cylinder{i};
        quiver(x_loc_vec_25_cylinder_nondimplot{counter_1},non_dim_y_loc_vec_25_cylinder{counter_1},non_dim_vel_deficit_vec_25_cylinder{counter_1},zeros(1,length(non_dim_vel_deficit_vec_25_cylinder{counter_1})),4)
        hold on
        counter_1 = counter_1 + 1;
    end 
end 
xlabel("$X-Distance\:[mm]$",'Interpreter','latex','FontSize',18)
ylabel("$Non-Dimensional\:Y-Distance$",'Interpreter','latex','FontSize',18)
title("$25\:m/s\:Cylinder\:Velocity\:Non-Dimensional\:Deficit$",'Interpreter','latex','FontSize',18)
%% Begin the code for the final bullet point
% First perform the calculation on the 15 m/s airfoil data
for i = 1:length(vel_deficit_vec_15_airfoil)
    if ~isempty(vel_deficit_vec_15_airfoil{i})
        drag_15_airfoil{i} = trapz(y_loc_vec_15_airfoil{i}.*10^(-3), vel_wake_vec_15_airfoil{i}.*vel_deficit_vec_15_airfoil{i});
        Cd_15_airfoil{i} = 2.*drag_15_airfoil{i}./((vel_tunnel_vec_15_airfoil{i}.^2).*.0889);
        mean_cd_15_airfoil(i) = mean(Cd_15_airfoil{i});
    end
end
final_mean_cd_15_airfoil = mean(nonzeros(mean_cd_15_airfoil));
% Repeat for 25 m/s airfoil
for i = 1:length(vel_deficit_vec_25_airfoil)
    if ~isempty(vel_deficit_vec_25_airfoil{i})
        drag_25_airfoil{i} = trapz(y_loc_vec_25_airfoil{i}.*10^(-3), vel_wake_vec_25_airfoil{i}.*vel_deficit_vec_25_airfoil{i});
        Cd_25_airfoil{i} = 2.*drag_25_airfoil{i}./((vel_tunnel_vec_25_airfoil{i}.^2).*.0889);
        mean_cd_25_airfoil(i) = mean(Cd_25_airfoil{i});
    end
end
final_mean_cd_25_airfoil = mean(nonzeros(mean_cd_25_airfoil));
% Now perform calculations for the cylinder
for i = 1:length(vel_deficit_vec_15_cylinder)
    if ~isempty(vel_deficit_vec_15_cylinder{i})
        drag_15_cylinder{i} = trapz(y_loc_vec_15_cylinder{i}./1000, vel_wake_vec_15_cylinder{i}.*vel_deficit_vec_15_cylinder{i});
        Cd_15_cylinder{i} = 2.*drag_15_cylinder{i}./((vel_tunnel_vec_15_cylinder{i}.^2).*(.0127));
        mean_cd_15_cylinder(i) = mean(Cd_15_cylinder{i});
    end
end
final_mean_cd_15_cylinder = mean(nonzeros(mean_cd_15_cylinder));
% Repeat for 25 m/s airfoil
for i = 1:length(vel_deficit_vec_25_cylinder)
    if ~isempty(vel_deficit_vec_25_cylinder{i})
        drag_25_cylinder{i} = trapz(y_loc_vec_25_cylinder{i}.*10^(-3), vel_wake_vec_25_cylinder{i}.*vel_deficit_vec_25_cylinder{i});
        Cd_25_cylinder{i} = 2.*drag_25_cylinder{i}./((vel_tunnel_vec_25_cylinder{i}.^2).*.0127);
        mean_cd_25_cylinder(i) = mean(Cd_25_cylinder{i});
    end
end
final_mean_cd_25_cylinder = mean(nonzeros(mean_cd_25_cylinder));
toc