%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2004 Bottle Rocket Simulation 3D
% This code exists to simulate the flight of a water/pressure rocket
% built with a 2 liter bottle. In its current form the code is totally self
% contained and will output the verification case, a full variable
% analysis, and then simulate a variety of control variable combinations
% within the critical range of each variable finally outputting the control
% conditions that output the shortest flight that lands with 0.1m of the
% target at 75m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions
% Phase 1:Isentropic and incompressible
% Phase 1.1:Heading is fixed
% Phase 2:Isentropic and compressible
% Phase 3:Ballistic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Trace Valade(102266015) && Delaney Jones(106081404)
% Created: 11/19/18 -- Edited: 3/20/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
global FID;
FID=fopen('YourPermanentRecord.txt','w+');
%% control variable testing
% Verification case conditions
water_vol_test=0.001; % m^3
p_init_test=358840; % Pa
launch_angle_test=pi/4; % radians
drag_coeff_test=0.3; % dimensionless

% Control Variable Analysis set
water_vol=1e-4:1e-4:2e-3; % m^3
p_init=1e5:5e4:10e5; % Pa
launch_at=0:pi/60:pi/2; % radians
drag_coeff=0:0.05:1; % dimensionless

% Vary each control variable while holding the others to the verification
% case
w_size=size(water_vol,2);
hw=zeros(1,w_size);
dw=zeros(1,w_size);
for i=1:w_size % volume varies
    [~,hw(i),dw(i),~]=Rocket_Analysis(water_vol(i),p_init_test,launch_angle_test,drag_coeff_test,0);
end

p_size=size(p_init,2);
hp=zeros(1,p_size);
dp=zeros(1,p_size);
for i=1:p_size % pressure varies
    [~,hp(i),dp(i),~]=Rocket_Analysis(water_vol_test,p_init(i),launch_angle_test,drag_coeff_test,0);
end

l_size=size(launch_at,2);
ht=zeros(1,l_size);
dt=zeros(1,l_size);
for i=1:size(launch_at,2) % launch angle varies
    [~,ht(i),dt(i),~]=Rocket_Analysis(water_vol_test,p_init_test,launch_at(i),drag_coeff_test,0);
end

d_size=size(drag_coeff,2);
hd=zeros(1,d_size);
dd=zeros(1,d_size);
for i=1:d_size % drag coefficient varies
    [~,hd(i),dd(i),~]=Rocket_Analysis(water_vol_test,p_init_test,launch_angle_test,drag_coeff(i),0);
end

close all; % some extra plots are generated during the previous routine
           % they are unnecessary.

% Plot the single variable variation results in max distance and max height
figure
hold on
subplot(2,2,1) % volume variation
plot(water_vol,dw,water_vol,hw)
    title('Water Volume Variation')
    xlabel('Water Volume (m^3)')
    ylabel('(m)')
    legend({'Max Distance','Max Height'},'Location','best');

subplot(2,2,2) % pressure variation
plot(p_init,dp,p_init,hp)
    title('Initial Pressure Variation')
    xlabel('Initial Pressure (Pa)')
    ylabel('(m)')
    
subplot(2,2,3) % launchh angle variation
plot(launch_at*180/pi,dt,launch_at*180/pi,ht)
    title('Launch Angle Variation')
    xlabel('Launch Angle (deg)')
    ylabel('(m)')
    
subplot(2,2,4) % drag coefficient variation
plot(drag_coeff,dd,drag_coeff,hd)
    title('Drag Coefficient Variation')
    xlabel('Drag Coefficient')
    ylabel('(m)')
    
% shrink control variable set to critical range where each variable is both
% realistic and optimalish.
water_vol=5e-4:1e-4:1.5e-3;
p_init=1e5:5e4:6.5e5;
launch_at=pi/6:pi/60:pi/3;
drag_coeff=0.3:0.05:0.5;

% set up parameter variation coordinates
% [V,P,T,D]=ndgrid(water_vol,p_init,launch_at,drag_coeff);
% l=size(V,1)*size(V,2)*size(V,3)*size(V,4); % length of coordinate combination set
% %FT=999999999999999; % max flight time initialization
% FD=-1;
% for i=1:l
%    [~,~,fd,~]=Rocket_Analysis(V(i),P(i),T(i),D(i),1); 
%    if fd>FD % extract current best state based on maximum flight distance of control set
%       FD=fd;
%       Best_State=[V(i),P(i),T(i),D(i)];
%    end
% end
% 
% fprintf(FID,'The Optimal State for minimum flight time is: Time: %5.2f.\n    WaterVol: %0.3f m^3, Pressure: %6.2f Pa, Angle: %2.0f radians, DragCoeff: %0.2f\n',...
%     FD,Best_State(1),Best_State(2),Best_State(3)*180/pi,Best_State(4));
Best_State=[5e-4,275790,0.837758040957278,0.5,0.37];
[~,~,~,path]=Rocket_Analysis(Best_State(1),Best_State(2),Best_State(3),Best_State(4),0);
figure(42)
colormap('jet')
plot3(path(:,1),path(:,2),path(:,3))
[tar_prob,landings]=MCescher(Best_State);

%compute eigenstats of landings
mean_x=mean(landings(:,1));
mean_y=mean(landings(:,2));
covM=cov(landings(:,1),landings(:,2));
[vec,val]=eig(covM);

%prep to plot error ellipis
n=100;
phi=0:1/n:2*pi;
xy=[cos(phi'),sin(phi')]*sqrt(val)*vec';
xy(:,3)=0;

figure(42)
hold on
cmap=[landings(:,1)-path(end,1), landings(:,2)-path(end,2), landings(:,3)-path(end,3)];
C=dot(cmap,cmap,2);
plot3(xy(:,1)+mean_x,xy(:,2)+mean_y,xy(:,3),'r')
plot3(2*xy(:,1)+mean_x,2*xy(:,2)+mean_y,xy(:,3),'y')
plot3(3*xy(:,1)+mean_x,3*xy(:,2)+mean_y,xy(:,3),'g')
scatter3(landings(:,1),landings(:,2),landings(:,3),0.1,C,'filled');
grid on
title('Optimal Flight w/MonteCarlo')
xlabel('Downrange Distance (m)')
ylabel('Crossrange Distance (m)')
zlabel('Height (m)')
%% Monte Carlo Analysis of Best State
%Monte Carlo Setup on Best State
function [tar_prob,landings]=MCescher(Best_State)
global FID;
N = 10000; %Number of tests to run
PT= randn(N,1)*3450+Best_State(2); %  0.5 psi variation assumed
VolWaterInit = randn(N,1)*1e-6+Best_State(1); % .001 liter variation assumed 
LaunchAngle = randn(N,1)*pi/180+Best_State(3); % 1 degree variation assumed
CoeffDrag = randn(N,1)*0.005+Best_State(4); % 0.005 variation assumed 

% initialize vector 
randomD=zeros(N,1);
landings=zeros(N,3);
fprintf("Simulating: [%-20s]\n", ''); %progress bar
for index = 1:N
    [~, ~, randomD(index),fp] = Rocket_Analysis( VolWaterInit(index),PT(index), LaunchAngle(index), CoeffDrag(index),1);
    landings(index,:)=fp(end,:);
    %progress bar animation
    if index ==1
        clc
        progstr = '#';
        fprintf("Simulating: [%-20s]\n", progstr);
    elseif (sum(index/N*20 == 1:1:20)>0)
        clc;
        progstr = '#';
        for j = 2:index/N*20
            progstr = cat(2,progstr,'#');
        end
        fprintf("Simulating: [%-20s]\n", progstr);
    end
end

%graph distribution
figure
histogram(randomD)
    title('Distance Distribution')
    xlabel('Distance (m)')
    ylabel('Amount (#)')

tar_prob=sum(randomD>=74.5 & randomD<=75.5)/N;
fprintf(FID,'Chance of Best State hitting within 0.5m of target is %0.3f\n',tar_prob);
end
%% ROCKET FUNCTION
function [flight_time,maxH,maxD,flight_path]=Rocket_Analysis(water_vol_init,p_init,launch_angle,drag_coeff,sim)
    global FID;
    % ACTUAL SIMULATION FUNCTION: simulates and plots rocket flights based
    % on control variables input, switcher controls gravity projection
    % onto the rail during phase 1.1
    %% STATE INITIALIZATION
    % constants
    mass_bottle=0.15; % kg
    rho_water=1000; % kg/m^3
    vol_bottle=0.002; % kg/m^3
    T_atm=300; % K
    R=287; % J/kgK
    target=75; % m 
    miss_by=0.1; % m 
    
    % initialization
    x=0;
    y=0;
    z=0.25;
    xdot=0;
    ydot=0;
    zdot=0;

    % init calculations
    mass_water=water_vol_init*rho_water; % kg
    T_init=T_atm; % K
    vol_air_init=vol_bottle-water_vol_init; % m^3
    mass_air=p_init*vol_air_init/R/T_init; % kg
    mass_init=mass_water+mass_air+mass_bottle; % kg
    windspeed=[0 0 0 ; 0 0 0];
    
    if sim==1
        windspeed=[0 1 0 ; 0 3 0];
    end
    % initial state for numerical integration
    state=[x y z xdot ydot zdot mass_init p_init vol_air_init].';


    %% ODE45 CALL
    % integration range
    tspan=[0 7];
    
    % option set, restrict relative tolerance and create phase change events  
    opts=odeset('RelTol',1e-5,'Events',@phase_swap);
    
    
    
    % CALL ODE45                       yay
    [t,s_vec,t_event,s_event,ind]=ode45(@(t,s_vec)State_function(s_vec,launch_angle,vol_air_init,p_init,mass_air,drag_coeff,windspeed),tspan,state,opts);
    
   
    
    % extract thrust from flight-state solution
    sv_size=size(s_vec,1);
    thrust=zeros(sv_size,3);
    for i=1:sv_size
        [dstate, thrust(i,:)]=State_function(s_vec(i,:),launch_angle,vol_air_init,p_init,mass_air,drag_coeff,windspeed);    
    end
    
    % extract max height and distance from flight-state solution
    maxH=max(s_vec(:,3));
    maxD=max(s_vec(:,1));
    
    %extract flight path
    flight_path=s_vec(:,1:3);
    
    % find thrust magnitudes
    thrust(:,4)=(thrust(:,1).^2+thrust(:,2).^2+thrust(:,3).^2).^(1/2);
    
    % Plotting routines, this could probably be made a function
    % FOR CONTROL VARIABLE VARIATION
    % if simulation is within range of target plot flight path and thrust
    % profile
    if maxD>=(target-miss_by) && maxD<=(target+miss_by) && sim==0
        figure
        hold on    
        subplot(1,2,1) % thrust profile
        plot(t,thrust(:,4)) 
            title('Thrust Profile')
            xlabel('Time (s)')
            ylabel('Thrust (N)')
            xl=xlim;
            xl(2)=t(find(thrust==0,1));
            xlim(xl);
            
        subplot(1,2,2) % flight path with phase transitions
        plot(s_vec(:,1),s_vec(:,3))
        hold on
            title('Flight Arc')
            xlabel('Distance (m)')
            ylabel('Height (m)')
            yl=ylim;
            yl(1)=0;
            ylim(yl);
            xl=xlim;
            xl(2)=s_vec(find(s_vec(:,2)<=0,1,'last'),1);
            xlim(xl)
            scatter(s_event(2:3,1),s_event(2:3,2),'filled')
            
        % print control variables to screen    
        fprintf(FID,'Max x: %3.2f, Max z: %3.2f, WaterVol: %0.3f, Pressure: %6.2f, Angle: %2.0f, DragCoeff: %0.2f',max(s_vec(:,1)),max(s_vec(:,3)),...
           water_vol_init, p_init,launch_angle*180/pi, drag_coeff );
        fprintf(FID,'\n');
        
        %extract flight time
        flight_time=max(t);
    else %need to assign something to flight time since the analysis function returns it
        flight_time='dontcare';
    end
    
    % initial condition variation of plot routine
    if p_init==4.281644000000000e5 && drag_coeff==0.5 && launch_angle==pi/4 && water_vol_init==0.001
        hold on    
        subplot(1,2,1) % thrust profile
        plot(t,thrust(:,4))
            title('Thrust Profile')
            xlabel('Time (s)')
            ylabel('Thrust (N)')
            xl=xlim;
            xl(2)=t(find(thrust==0,1));
            xlim(xl);
            
        subplot(1,2,2) % Flight path

            plot(s_vec(:,1),s_vec(:,3),'DisplayName','Rail Projection')
            hold on
                title('Flight Arc')
                xlabel('Distance (m)')
                ylabel('Height (m)')
                yl=ylim;
                yl(1)=0;
                ylim(yl);
                legend('show','Location','southoutside')
                scatter(s_event(2:3,1),s_event(2:3,3),'filled','DisplayName','Phase Transitions')
                
        fprintf(FID,'Max x: %3.2f, Max z: %3.2f, WaterVol: %0.3f, Pressure: %6.2f, Angle: %2.0f, DragCoeff: %0.2f',max(s_vec(:,1)),max(s_vec(:,2)),...
           water_vol_init, p_init,launch_angle*180/pi, drag_coeff );
        fprintf('\n');
        flight_time='dontcare';
    end
end
%% STATE FUNCTION
% 2012 Bottle Rocket State Function for ODE45
% State = [x z xdot zdot mass pressure volume_air];
% lots of math below, output is the first time derivative of the state

function [dState,thrust]=State_function(state,init_angle,init_vol,init_P,m_air_i,drag_coeff,windspeed)
    % CONSTANTS
    discharge_coeff=0.8; % dimensionless
    bottle_area=0.008659; % m^2
    throat_area=3.4636e-04; % m^2
    p_atm=83054.2464; % Pa
    rho_w= 1000; %kg/m^3
    gamma=1.4; % dimensionless
    rho_air_atm=0.961; % kg/m^3
    vol_bottle=0.002; % m^3
    %drag_coeff=0.5; % dimensionless
    stand_length=0.5; % m
    R=287; % J/kgK
    stand_height=0.25; % m
    
    % STATE EXTRACTION
    x=state(1);
    y=state(2);
    z=state(3);
    xdot=state(4);
    ydot=state(5);
    zdot=state(6);
    mass=state(7);
    P=state(8);
    Vol_air=state(9);
    
    % phase independent calculations
    rho_init=m_air_i/init_vol; % kg/m^3
    rho=rho_init*(P/init_P)^(1/gamma); % kg/m^3
    m_air=rho*Vol_air; % kg
    velmag=sqrt(xdot^2 + ydot^2 + zdot^2); % m/s
    heading=[xdot/velmag ydot/velmag zdot/velmag]; % direction unit vector   
    g=[0 0 -9.81]; % m/s^2
    WS=linspeed(windspeed(1,:),windspeed(2,:),30,z);
    rel_wind=velmag*heading-WS;
    rel_wind_mag=sqrt(dot(rel_wind,rel_wind));
    
    % pressure at end of phase 1, updated on the fly
    persistent P_end
    
    % heading control on rail
    if sqrt(x^2 + y^2 +(z-stand_height)^2)<stand_length
        heading=[cos(init_angle) 0 sin(init_angle)];           
        g=-9.81*sin(init_angle)*heading;
        rel_wind=velmag*heading;
        rel_wind_mag=sqrt(dot(rel_wind,rel_wind));
    end
    
    % phase dependent calculations
    if P<=p_atm % ballistic phase
        % mass, pressure, and volume are static.
        mdot=0;
        Pdot=0;
        Voldot=0;
        
       % force - acceleration calculations
       thrust=0;
       drag=rho_air_atm*rel_wind_mag*drag_coeff*bottle_area*rel_wind/2;
       accel=-drag/mass+g;
       xdd=accel(1);
       ydd=accel(2);
       zdd=accel(3);
                
    elseif Vol_air>=vol_bottle % thrust w/o water
       % mass and pressure vary, volume is constant, choked flow check must
       % happen as air is compressible in the phase
       T=P/(rho*R); % gas temperature       
       ps=P*(2/(gamma+1))^(gamma/(gamma-1)); % critical choking pressure
       if ps>p_atm % choked flow -- mach 1 exit conditions
           P_exit=ps;
           T_exit=T*(2/(gamma+1));
           vel_exit=sqrt(gamma*R*T_exit);
       else % compressed flow -- mach<1 exit conditions
           P_exit=p_atm;
           Mach_exit=sqrt(((P/p_atm)^((gamma-1)/gamma)-1)*(2/(gamma-1)));
           T_exit=T/(1+((gamma-1)/2)*Mach_exit^2);
           vel_exit=Mach_exit*sqrt(gamma*R*T_exit);
       end
       rho_exit=P_exit/R/T_exit; % exit density
       mdot=-discharge_coeff*rho_exit*throat_area*vel_exit; % dm/dt
       Pdot=mdot*gamma*P_end/m_air_i*(m_air/m_air_i)^(gamma-1); % dP/dm*dm/dt
       Voldot=0;
       
       %force - acceleration calculations
       thrust=(-mdot*vel_exit+(P_exit-p_atm)*throat_area)*heading;
       drag=rho_air_atm*rel_wind_mag*drag_coeff*bottle_area*rel_wind/2;
       accel=thrust/mass-drag/mass+g;
       xdd=accel(1);
       ydd=accel(2);
       zdd=accel(3);
    else % thrust w/water
       vel_exit=sqrt(2*(P-p_atm)/rho_w); % exit velocity
       mdot=-discharge_coeff*rho_w*throat_area*vel_exit; % mass change
       Voldot=discharge_coeff*throat_area*vel_exit; % air volume change
       Pdot=-Voldot*init_P*gamma*init_vol/Vol_air^2*(init_vol/Vol_air)^(gamma-1); % pressure change from volume
       % force - acceleration calculations
       thrust=-mdot*vel_exit*heading;
       thrust = (2*discharge_coeff*throat_area*(P-p_atm)).*heading;
       drag=rho_air_atm*rel_wind_mag*drag_coeff*bottle_area*rel_wind/2;
       accel=thrust/mass-drag/mass-g;
       xdd=accel(1);
       ydd=accel(2);
       zdd=accel(3);
       
       % persistant end pressure update, final value will be final pressure
       % of phase 1
       P_end=P; 
    end
    
    dState=[xdot ydot zdot xdd ydd zdd mdot Pdot Voldot].';
end

%% EVENTS FUNCTION
% Helper function for ode45 call events are set as
% 1) Leaving launch rail
% 2) Phase change 1->2
% 3) Phase change 2->3
% 4) Landing
function [value,isterminal,direction]=phase_swap(t,state)
    %CONSTANTS
    vol_bottle=0.002; % m^3
    p_atm=83426.56; % Pa
    stand_length=0.5; % m
    stand_height=0.25; % m

    % event value state extraction
    x=state(1);
    y=state(2);
    z=state(3);
    P=state(8);
    Vol_air=state(9);
    
    % Event functions
    E1=stand_length-sqrt(x^2+y^2+(z-stand_height)^2);% position at end of rail
    E2=vol_bottle-Vol_air; % volume of water depleted
    E3=P-p_atm; % pressure equalizes
    E4=z; % landing
    
    % function outputs events -- end of integration -- approach gradient
    value=[E1,E2,E3,E4];
    isterminal=[0,0,0,1];
    direction=[-1,-1,-1,-1];
end
%% linear interpolation of windspeed
function ws=linspeed(min_speed,max_speed,delh,currh)
ws=(max_speed-min_speed)/delh*currh;
end