function [isp,delta_V] = isp_process(temp,pressure,water_mass,...
    sample_rate, filename)
%ISP_PROCESS Summary of this function goes here
%   Detailed explanation goes here
%Convert rate to 1/s
sample_rate = sample_rate*1000;
%Convert water mass to kg
water_vol = water_mass/(1*10^6);
water_mass = water_mass/1000;
%Convert pressure to Pa
pressure = pressure * 6894.76;
%Convert temp to kelvin
temp = temp + 273.15;

% Pull data from file
data = load(filename);
% get thrust in Newtons
thrust = data(:,3) * 4.44822;
[value, index] = max(thrust);
pre_fire = find(thrust(1:index) < 0.1 ,1,'last');
thrust = thrust(pre_fire:end);
[~,min_thrust] = min(thrust);
post_fire = find(thrust<0,1,'first');
thrust = thrust(1:min_thrust);
% thrust = thrust(1936:2315);
time = linspace(0,length(thrust),length(thrust)) ./sample_rate;

plot(time,thrust)
isp = trapz(time,thrust)/(water_mass*9.81);
% Get delta V
R = 287;
P_initial = pressure + 83426.56;
m_bottle = 0.15;
Vol_air_initial = .002 - water_vol;
m_air_initial = (P_initial*Vol_air_initial)/(R*temp);
m_initial_total = water_mass + m_bottle + m_air_initial;
delta_V = isp*9.81*log(m_initial_total/m_bottle);
end

