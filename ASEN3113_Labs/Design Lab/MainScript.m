% Thermo Design Lab 

%% Solar Flux on Radiator over time
clc;clear all; close all;

WinterIR = 88; %W/m^2
SummerIR = 63;
EquinoxIR = (88+63)/2;
EclipseIR = 11;

sigma = 5.67*10^-8;

epsilon = 0.85;
alpha = 0.2;

Ts = 30 + 273.15;

Gs = 1361; 

syms A
eqn = sigma*epsilon*A*(Ts^4) ==  20 + A*EquinoxIR*epsilon + A*Gs*alpha;

Area = solve(eqn,A);
Area = double(Area);
%Area = .2863;

% Area = 0.283


Gs_Winter = 1361/(1-0.0167)^2 * cosd(23.5);
Gs_Summer = 1361/(1+0.0167)^2 * cosd(23.5);
Gs_Equinox = 1361;

time = [0:1:24*3600];

theta = 360/86400;

% flux for each
GsWinter = sind(theta*time)*Gs_Winter*Area;

GsSummer = sind(theta*time)*Gs_Summer*Area;

GsEquinox = sind(theta*time)*Gs_Equinox * Area;


%time = [12:1/3600:24];

time = linspace(0,24,24*3600 + 1);

for i = 1:length(GsWinter)
if GsWinter(i) < 0
    GsWinter(i) = 0;
end
end

for i = 1:length(GsSummer)
if GsSummer(i) < 0
    GsSummer(i) = 0;
end
end

for i = 1:length(GsEquinox)
if GsEquinox(i) < 0
    GsEquinox(i) = 0;
end
end

% Compensate for the eclipse in the Equinox
GsEquinox(39600:43200) = 0;

figure
plot(time,GsWinter)
hold on
plot(time,GsSummer)
plot(time,GsEquinox)
xticks([0 3 6 9 12 15 18 21 24])
title('Solar Flux on Radiator');
xlabel('Time (hours)');
ylabel('Solar Flux (W/m^2)');
legend('Winter','Summer','Equinox');

%% Unheated temperature

IRBack = [88,66,(88+63)/2];

%Gs = [sind(theta*time)*Gs_Winter, sind(theta*time)*Gs_Summer, sind(theta*time)*Gs_Equinox];


Power = 20;
area = .283;

% surface temp 

TsWinter = (((IRBack(1)*epsilon*Area + Power + alpha.*GsWinter)./ (epsilon*sigma*Area)) ).^(0.25) - 273;
TsSummer = (((IRBack(2)*epsilon*Area + Power + alpha.*GsSummer)./ (epsilon*sigma*Area)) ).^(0.25) - 273;
TsEquinox = (((IRBack(3)*epsilon*Area + Power + alpha.*GsEquinox)./ (epsilon*sigma*Area)) ).^(0.25) - 273;



% Surface temp during Eclipse 

% from 11pm to 12pm need to change the Ts for the equinox case onle change
% is the IR backload from the space craft
IRBackEclipse = 11;

TsEquinox(39600:43200) = (((IRBackEclipse*epsilon*Area + Power + alpha.*GsEquinox(39600:43200))./ (epsilon*sigma*Area)) ).^(0.25) - 273;

% Change the IRBackload from the same index as the eclipse for the power
% calculations



% TsEclipse = (QinE./sigma*epsilon*Area).^(.25);


%% Operational heat power drawn

% Need to keep the temperature between 20 and 30 degrees C
Ts = 20 + 273.15;

PowerAdd_Winter = (epsilon*sigma*Area*Ts^4) - (IRBack(1)*epsilon*Area + Power + alpha.*GsWinter) ;

for i = 1:length(GsWinter)
if PowerAdd_Winter(i) < 0
    PowerAdd_Winter(i) = 0;
end
end

PowerAdd_Summer = (epsilon*sigma*Area*Ts^4) - (IRBack(2)*epsilon*Area + Power + alpha.*GsSummer) ;

for i = 1:length(GsSummer)
if PowerAdd_Summer(i) < 0
    PowerAdd_Summer(i) = 0;
end
end

PowerAdd_Equinox = (epsilon*sigma*Area*Ts^4) - (IRBack(3)*epsilon*Area + Power + alpha.*GsEquinox) ;

for i = 1:length(GsEquinox)
if PowerAdd_Equinox(i) < 0
    PowerAdd_Equinox(i) = 0;
end
end

% Change for eclipse case
PowerAdd_Equinox(39600:43200) = (epsilon*sigma*Area*Ts^4) - (IRBackEclipse*epsilon*Area + Power + alpha.*GsEquinox(39600:43200) ) ;


%%  Survival Power
% Need to keep the temperature above -40 degrees C
Ts = -40 + 273.15;

PowerSurvival_Winter = (epsilon*sigma*Area*Ts^4) - (IRBack(1)*epsilon*Area + alpha.*GsWinter) ;

for i = 1:length(GsWinter)
if PowerSurvival_Winter(i) < 0
    PowerSurvival_Winter(i) = 0;
end
end


PowerSurvival_Summer = (epsilon*sigma*Area*Ts^4) - (IRBack(2)*epsilon*Area + alpha.*GsSummer) ;

for i = 1:length(GsSummer)
if PowerSurvival_Summer(i) < 0
    PowerSurvival_Summer(i) = 0;
end
end


PowerSurvival_Equinox = (epsilon*sigma*Area*Ts^4) - (IRBack(3)*epsilon*Area + alpha.*GsEquinox) ;

for i = 1:length(GsEquinox)
if PowerSurvival_Equinox(i) < 0
    PowerSurvival_Equinox(i) = 0;
end
end


% Need to change the IRBack load for the eclipse case 

PowerSurvival_Equinox(39600:43200) = (epsilon*sigma*Area*Ts^4) - (IRBackEclipse*epsilon*Area + alpha.*GsEquinox(39600:43200) ) ;
    
%% Plots

% Winter
figure
grid on
hold on 

yyaxis left
plot(time,TsWinter)
xlabel('Time (hours)');
ylabel('Temperature (\circC)')

xticks([0 3 6 9 12 15 18 21 24])

yyaxis right
plot(time, PowerAdd_Winter)
plot(time, PowerSurvival_Winter,'--k')
xlabel('Time (hours)');
ylabel('Power (W)') 
title('Winter Solstice');
legend('Unheated Temperature','Added Power','Survival Power');



%Summer
figure
grid on
hold on
yyaxis left
plot(time,TsSummer)
xlabel('Time (hours)');
ylabel('Temperature (\circC)')

xticks([0 3 6 9 12 15 18 21 24])

yyaxis right
plot(time, PowerAdd_Summer)
plot(time, PowerSurvival_Summer,'--k')
xlabel('Time (hours)');
ylabel('Power (W)') 
title('Summer Solstice');

legend('Unheated Temperature','Added Power','Survival Power');





% Equinox
figure
grid on
hold on
yyaxis left
plot(time,TsEquinox)
xlabel('Time (hours)');
ylabel('Temperature (\circC)')

xticks([0 3 6 9 12 15 18 21 24])

yyaxis right
plot(time, PowerAdd_Equinox)
plot(time, PowerSurvival_Equinox,'--k')
xlabel('Time (hours)');
ylabel('Power (W)') 
title('Equinox');

legend('Unheated Temperature','Added Power','Survival Power');