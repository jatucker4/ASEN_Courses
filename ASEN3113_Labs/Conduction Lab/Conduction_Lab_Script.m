clc;
clear all;
close all;
%% Create constants
alpha_steel = 16.2*10000/(500*8000);
alpha_brass = 115*10000/(8500*380);
alpha_al = 130*10000/(960*2810);

%% Parse the data into matrices
steel_data = readmatrix('Steel_20V_185mA.csv');
brass_data = readmatrix('Brass_30V_285mA.csv');
aluminum_data = readmatrix('Aluminum_25V_240mA.csv');

steel_data = steel_data(2:end,1:9);
brass_data = brass_data(2:end,1:9);
aluminum_data = aluminum_data(2:384,1:9);
%% Analytical solution
v_al = 25;
v_brass = 30;
v_steel = 20;

q_dot_al = .240*v_al;
q_dot_brass = .285*v_brass;
q_dot_steel = .185*v_steel;

A = pi*(0.5*2.54)^2;

k_al = 130/100;
k_brass = 115/100;
k_steel = 16.2/100;

H_an_al = q_dot_al/(k_al*A);
H_an_brass = q_dot_brass/(k_brass*A);
H_an_steel = q_dot_steel/(k_steel*A);
%% Create constants
d_1 = 0.5;
pos_vec =[d_1, d_1+0.5 , d_1+1, d_1+1.5, d_1+2, d_1+2.5, d_1+3, d_1+3.5] .* 2.54;
L = pos_vec(end) + 2.54;
% Solving for T_0 in each case
T_0_vec_steel = steel_data(1:end,2:end);
T_0_vec_brass = brass_data(1:end,2:end);
T_0_vec_aluminum = aluminum_data(1:end,2:end);

p_steel_vec = polyfit(pos_vec, T_0_vec_steel(end,:),1);

p_brass_vec = polyfit(pos_vec, T_0_vec_brass(end,:),1);

p_aluminum_vec = polyfit(pos_vec, T_0_vec_aluminum(end,:),1);

H_steel = p_steel_vec(1);
H_brass = p_brass_vec(1);
H_aluminum = p_aluminum_vec(1);

T_0_steel = p_steel_vec(2);
T_0_brass = p_brass_vec(2);
T_0_aluminum = p_aluminum_vec(2);

steel_line_exp = H_steel.*pos_vec + T_0_steel;
brass_line_exp = H_brass.*pos_vec + T_0_brass;
aluminum_line_exp = H_aluminum.*pos_vec + T_0_aluminum;

steel_line_an = H_an_steel.*pos_vec + T_0_steel;
brass_line_an = H_an_brass.*pos_vec + T_0_brass;
aluminum_line_an = H_an_al.*pos_vec + T_0_aluminum;
% Experimental slope with steady state temp distribution
figure(1)
plot(pos_vec,steel_line_exp)
hold on
plot(pos_vec,steel_line_an)
hold on
plot(pos_vec,T_0_vec_steel(end,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Analytical\:and\:Experimental\:Slopes\:for\:Steel$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Analytical Slope", "Experimental Temp Distribution")

figure(2)
plot(pos_vec,brass_line_exp)
hold on
plot(pos_vec,brass_line_an)
hold on
plot(pos_vec,T_0_vec_brass(end,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Analytical\:and\:Experimental\:Slopes\:for\:Brass$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Analytical Slope", "Experimental Temp Distribution")

figure(3)
plot(pos_vec,aluminum_line_exp)
hold on
plot(pos_vec,aluminum_line_an)
hold on
plot(pos_vec,T_0_vec_aluminum(end,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Analytical\:and\:Experimental\:Slopes\:for\:Aluminum$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Analytical Slope", "Experimental Temp Distribution")
%% Analytical solution over time for brass
t = 1:5840;
% syms n
n = 1;
for j = 1:length(pos_vec)
    for i = 1:length(t)
        sum = (T_0_brass + H_an_brass*pos_vec(j)) ;
        sum1 = 0;
        temp = -8*H_an_brass*L;
        temp2 = (pi*pos_vec(j))/(2*L);   
        temp3 = (pi)/(2*L);
        for k = 1:length(n)
            u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
            sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_brass * t(i));
            sum1 = sum1+u_2_temp;
        end
        u_2(i) = sum+sum1;
    end
    brass_final_temps(j,:) = double(u_2(end));
    
    figure(4)
    plot(t,u_2)
    hold on
    plot(brass_data(1:end,1),T_0_vec_brass(1:end,j),'--')
    xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
    ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
    title("$Analytical\:Solution\:For\:Each\:Thermocouple\:Over\:Time\:for\:Brass$",'Interpreter','latex','FontSize',26)
end
legend("Th1", "ExpTh1","Th2", "ExpTh2","Th3","ExpTh3","Th4","ExpTh4","Th5","ExpTh5","Th6","ExpTh6","Th7","ExpTh7","Th8","ExpTh8")
percent_diff_brass = 100*abs((T_0_vec_brass(end,:)') - brass_final_temps) ./brass_final_temps;
fprintf("Error in the steady state solutions for Brass is: %f\n",mean(percent_diff_brass))
%%  Analytical solution over time for aluminum
t = 1:3820;
% syms n
n = 1;
for j = 1:length(pos_vec)
    for i = 1:length(t)
        sum = (T_0_aluminum + H_aluminum*pos_vec(j)) ;
        sum1 = 0;
        temp = -8*H_aluminum*L;
        temp2 = (pi*pos_vec(j))/(2*L);   
        temp3 = (pi)/(2*L);
        for k = 1:length(n)
            u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
            sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_al * t(i));
            sum1 = sum1+u_2_temp;
        end
        u_2_al(i) = sum+sum1;
    end
    aluminum_final_temps(j,:) = double(u_2_al(end));
    
    figure(5)
    plot(t,u_2_al)
    hold on
    plot(aluminum_data(1:end,1),T_0_vec_aluminum(1:end,j),'--')
    xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
    ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
    title("$Analytical\:Solution\:For\:Each\:Thermocouple\:Over\:Time\:for\:Aluminum$",'Interpreter','latex','FontSize',26)
end
legend("Th1", "ExpTh1","Th2", "ExpTh2","Th3","ExpTh3","Th4","ExpTh4","Th5","ExpTh5","Th6","ExpTh6","Th7","ExpTh7","Th8","ExpTh8")
percent_diff_aluminum = 100*abs((T_0_vec_aluminum(end,:)') - aluminum_final_temps) ./aluminum_final_temps;
fprintf("Error in the steady state solutions for Aluminum is: %f\n",mean(percent_diff_aluminum))
%%  Analytical solution over time for steel
t = 1:14540;
% syms n
n = 1;
for j = 1:length(pos_vec)
    for i = 1:length(t)
        sum = T_0_steel + H_steel*pos_vec(j) ;
        sum1 = 0;
        temp = -8*H_steel*L;
        temp2 = (pi*pos_vec(j))/(2*L);   
        temp3 = (pi)/(2*L);
        for k = 1:length(n)
            u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
            sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_steel * t(i));
            sum1 = sum1+u_2_temp;
        end
        u_2_steel(i) = sum+sum1;
    end
    steel_final_temps(j,:) = double(u_2_steel(end));
    
    figure(6)
    plot(t,u_2_steel)
    hold on
    plot(steel_data(1:end,1),T_0_vec_steel(1:end,j),'--')
    xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
    ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
    title("$Analytical\:Solution\:For\:Each\:Thermocouple\:Over\:Time\:for\:Steel$",'Interpreter','latex','FontSize',26)
end
legend("Th1", "ExpTh1","Th2", "ExpTh2","Th3","ExpTh3","Th4","ExpTh4","Th5","ExpTh5","Th6","ExpTh6","Th7","ExpTh7","Th8","ExpTh8")
percent_diff_steel = 100*abs((T_0_vec_steel(end,:)') - steel_final_temps) ./steel_final_temps;
fprintf("Error in the steady state solutions for Steel is: %f\n",mean(percent_diff_steel))
%% Question 3: Slope of a line across the different steady state cold temps
p_cold_steel_vec = polyfit(pos_vec, T_0_vec_steel(1,:),1);

p_cold_brass_vec = polyfit(pos_vec, T_0_vec_brass(1,:),1);

p_cold_aluminum_vec = polyfit(pos_vec, T_0_vec_aluminum(49,:),1);

H_cold_steel = p_cold_steel_vec(1);
H_cold_brass = p_cold_brass_vec(1);
H_cold_aluminum = p_cold_aluminum_vec(1);

T_cold_steel = p_cold_steel_vec(2);
T_cold_brass = p_cold_brass_vec(2);
T_cold_aluminum = p_cold_aluminum_vec(2);

steel_cold_line_exp = H_cold_steel.*pos_vec + T_cold_steel;
brass_cold_line_exp = H_cold_brass.*pos_vec + T_cold_brass;
aluminum_cold_line_exp = H_cold_aluminum.*pos_vec + T_cold_aluminum;

% Experimental slope with steady state temp distribution
figure(7)
plot(pos_vec,steel_cold_line_exp)
hold on
plot(pos_vec,T_0_vec_steel(1,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Cold\:Steady\:State\:Experimental\:Slopes\:for\:Steel$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Experimental Temp Distribution")

figure(8)
plot(pos_vec,brass_cold_line_exp)
hold on
plot(pos_vec,T_0_vec_brass(1,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Cold\:Steady\:State\:Experimental\:Slopes\:for\:Brass$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Experimental Temp Distribution")

figure(9)
plot(pos_vec,aluminum_cold_line_exp)
hold on
plot(pos_vec,T_0_vec_aluminum(49,:))
xlabel("$Thermocouple\:Position\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Cold\:Steady\:State\:Experimental\:Slopes\:for\:Aluminum$",'Interpreter','latex','FontSize',26)
legend("Experimental Slope", "Experimental Temp Distribution")

%% Question 4: Finding the optimal alpha values for brass at thermocouple 1
tolerance = 0.25;
alpha_brass_vec = linspace(alpha_brass,2*alpha_brass);
u2_brass_temp = [0];
m = 1;
while abs(u2_brass_temp(end) - T_0_vec_brass(end,1)) > tolerance
    t = 1:5840;
    % syms n
    n = 1;
        for i = 1:length(t)
            sum = (T_0_brass + H_brass*pos_vec(1)) ;
            sum1 = 0;
            temp = -8*H_brass*L;
            temp2 = (pi*pos_vec(1))/(2*L);   
            temp3 = (pi)/(2*L);
            for k = 1:length(n)
                u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
                sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_brass_vec(m) * t(i));
                sum1 = sum1+u_2_temp;
            end
            u2_brass_temp(i) = sum+sum1;
        end
    m = m + 1;
end
fprintf("Optimal alpha for brass is: %f\n",alpha_brass_vec(m))

%% Question 4: Finding the optimal alpha values for aluminum
tolerance = 0.25;
alpha_al_vec = linspace(alpha_al,2*alpha_al);
u2_aluminum_temp = [0];
m = 1;
while abs(u2_aluminum_temp(end) - T_0_vec_aluminum(end,1)) > tolerance
    t = 1:3820;
    % syms n
    n = 1;
        for i = 1:length(t)
            sum = (T_0_aluminum + H_aluminum*pos_vec(1)) ;
            sum1 = 0;
            temp = -8*H_aluminum*L;
            temp2 = (pi*pos_vec(1))/(2*L);   
            temp3 = (pi)/(2*L);
            for k = 1:length(n)
                u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
                sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_al_vec(m) * t(i));
                sum1 = sum1+u_2_temp;
            end
            u2_aluminum_temp(i) = sum+sum1;
        end
    m = m + 1;
end
fprintf("Optimal alpha for aluminum is: %f\n",alpha_al_vec(m))

%% Question 4: Finding the optimal alpha values for steel
tolerance = 0.25;
alpha_steel_vec = linspace(alpha_steel,2*alpha_steel);
u2_steel_temp = [0];
m = 1;
while abs(u2_steel_temp(end) - T_0_vec_steel(end,1)) > tolerance
    t = 1:14540;
    % syms n
    n = 1;
        for i = 1:length(t)
            sum = (T_0_steel + H_steel*pos_vec(1)) ;
            sum1 = 0;
            temp = -8*H_steel*L;
            temp2 = (pi*pos_vec(1))/(2*L);   
            temp3 = (pi)/(2*L);
            for k = 1:length(n)
                u_2_temp = ((temp*((-1)^(n(k) + 1)))/(((2*n(k) - 1)*pi)^2))* ...
                sin((2*n(k) - 1)*temp2) * exp(-(((2*n(k) - 1)*temp3)^2) * alpha_steel_vec(m) * t(i));
                sum1 = sum1+u_2_temp;
            end
            u2_steel_temp(i) = sum+sum1;
        end
    m = m + 1;
end
fprintf("Optimal alpha for steel is: %f\n",alpha_steel_vec(m))
%% Question 5
% Looking at the plots the time to steady states are:
tss_brass = brass_data(find(T_0_vec_brass == T_0_vec_brass(end,1),1),1);
tss_al = aluminum_data(find(T_0_vec_aluminum == T_0_vec_aluminum(end,1),1),1);

L = 19.05;

exp_brass_0 = tss_brass*alpha_brass / (L^2);
exp_al_0 = tss_al*alpha_al/ (L^2);
% See how these exponent values change with length
L_vec = linspace(L,(L+2));
alpha_al_vec = linspace(alpha_al,(alpha_al+0.2));
alpha_brass_vec = linspace(alpha_brass,(alpha_brass+0.2));

tss_vary_L_al = (exp_al_0.*(L_vec.^2))./alpha_al;
tss_vary_alpha_al = (exp_al_0.*(L.^2))./alpha_al_vec;

tss_vary_L_brass = (exp_brass_0.*(L_vec.^2))./alpha_brass;
tss_vary_alpha_brass = (exp_brass_0.*(L.^2))./alpha_brass_vec;
% Plot the results
figure(10)
plot(L_vec,tss_vary_L_al)
xlabel("$Rod\:Length\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Time\:to\:Steady\:State[s]$",'Interpreter','latex','FontSize',26)
title("$Time\:to\:Steady\:State\:vs\:Variation\:in\:Rod\:Length-Aluminum$",'Interpreter','latex','FontSize',26)

figure(11)
plot(L_vec,tss_vary_L_brass)
xlabel("$Rod\:Length\:[cm]$",'Interpreter','latex','FontSize',26)
ylabel("$Time\:to\:Steady\:State[s]$",'Interpreter','latex','FontSize',26)
title("$Time\:to\:Steady\:State\:vs\:Variation\:in\:Rod\:Length-Brass$",'Interpreter','latex','FontSize',26)

figure(12)
plot(alpha_al_vec,tss_vary_alpha_al)
xlabel("$Thermal\:Diffusivity\:[\frac{cm^2}{s}]$",'Interpreter','latex','FontSize',26)
ylabel("$Time\:to\:Steady\:State[s]$",'Interpreter','latex','FontSize',26)
title("$Time\:to\:Steady\:State\:vs\:Variation\:in\:Thermal\:Diffusivity-Aluminum$",'Interpreter','latex','FontSize',26)

figure(13)
plot(alpha_brass_vec,tss_vary_alpha_brass)
xlabel("$Thermal\:Diffusivity\:[\frac{cm^2}{s}]$",'Interpreter','latex','FontSize',26)
ylabel("$Time\:to\:Steady\:State[s]$",'Interpreter','latex','FontSize',26)
title("$Time\:to\:Steady\:State\:vs\:Variation\:in\:Thermal\:Diffusivity-Brass$",'Interpreter','latex','FontSize',26)