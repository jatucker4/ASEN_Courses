function [Nnominal] = CA3Prob2Function()
% CA3Prob2Function.m completes all of the computation and plotting required
% for question 2

%CA3Prob2Function.m Evaluates the nominal number of panels required the
%vortex panel method to produce a desired level of accuracy. This error
%analysis method is specifically discussed below in the function. Then,
%using the computed nominal number of panels, the function computes the
%lifting flow over a NACA0012 at four different angles of attack. For each
%angle of attack a plot of c_l vs. alpha and Cp vs. X/c is produced. 

% Author: Michael Martinson
% Collaborators: Trace Valade
% Date: 3/28/2020

%% Airfoils 
%Define required parameters to compute x and y values for all of the 
%required airfoils for later problems: 

%NACA 0012
m0012 = 0/100; 
p0012 = 0/10; 
t0012 = 12/100; 
Digits0012 = 0012; 

%% Error Analysis Section: 
%Error analysis was conducted by considering the c_l value computed for the
%NACA0012 at 0deg angle of attack for increasing values of N. Because the
%NACA0012 is symmetric, a perfectly accurate computation would yield c_l =
%0 for the NACA0012 at 0deg = alpha. 

%The desired quantitative level of accuracy for the analysis was set as the
%number of panels required to produce a c_l result under 0.02 from 0. The N
% determined to meet this requirement would be considered the "nominal" N
% value. 

%Considering c_l vs N for alpha = 0 deg.: 
alphadeg = 0; %deg
alpharad = alphadeg; %rad

% Assign parameters that apply to all test cases: 
V_inf = 20; %m/s
CpSwitch = 0; %Cp plotting not used
c = 1; %m 

%Assign vector of N values to test:
%Use N from 20 to 640 in increments of 20 panels:
Nerrorvec = [20:20:640]; 

%For each value of N considered: 
for i = 1:length(Nerrorvec)

% NACA 0012: Assign the x and y coordinates for the boundary points:
[x0012,y0012] = NACA_Airfoils(m0012,p0012,t0012,c,Nerrorvec(i));

%Fix the NaN at 0 in x0012 and y0012 issue (both NaN values should be 0): 
xfixindex = find(isnan(x0012) == 1); 
yfixindex = find(isnan(y0012) == 1); 
x0012(xfixindex) = 0; %m
y0012(yfixindex) = 0; %m

% Compute the c_l:
c_lErrorVec(i,1) =...
    Vortex_Panel(x0012,y0012,V_inf,alpharad,CpSwitch,Digits0012,0); 

end 

%Change in error vs. N plotting: 
figure; 
plot(Nerrorvec,c_lErrorVec,'r')
hold on; 

%Plot the tolerance line at c_l = 0.02: 
plot([0, Nerrorvec],0.02.*ones(length(Nerrorvec)+1,1),'b'); 

%Label: 
title({'Prob. 2 Error Analysis: $c_l$ deviation from 0;',...
    'NACA0012, $\alpha$ = 0$^{\circ}$'},...
    'FontSize',14)
xlabel('N (Total Panels)','FontSize',14)
ylabel('$c_l$','FontSize',14)
grid on; grid minor; 
xlim([0 max(Nerrorvec)]); 

%Find nominal N where the c_l is under 0.02: 
Nindex = min(find(c_lErrorVec <= 0.02)); 
Nnominal = Nerrorvec(Nindex);

%Add this point clearly marked to the plot: 
hold on; 
scatter(Nnominal,c_lErrorVec(Nindex),15,'filled','g'); 
legend('$c_l$ vs. N','$c_l$ = 0.02','$N_{Nominal}$'); 

%Print the nominal panel count to the command window:
fprintf('\n\nProblem 2 Results:\n'); 
fprintf('Panels required for c_l under 0.02 from 0 for NACA0012...\n')
fprintf('at alpha = 0 deg. = %d\n\n',Nnominal)

%% Using the nominal number of panels, compute lifting flow over a NACA 
%0012 airfoil and plot the results for a range of alpha values: 

%Assign a vector of alpha values to use: 
alphavecdeg = [-5 0 5 10]'; %deg
alphavec = alphavecdeg.*(pi/180); %rad

% Assign parameters that apply to all cases: 
V_inf = 20; %m/s
CpSwitch = 2; %Use problem 2 plotting setting 
c = 1; %m 

%Assign the N value to use and the nominal number of panels: 
N2 = Nnominal; 

% NACA 0012: Assign the x and y coordinates for the boundary points:
[x0012,y0012] = NACA_Airfoils(m0012,p0012,t0012,c,N2);

%Fix the NaN at 0 in x0012 and y0012 issue (both NaN values should be 0): 
xfixindex = find(isnan(x0012) == 1); 
yfixindex = find(isnan(y0012) == 1); 
x0012(xfixindex) = 0; %m
y0012(yfixindex) = 0; %m

%% For every value of alpha considered, find c_l as a function of alpha
%IN RAD for each airfoil; also conduct the Cp plotting for each alpha: 

figure; %Initialize the figure for the Cp plotting
for i = 1:length(alphavec) 
% NACA 0012
c_l0012(i,1) =...
    Vortex_Panel(x0012,y0012,V_inf,alphavec(i),CpSwitch,Digits0012,i); 
end

%% Finalize Plotting Results: Cp
%Add a title to the figure: 
sgtitle(...
 ['Prob. 2: $C_p$ vs. $\frac{x}{c}$, Assorted $\alpha$; NACA0012, N = ',...
 num2str(Nnominal)],'FontSize',14); 

%% Plotting Results: c_l 
%Plot each c_l vs. alpha data point: 
figure;
for i = 1:length(alphavec)
hold on; 
scatter(alphavecdeg(i),c_l0012(i),25,'filled')
end

%labeling: 
xlabel('$\alpha^{\circ}$','FontSize',14); 
ylabel('$c_l$','FontSize',14); 
title(['Prob. 2: $c_l$ vs. $\alpha$; NACA0012, N = ',num2str(Nnominal)],...
    'FontSize',14); 
grid on; grid minor; 

%Mark the x and y axes: 
hold on; 
xline(0,'k','LineWidth',0.5); 
plot(alphavecdeg,zeros(length(alphavecdeg),1),'k','LineWidth',0.5); 

%legend: 
legend('$\alpha$ = -5$^{\circ}$','$\alpha$ = 0$^{\circ}$',...
    '$\alpha$ = 5$^{\circ}$','$\alpha$ = 10$^{\circ}$',...
    'Location','NorthWest')

end 