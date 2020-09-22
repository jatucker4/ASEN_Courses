function []= CA3Prob3Function(N3)
%CA3Prob3Function.m performs all of the computation, plotting, and printing
%required in question 3. 

%CA3Prob3Function.m is used to obtain plots of c_l vs. alpha for NACA
%0012, 2412, and 2424 airfoils. The computations to obtain these plots are
%then used to estimate the lift slope and zero lift angle of attack for
%each airfoil; additional computations are then performed to compare these
%results to predictions from thin airfoil theory (TAT).

% Author: Michael Martinson
% Collaborators: Trace Valade
% Date: 3/27/2020

%% Airfoils 
%Define required parameters to compute x and y values for all of the 
%required airfoils for later problems: 

%NACA 0012
m0012 = 0/100; 
p0012 = 0/10; 
t0012 = 12/100; 
Digits0012 = 0012; 

%NACA 2412
m2412 = 2/100;
p2412 = 4/10; 
t2412 = 12/100;
Digits2412 = 2412; 

%NACA 4412
m4412 = 4/100;  
p4412 = 4/10;
t4412 = 12/100;
Digits4412 = 4412; 

%NACA 2424
m2424 = 2/100;
p2424 = 4/10;
t2424 = 24/100;
Digits2424 = 2424; 

%%
%Assign a vector of alpha values to use: 
%Alpha values between -6 and 15 deg in increments of 1deg. 
alphavecdeg = [-6:1:15]'; %deg
alphavec = alphavecdeg.*(pi/180); %rad

% Assign parameters that apply to all cases: 
V_inf = 20; %m/s
CpSwitch = 0; 
c = 1; %m 

%% Using N3 in each case, find the x and y points for each airfoil: 
%N3 is the number of panels set to be used for all question 3 analysis. 

%% NACA 0012
[x0012,y0012] = NACA_Airfoils(m0012,p0012,t0012,c,N3);

%Fix the NaN at 0 in x0012 and y0012 issue (both NaN values should be 0): 
xfixindex = find(isnan(x0012) == 1); 
yfixindex = find(isnan(y0012) == 1); 
x0012(xfixindex) = 0; 
y0012(yfixindex) = 0; 

%% NACA 2412
[x2412,y2412] = NACA_Airfoils(m2412,p2412,t2412,c,N3);

%% NACA 4412
[x4412,y4412] = NACA_Airfoils(m4412,p4412,t4412,c,N3);

%% NACA 2424
[x2424,y2424] = NACA_Airfoils(m2424,p2424,t2424,c,N3);

%% For every value of alpha considered, find c_l as a function of alpha
%IN RAD for each airfoil:
for i = 1:length(alphavec) 
    
% NACA 0012
c_l0012(i,1) =...
    Vortex_Panel(x0012,y0012,V_inf,alphavec(i),CpSwitch,Digits0012,0); 

% NACA 2412
c_l2412(i,1) =...
    Vortex_Panel(x2412,y2412,V_inf,alphavec(i),CpSwitch,Digits2412,0); 

% NACA 4412
c_l4412(i,1) =...
    Vortex_Panel(x4412,y4412,V_inf,alphavec(i),CpSwitch,Digits4412,0); 

% NACA 2424
c_l2424(i,1) =...
    Vortex_Panel(x2424,y2424,V_inf,alphavec(i),CpSwitch,Digits2424,0); 

end 

%% Plot analysis for each airfoil: 

%% Determine the lift curve slopes:
%% Fit a linear regression to each c_l vs. alpha curve: 
%Assume that the slope of the linear regression provides an effective
%approximation of the lift curve slope in each case: 

% NACA 0012
PolyRes0012 = polyfit(alphavec,c_l0012,1); 

% NACA 2412
PolyRes2412 = polyfit(alphavec,c_l2412,1); 

% NACA 4412
PolyRes4412 = polyfit(alphavec,c_l4412,1); 

% NACA 2424
PolyRes2424 = polyfit(alphavec,c_l2424,1); 

%% Identify the approximate lift curve slope for each case: 

% NACA 0012
slope0012 = PolyRes0012(1); %/rad

% NACA 2412
slope2412 = PolyRes2412(1); %/rad

% NACA 4412
slope4412 = PolyRes4412(1); %/rad

% NACA 2424
slope2424 = PolyRes2424(1); %/rad

%% Determine the approx. zero lift angles of attack:
%Set cl to 0, solve for alpha using the linear regression equation for each
%airfoil case: 

% NACA 0012
alphaL00012 = -PolyRes0012(2)/PolyRes0012(1)*(180/pi); %deg

% NACA 2412
alphaL02412 = -PolyRes2412(2)/PolyRes2412(1)*(180/pi); %deg

% NACA 4412
alphaL04412 = -PolyRes4412(2)/PolyRes4412(1)*(180/pi); %deg

% NACA 2424
alphaL02424 = -PolyRes2424(2)/PolyRes2424(1)*(180/pi); %deg


%% Compare with thin airfoil theory (TAT)
%Lift curve slope for TAT: 
slopeTAT = 2*pi; %/rad

%Zero lift alpha for TAT:  
%Sym airfoil TAT zero lift angle of attack (NACA0012): 
alphaL0TAT0012 = 0; %deg

%% Cambered airfoil TAT zero lift angles of attack:
%(apply Anderson equation 4.61); to do so compute theta and dycdc in terms
%of theta vectors for each cambered airfoil: 

%% NACA 2412
[thetavals2412, dycdx2412] = Prob3NACA_Airfoils(m2412,p2412,t2412,c,N3);

%Assign the function to integrate: 
intgr2412 = dycdx2412.*(cos(thetavals2412 - 1)); 

%Compute the zero lift angle of attack prediction: 
alphaL0TAT2412 = -(1/pi).*trapz(thetavals2412,intgr2412)*(180/pi); 

%% NACA 4412
[thetavals4412, dycdx4412] = Prob3NACA_Airfoils(m4412,p4412,t4412,c,N3);

%Assign the function to integrate: 
intgr4412 = dycdx4412.*(cos(thetavals4412 - 1)); 

%Compute the zero lift angle of attack prediction: 
alphaL0TAT4412 = -(1/pi).*trapz(thetavals4412,intgr4412)*(180/pi); 

%% NACA 2424
[thetavals2424, dycdx2424] = Prob3NACA_Airfoils(m2424,p2424,t2424,c,N3); 

%Assign the function to integrate: 
intgr2424 = dycdx2424.*(cos(thetavals2424 - 1)); 

%Compute the zero lift angle of attack prediction: 
alphaL0TAT2424 = -(1/pi).*trapz(thetavals2424,intgr2424)*(180/pi); 

%% Plotting Problem 3 Results:  
%Compare the results for each airfoil and the TAT results: 

%% Airfoil results: 

%Plot the computed results for each airfoil:
figure;
subplot(2,1,1)

% NACA 0012
plot(alphavecdeg,c_l0012,'b')

% NACA 2412
hold on; 
plot(alphavecdeg,c_l2412,'r')

% NACA 4412
hold on;
plot(alphavecdeg,c_l4412,'m')

% NACA 2424
hold on;
plot(alphavecdeg,c_l2424,'g')

%labeling: 
xlabel('$\alpha^{\circ}$','FontSize',14); 
ylabel('$c_l$','FontSize',14); 
title(['Prob. 3: $c_l$ vs. $\alpha$; Assorted Airfoils, N = ',...
    num2str(N3)],'FontSize',14); 
grid on; grid minor; 

%Mark the x and y axes: 
hold on; 
xline(0,'k','LineWidth',1.5); 
plot(alphavecdeg,zeros(length(alphavecdeg),1),'k','LineWidth',1.5); 

%Additional plot labeling: 
xlim([min(alphavecdeg) max(alphavecdeg)]); 
legend('NACA 0012','NACA 2412','NACA 4412','NACA 2424',...
    'Location','Northwest'); 

%% Airfoil Results compared to TAT: 
% compute cl vs alpha curves according to TAT for each airfoil: 
%(use Anderson eq 4.60)

%NACA0012: 
c_lvsalphaTAT0012 = 2*pi.*(alphavec - 0); 

%NACA2412: 
c_lvsalphaTAT2412 = 2*pi.*(alphavec - alphaL0TAT2412.*(pi/180)); 

%NACA4412: 
c_lvsalphaTAT4412 = 2*pi.*(alphavec - alphaL0TAT4412.*(pi/180)); 

%NACA2424: 
c_lvsalphaTAT2424 = 2*pi.*(alphavec - alphaL0TAT2424.*(pi/180)); 

%% Plotting portion: 

%Create a second plot where the result can be compared with their TAT
%counterparts: 
subplot(2,1,2)
% NACA 0012
plot(alphavecdeg,c_l0012,'b')
hold on; 
plot(alphavecdeg,c_lvsalphaTAT0012,'--b')

% NACA 2412
hold on; 
plot(alphavecdeg,c_l2412,'r')
hold on; 
plot(alphavecdeg,c_lvsalphaTAT2412,'--r')

% NACA 4412
hold on;
plot(alphavecdeg,c_l4412,'m')
hold on; 
plot(alphavecdeg,c_lvsalphaTAT4412,'--m')

% NACA 2424
hold on;
plot(alphavecdeg,c_l2424,'g')
hold on; 
plot(alphavecdeg,c_lvsalphaTAT2424,'--g')

%labeling: 
xlabel('$\alpha^{\circ}$','FontSize',14); 
ylabel('$c_l$','FontSize',14); 
title(['Prob. 3: $c_l$ vs. $\alpha$; Assorted Airfoils, N = ',...
    num2str(N3),'; Comparison with T.A.T'],'FontSize',14); 
grid on; grid minor; 

%Mark the x and y axes: 
hold on; 
xline(0,'k','LineWidth',1.5); 
plot(alphavecdeg,zeros(length(alphavecdeg),1),'k','LineWidth',1.5); 

%Additional plot labeling: 
xlim([min(alphavecdeg) max(alphavecdeg)]); 
legend('NACA 0012','NACA 0012 T.A.T',...
'NACA 2412','NACA 2412 T.A.T',...
'NACA 4412','NACA 4412 T.A.T',...
'NACA 2424','NACA 2424 T.A.T',...
'Location','Northwest');

%Overall title: 
sgtitle('Prob. 3 Results','FontSize',14);

%% Printing results per airfoil case: assign values to print: 
%Vector of airfoil digits: 
Digitsvec = [0012 2412 4412 2424]'; 

%Vector of slope estimates (/rad)
slopevec = [slope0012; slope2412; slope4412; slope2424]; %/rad

%Vector of zero lift angle of attack: 
alphaL0vec = [alphaL00012 alphaL02412 alphaL04412 alphaL02424]'; %deg

%TAT vector of zero lift angle of attack: 
alphaL0TATvec = [0 alphaL0TAT2412 alphaL0TAT4412 alphaL0TAT2424]'; %deg

%% Print the problem 3 comparison with thin airfoil theory results:
fprintf('\nProblem 3 Results:')
%For each airfoil: 
for i = 1:length(Digitsvec)
    
    %Print which airfoil the results are for: 
if i == 1
    fprintf('\n\nNACA 00%s Results:\n',num2str(Digitsvec(i))); 
else
    fprintf('\n\nNACA %s Results:\n',num2str(Digitsvec(i))); 
end 

%Print lift curve results and comparison with TAT:
fprintf('Lift Curve Slope Estimate: %f /rad.\n',slopevec(i))
fprintf('Thin Airfoil Theory Lift Curve Slope: %f /rad.\n',slopeTAT)
fprintf('Percent Error From TAT: %f%%',...
    (abs(slopevec(i) - slopeTAT)/(slopeTAT))*100)

%Print zero lift angle of attack results and comparison with TAT: 
fprintf('\n\nZero Lift Angle of Attack Estimate: %f deg.\n',alphaL0vec(i))
fprintf('Thin Airfoil Theory Estimate: %f deg.\n',alphaL0TATvec(i))
fprintf('Absolute Error from TAT Estimate: %f deg.\n',...
    abs(alphaL0vec(i) - alphaL0TATvec(i)));

end 