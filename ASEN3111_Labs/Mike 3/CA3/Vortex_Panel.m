function c_l = Vortex_Panel(x,y,V_inf,alpha,CpSwitch,Digits,PlotLoop)
%Vortex_Panel: This function completes all of the tasks required in
%question 1. 

%This function executes the vortex panel method as detailed in the in
%provided Kuethe and Chow (KC) Fortran code. The result of this method 
%is then used to compute and plot the coefficient of pressure over the 
%upper and lower surfaces for the airfoil being evaluated. 

% Author: Michael Martinson
% Collaborators: Isaac Goldner, Trace Valade, Ruben Hinojosa Torres
% Date: 3/27/2020

%% Notes
%alpha entered in rad. 

%% In order to better translate Kuethe and Chow, convert variables to the 
%same naming conventions they use: 

%Panel count: 
N = length(x) - 1;
M = N; 
MP1 = M + 1; 

%Boundary points:
XB = x; %m 
YB = y; %m

%% Implement the first "Do" loop in the KC code:
for I = 1:M
    
    IP1 = I + 1;
    X(I) = 0.5*(XB(I) + XB(IP1)); 
    Y(I) = 0.5*(YB(I) + YB(IP1));
    S(I) = sqrt(((XB(IP1) - XB(I))^2) + ((YB(IP1) - YB(I))^2)); 
 
    THETA(I) = atan2( (YB(IP1) - YB(I)),  (XB(IP1) - XB(I))); 
    SINE(I) = sin(THETA(I)); 
    COSINE(I) = cos(THETA(I)); 
    
    %Matrix equation right hand side:
    RHS(I) = sin(THETA(I) - alpha);  
end

%% Nested "Do" loops in the KC Code: 
for I = 1:M
    for J = 1:M
        
        if I == J
            CN1(I,J) = -1.0; 
            CN2(I,J) = 1.0; 
            
            CT1(I,J) = 0.5*pi; 
            CT2(I,J) = 0.5*pi; 

        else
            A = -(X(I) - XB(J))*COSINE(J) - (Y(I) - YB(J))*SINE(J); 
            B = ((X(I) - XB(J))^2) + ((Y(I) - YB(J))^2);
            C = sin(THETA(I) - THETA(J)); 
            D = cos(THETA(I) - THETA(J)); 
            E = (X(I) - XB(J))*SINE(J) - (Y(I) - YB(J))*COSINE(J); 
            F = log( 1.0 + (S(J)*(S(J) + 2.*A)/(B)) ); 
            G = atan2( E*S(J), B + A*S(J) ); 
            
            P = (X(I) - XB(J))*sin(THETA(I) - 2*THETA(J))...
                + (Y(I) - YB(J))*cos(THETA(I) - 2*THETA(J)); 
            
            Q = (X(I) - XB(J))*cos(THETA(I) - 2*THETA(J))...
                - (Y(I) - YB(J))*sin(THETA(I) - 2*THETA(J)); 
            
            CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C + D*E)*G/S(J); 
            CN1(I,J) = 0.5*D*F + C*G - CN2(I,J); 
            
            CT2(I,J) = C + 0.5*P*F/S(J) + (A*D - C*E)*G/S(J); 
            CT1(I,J) = 0.5*C*F - D*G - CT2(I,J); 
            
        end 
        
    end 
    
end 

%% Compute influence coefficients in eqns 5.47 and 5.49, respectivly 
%(equations in the KC code documentation)

for I = 1:M
  AN(I,1) = CN1(I,1); 
  AN(I,MP1) = CN2(I,M); 
  AT(I,1) = CT1(I,1); 
  AT(I,MP1) = CT2(I,M); 
  
  for J = 2:M
      AN(I,J) = CN1(I,J) + CN2(I,J-1); 
      AT(I,J) = CT1(I,J) + CT2(I,J-1); 
  end 
  
end 

AN(MP1,1) = 1.0; 
AN(MP1,MP1) = 1.0; 

for J = 2:M
   AN(MP1,J) = 0; 
end

RHS(MP1) = 0; 

%% Checking work (commented out intentionally)
% CheckMat = horzcat(X',Y',THETA',S'); 

%% Solving eq 5.47 for dimensionless strengths Gamma: 
%Solve for Gamma' = GAMA (KC code convention)
GAMA = AN\RHS'; 

%% Solve for dimensionless Cp (labeled CP) and dimensionless velocity
%Vi (Labled V): 
%Use equations provided in CA3 lab notes for Cp and Vi: 
for I = 1:M 
    V(I) = cos(THETA(I) - alpha); 
    
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
        CP(I) = 1 - (V(I)^2) ; 
    end 
end

%Sort results into column vectors, rename results: 
Vi = V'; %unitless
Cp = CP'; %unitless
gammaPrime = GAMA; 

%% find cl for the desired inputs: 
%Use equation process outlined in CA3 in lab notes:

%find c, airfoil chord:
c = max(x); %m 

%Compute gamma: 
gamma = (gammaPrime)*(2*pi)*(V_inf); 

%Compute Gamma (circulation):
%Assign S as a column vector:
Scol = S'; 

%For j = 1:M, find gamma(j)*S(j) (CA3 Class Notes Equation)
SumInner = gamma(1:M,1).*Scol(1:M,1); 

%Compute the sum from 1 to M: 
SumRes = sum(SumInner); 

%Compute Gamma (circulation): 
Gamma =  SumRes; 

%calculate sectional lift coeff cl for the airfoil: 
c_l = ((2*Gamma)/(V_inf*c)); %unitless
%This value is returned from the function. 

%% Plotting: Cp vs. x/c 
%Assign the leading zeros to the NACA number if it is the NACA0012 airfoil:
if Digits == 12
    AirfoilName = ['00',num2str(Digits)]; 
else 
    AirfoilName = num2str(Digits); 
end 


% %(If Cp plotting switch is set for single plot) (Commented out, only used
% %for testing)
% if CpSwitch == 1
% 
% %Plot and label: 
% figure; 
% scatter(X/c,Cp,8,'filled'); %plotting both surfaces together (this config
% %is only used for testing.
% 
% set(gca, 'YDir','reverse'); 
% grid on; grid minor; 
% xlabel('$\frac{x}{c}$','FontSize',14)
% ylabel('$C_p$','FontSize',14)
% 
% %Title: set to include the airfoil, alpha, and N: 
% title(['$C_p$ vs. $\frac{x}{c}$: ',...
%     '$\alpha$ = ',num2str(alpha*(180/pi)),'$^{\circ}$',...
%     ', N = ',num2str(M),...
%     ', NACA',AirfoilName],...
%     'FontSize',14)

%If CpSwitch is on the problem 2 setting to subplot the figures for 
%multiple angles of attack:
if CpSwitch == 2

%Plot and label: 
hold on; 
subplot(2,2,PlotLoop)

%Assign X/c vector
Xoverc = X./c; %unitless

%Plot top and bottom surfaces same color (testing only): 
%scatter(X/c,Cp,8,'filled'); %Plotting both surfaces together (Only used in
%testing) 

%Plot the Cp values vs. X/c for the lower surface: 
scatter(Xoverc(1:(length(Xoverc)/2)),Cp(1:(length(Cp)/2)),...
    8,'filled','b');

%Plot the Cp values vs. X/c for the upper surface: 
hold on;
scatter(Xoverc((length(Xoverc)/2)+1:end),Cp((length(Cp)/2)+1:end),...
    8,'filled','r');

%Flip the y axis:
set(gca, 'YDir','reverse'); 

%Label:
grid on; grid minor; 
xlabel('$\frac{x}{c}$','FontSize',14)
ylabel('$C_p$','FontSize',14)

%Title: set to include the airfoil, alpha, and N: 
title(['$C_p$ vs. $\frac{x}{c}$: ',...
    '$\alpha$ = ',num2str(alpha*(180/pi)),'$^{\circ}$',...
    ', N = ',num2str(M),...
    ', NACA',AirfoilName],...
    'FontSize',14) 

%Mark the surfaces: 
legend('Lower Surf.','Upper Surf.')

%else if plotting is not desired: 
else
    %Do not plot Cp vs x/c
    
end 
