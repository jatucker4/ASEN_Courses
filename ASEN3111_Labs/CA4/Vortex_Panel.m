function [Cl, CP] = Vortex_Panel(x,y,V_inf,alpha,plotcp,increment)
%Vortex_Panel Performs the calculations detailed in the Kuethe and Chow
%document
%
% Author: Johnathan Tucker
%
% Collaborators: N/A
% 
% Corrected using the reference solution posetd on canvas.
%
% This function takes in the x and y values from the NACA airfoil function,
% the free stream velocity, the angle of attack, plotcp flag, and plot
% increments. It outputs the coefficient of pressure, coefficient of lift,
% and displays a plot of coefficient of pressure vs x/c
%
% Last Revised: 3/26/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the chord length for calculating Cl
c = max(x) - min(x);

% Begin the translation from 4tran to MATLAB
% Create necessary variables
M = length(x) - 1;
MP1 = M + 1;
alpha = alpha *(pi/180);


for I = 1:M
    IP1 = I + 1;
    X(I) = 0.5*(x(I) + x(IP1));
    Y(I) = 0.5*(y(I) + y(IP1));
    S(I) = sqrt( (x(IP1) - x(I))^2 + (y(IP1) - y(I))^2 );
    theta(I) = atan2( (y(IP1) - y(I)) , (x(IP1) - x(I)));
    sine(I) = sin(theta(I));
    cosine(I) = cos(theta(I)); 
end
RHS = sin(theta - alpha);

for I = 1:M
    for J = 1:M     
        if I == J
            CN1(I,J) = -1;
            CN2(I,J) = 1;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
        else
            A = -(X(I) - x(J))*cosine(J) - (Y(I) - y(J))*sine(J);
            B = (X(I) - x(J))^2 + (Y(I) - y(J))^2;
            C = sin(theta(I) - theta(J));
            D = cos(theta(I) - theta(J));
            E = (X(I) - x(J))*sine(J) - (Y(I) - y(J))*cosine(J);
            F = log(1 + S(J)*(S(J) + 2.*A)/B);
            G = atan2(E*S(J), B + A*S(J));
            P = (X(I) - x(J))*sin(theta(I) - 2.*theta(J)) +...
                (Y(I) - y(J))*cos(theta(I) - 2.*theta(J));
            Q = (X(I) - x(J))*cos(theta(I) - 2.*theta(J)) -...
                (Y(I) - y(J))*sin(theta(I) - 2.*theta(J));
            CN2(I,J) = D + .5*Q*F/S(J) - (A*C + D*E)*G/S(J);
            CN1(I,J) = .5*D*F + C*G - CN2(I,J);
            CT2(I,J) = C + .5*P*F/S(J) + (A*D - C*E)*G/S(J);
            CT1(I,J) = .5*C*F - D*G - CT2(I,J);
        end
    end
end

for  I = 1:M
     AN(I,1) = CN1(I,1);
     AN(I,MP1) = CN2(I,M);
     AT(I,1) = CT1(I,1);
     AT(I,MP1) = CT2(I,M);
     for J = 2:M
         AN(I,J) = CN1(I,J) + CN2(I,J-1);
         AT(I,J) = CT1(I,J) + CT2(I,J-1);
     end
end

AN(MP1,1) = 1;
AN(MP1,MP1) = 1;
for J = 2:M
    AN(MP1,J) = 0;
end

RHS(MP1) = 0;

% Solve the system of equations using A\b instead of Cramers
GAMA = AN\RHS';

for I = 1:M
    V(I) = cos(theta(I) - alpha);
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
        CP(I) = 1 - V(I)^2;
    end
end
% Change from gamma prime to gamma via Kuethe and Chow
Gamma = 0;
for i = 1:M
    Gamma = Gamma + 2*pi*0.5*(GAMA(i) + GAMA(i+1))*S(i);
end
% Solve for CL by calculating capital GAMMA inline using CA3 Notes formulas
Cl = 2*Gamma;
%% Create the pressure plot
% This flag is if only one cp plot is wanted
if plotcp == 1
    figure
    cp_lower = CP(1:(length(x)+1)/2);
    cp_upper = CP((length(x)+1)/2:end);
    scatter(x((length(x)+1)/2:end-1)./c,cp_upper,'r')
    hold on
    scatter(x(1:(length(x)+1)/2)./c,cp_lower,'b')
    title('$Coefficeint\:of\:Pressure\:vs\:\\frac{x}{c}$','Interpreter','latex')
    xlabel('$x-distance\:[\% Chord]$','Interpreter','latex')
    ylabel('$Coefficient\:of\:Presure$','Interpreter','latex')
    legend('$Upper\:Surface$','$Lower\:Surface$','Interpreter','latex')
    
% This flag is for the cp subplots required for question two
elseif plotcp == 2
    hold on
    subplot(2,2,increment)
    cp_lower = -CP(1:(length(x)+1)/2);
    cp_upper = -CP((length(x)+1)/2:end);
    scatter(x(1:(length(x)+1)/2)./c,cp_lower,'b')
    hold on
    scatter(x((length(x)+1)/2:end-1)./c,cp_upper,'r')
    title(strcat('CP vs $\frac{x}{c}$ with: $\alpha$ = ',num2str(alpha*180/pi),'$^\circ$'),'Interpreter','latex');
    sgtitle('CP vs $\frac{x}{c}$ at Different $\alpha$ Values','Interpreter','latex');
    xlabel('$\frac{x}{c}$','Interpreter','latex')
    ylabel('$-CP$','Interpreter','latex')
    legend('$Lower\:Surface$','$Upper\:Surface$','Interpreter','latex')
end
    
end

