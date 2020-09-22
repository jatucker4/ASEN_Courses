function CL_3d = C_L_Lab1(alpha, CL, efficiency, AR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
alpha_0 = alpha(1,27) + (0 - CL(1,27))*(alpha(1,28)-alpha(1,27))/(CL(1,28)-CL(1,27));
x = alpha(1,48:56);
y = CL(1,48:56);
p = polyfit(x,y,1);
a_0 = p(1);
a = a_0/(1+ (57.3*a_0/(pi*efficiency*AR)));
CL_3d = a*(alpha - alpha(32));
end

