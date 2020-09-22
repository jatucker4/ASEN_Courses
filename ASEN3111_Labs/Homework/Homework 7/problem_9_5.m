M1 = 3.2;
theta = 18.2*pi/180;

beta1 = get_beta(theta,M1);

Mn1 = M1*sin(beta1);

Mn2 = get_Mn2(Mn1);

M2 = Mn2/sin(beta1 - theta);

beta2 = get_beta(theta,M2);

Mn2_2 = M2*sin(beta2);

Mn3 = get_Mn2(Mn2_2);

M3 = Mn3/sin(beta2 - theta);

rho2orho1 = get_rho2orho1(Mn1);
p2op1 = get_p2op1(Mn1);
T2oT1 = get_T2oT1(rho2orho1,p2op1);

rho3orho2 = get_rho2orho1(Mn2_2);
p3op2 = get_p2op1(Mn2_2);
T3oT1 = get_T2oT1(rho3orho2,p3op2);
