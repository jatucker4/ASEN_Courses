%Load the data from the ductile sample and save it to a variable
Brittle = load('Brittle2');
%Retrieve the force from the ductile data
force = Brittle(27:147 , 2);
%Retrieve the elongation from the ductile data
elongation = Brittle(27:147 , 3);
subplot(2,2,1);
scatter(elongation+0.000976,force-117.8);
xlabel("Elongation (in)");
ylabel("Force (lbf)");
%Calculating Stress and Strain including unit conversion from psi to Pa
Stress = (force/0.0995)*6894.75729;
Strain = elongation/1;
subplot(2,2,2);
scatter(Strain+0.000976,Stress-8.164e+06);
xlabel("Strain (in/in)");
ylabel("Stress (Pa)");
youngs = (5.567e+07 - 4.856e+07)/(0.0009911-0.0008813);
%This results in a youngs modulus of 64.75GPa therefore it's most likey
%aluminum 7000
x = linspace(0,.005,20);
y = youngs*x - youngs*.002;
subplot(2,2,3);
scatter(Strain + .000976,Stress-8.164e+06);
xlabel("Strain (in/in)");
ylabel("Stress (Pa)");
hold on;
plot(x,y);

slope = ((1.395e+08)-(1.385e+08))/(0.004205-0.004103);
intercept = (1.385e+08)-slope*(0.004103);

%Find x value where the two lines intersect
xIntersection = ((youngs*0.002)+ intercept)/(youngs-slope);
%From this we find where the yield strength
yieldStrength = youngs*xIntersection - youngs*0.002;

%Ultimate Tensile Strength
UTS =(max(force(:,1))/0.0995)*6894.75729;
%Fracture Stress is 149.27MPa
