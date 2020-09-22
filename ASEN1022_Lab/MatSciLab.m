%Load the data from the ductile sample and save it to a variable
Ductile = load('Ductile');
%Retrieve the force from the ductile data
force = Ductile(: , 2);
%Retrieve the elongation from the ductile data
elongation = Ductile(: , 3);
%Plot the Force vs Elongation plot
subplot(2,2,1);
scatter(elongation + 0.0050 , force);
xlabel("Elongation (in)");
ylabel("Force (lbf)");
%Calculating Stress and Strain including unit conversion from psi to Pa
Stress = (force/.0995)*6894.75729;
Strain = elongation/1;
subplot(2,2,2);
scatter(Strain + .005,Stress);
xlabel("Strain (in/in)");
ylabel("Stress (Pa)");
youngs = (Stress(13,1)./(Strain(13,1)+0.005));
%This results in a youngs modulus of 68.7GPa therefore it's Al 6061-T6
%because the actual modulus of elasticity for that sample is 68.9GPa
x = linspace(0,.007,20);
y = youngs*x - youngs*.002;
subplot(2,2,3);
scatter(Strain + .005,Stress);
xlabel("Strain (in/in)");
ylabel("Stress (Pa)");
hold on;
plot(x,y);

slope = ((2.529e+08)-(2.48e+08))/(0.00708-0.0052);
intercept = (2.48e+08)-slope*(0.0052);

%Find x value where the two lines intersect
xIntersection = ((youngs*0.002)+ intercept)/(youngs-slope);
%From this we find where the yield strength
yieldStrength = youngs*xIntersection - youngs*0.002;

%Ultimate Tensile Strength
UTS =(max(force(1:24))/0.0995)*6894.75729;
%Fracture Stress is 277.7MPa





