function [cl, cp] = Vortex_Panel(x,y,V_inf,alpha,doPlot)
%% Computational Lab 4 - Flow Over Thick Airfoils
% This function finds sectional coefficient of lift
%
% Author: Duncan McGough
% CUID: 104708558
% Date Created: 10/6/17


%% NOTE- PLEASE READ
% doPlot is a boolean input that tells the function whether or not to plot the airfoil. 
% If you are testing the function standalone, you are welcome to use either one of these choices:
%
% doPlot = true: will produce the Cp vs x location plots
% doPlot = false: won't produce any plots


%% Convert alpha to radians
alpha = alpha*pi/180;

%% Find chord length
c = x(1); 

%% FORTRAN Translation
m = length(x)-1;
mp1 = length(x);

for i=1:m
	xavg(i) = 0.5*(x(i) + x(i+1));
	yavg(i) = 0.5*(y(i) + y(i+1));
	s(i) = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
	theta(i) = atan2( (y(i+1)-y(i)), (x(i+1)-x(i)) );
end
sine = sin(theta);
cosine = cos(theta);
RHS = sin(theta-alpha);

for i=1:m
	for j=1:m
		if i == j
			CN1(i,j) = -1.0;
			CN2(i,j) = 1.0;
			CT1(i,j) = 0.5*pi;
			CT2(i,j) = 0.5*pi;
		else
			A = -(xavg(i)-x(j))*cosine(j) - (yavg(i)-y(j))*sine(j);
			B = (xavg(i)-x(j))^2 + (yavg(i) - y(j))^2;
			C = sin(theta(i)-theta(j));
			D = cos(theta(i)-theta(j));
			E = (xavg(i)-x(j))*sine(j) - (yavg(i)-y(j))*cosine(j);
			F = log( 1.0 + s(j)*(s(j)+2.*A)/B );
			G = atan2( E*s(j), B+A*s(j) );
			P = (xavg(i)-x(j)) * sin( theta(i)-2.*theta(j) ) + (yavg(i)-y(j)) * cos( theta(i)-2.*theta(j) ); 
			Q = (xavg(i)-x(j)) * cos( theta(i)-2.*theta(j) ) - (yavg(i)-y(j)) * sin( theta(i)-2.*theta(j) ); 
			CN2(i,j) = D + 0.5*Q*F/s(j) - (A*C+D*E)*G/s(j);
			CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
			CT2(i,j) = C + 0.5*P*F/s(j) + (A*D-C*E)*G/s(j);
			CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
		end
	end
end


for i=1:m
	AN(i,1) = CN1(i,1);
	AN(i,mp1) = CN2(i,m);
	AT(i,1) = CT1(i,1);
	AT(i, mp1) = CT2(i,m);

	for j=2:m
		AN(i,j) = CN1(i,j) + CN2(i,j-1);
		AT(i,j) = CT1(i,j) + CT2(i,j-1);
	end 
end

AN(mp1,1) = 1.0;
AN(mp1,mp1) = 1.0;

for j=2:m
	AN(mp1,j) = 0.0;
end

RHS(mp1) = 0.0;

% Can ignore the determinant and cramers rule code as MATLAB is great at
% matrix operations 
gama = AN\(RHS');

for i=1:m
	v(i) = cos( theta(i) - alpha );
	for j=1:mp1
		v(i) = v(i) + AT(i,j)*gama(j);
		cp(i) = 1.0 - v(i).^2;
	end
end

%% Calculate Cl
cl = sum(2.*v.*s/c);

if doPlot == true
	%% Plot Cp
	figure
	hold on
	grid on
	h = gca;
	set(h,'YDir','reverse');
	set(h,'FontSize',18);
	plot(xavg,cp', 'LineWidth',2)
	% plot(xavg,cp', '-o', 'LineWidth',2)
	title('Cp vs Location')
	xlabel('x')
	ylabel('Cp')
	hold off
end

end