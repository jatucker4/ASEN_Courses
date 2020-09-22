function [x, y] =  NACA_Airfoil(string,c,N)
%% Author: Duncan McGough
% CUID: 104708558
% Date Created: 10/6/17

%% Convert the NACA string into relevant information - Lucas Calvert
m = str2num(string(1))/100; % maximum camber
p = str2num(string(2))/10; % location of maximum camber
t = str2num(string(3:4))/100;



%% Set up x vector
x = linspace(0,c,N);
yt = (t/0.2)*c*(0.2969*sqrt(x./c) - 0.1260*(x./c) - 0.3516*(x./c).^2 + 0.2843*(x./c).^3 - 0.1036*(x./c).^4); % NACA airfoil equation for 12% thickness to chord


%% Find the location of pc
pc_exact = p*c; % exact location
[min_diff, pc_index] = min(abs(pc_exact-x)); % find the closest value to the exact value in the x array


%% Find the mean camber line 
if pc_index ~= 1
	for i=1:pc_index
		yc(i) = m*x(i)/p^2 * (2*p - x(i)/c);
	end

	for i=(pc_index+1):length(x)
		yc(i) = m*(c - x(i))/(1 - p)^2 * (1 + x(i)/c - 2*p);
	end
else
	for i=1:length(x)
		yc(i) = m*(c - x(i))/(1 - p)^2 * (1 + x(i)/c - 2*p);
	end
end


%% Define zeta
Zeta = atan2(diff(yc), diff(x));
Zeta(N) = 0;


%% Find upper and lower x and y
xu = x - yt.*sin(Zeta);
xl = x + yt.*sin(Zeta);
yu = yc + yt.*cos(Zeta);
yl = yc - yt.*cos(Zeta);



%% Flip some vectors so we move from trailing edge all the way around the airfoil back to the airfoil trailing edge
xu = fliplr(xu);
yu = fliplr(yu);


%% Find x and y, remove extra zero in middle
x = [xu xl(2:end)];
y = -[yu yl(2:end)];

%% Plot airfoil for checking
% plot(x,y)
%axis equal

end