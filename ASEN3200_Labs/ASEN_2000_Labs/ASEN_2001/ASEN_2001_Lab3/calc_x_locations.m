function [x_locations,centroid,distances,tier_2] = calc_x_locations(p_naught)
% Declare constant
L = 0.9144;
% Symbolically declare x for integration
syms x  
% Create q equation
q = 4*.0254 * p_naught * sqrt(1-(2*x/L)^2)
% Splitting the bar up into 8 equal length;s (L/8)
vector = [-.4572, -.3429, -.2286,-.1143,0,.1143,.2286,.3429,.4572];
% Solve for the forces and centroid between each of the lengths in the
% vector
for i = 1:8
    force(i) = double(int(q,[vector(i),vector(i+1)]));
    centroid(i) = double(int(q*x,[vector(i),vector(i+1)]))/force(i);
end
for i = 1:3
    x_locations(i) = (centroid(i)-centroid(i+1))*(force(i+1)+.49+.5*1.77)/(force(i)+force(i+1)+2*.49+1.77);
    distances(i) = centroid(i)-centroid(i+1);
    
end
tier_2(1,1) = (x_locations(3)-x_locations(1))*(force(3)+force(4)+2*.49+1.77+.5*2.94)/((force(3)+force(4)+2*.49+1.77+2.94+force(1)+force(2)));

end