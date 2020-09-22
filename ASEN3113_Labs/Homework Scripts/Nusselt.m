function Nu = Nusselt(Ra,Pr)
%NUSSELT Summary of this function goes here
%   Detailed explanation goes here
Nu = (0.825 + (0.387*(Ra)^(1/6))/((1 + (.492/Pr)^(9/16))^(8/27)) )^2;
end

