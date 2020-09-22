%This function will return a unit vector in component form when v =
%[dx,dy,dx] are input
function unit = unitV(v)

%converts input v to a matrix if it is a cell
if iscell(v) == 1
v = cell2mat(v);
end 

%computes the unit vector by dividing v by the magnitude of v
unit = v/norm(v);

end