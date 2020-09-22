%This function will return the type, dx, dy, dz of a reaction joint as a
%cell array
function C = React(cell_input)
[A,B] = strtok(cell_input) ; %separates non-numerical values before converting to a number vector
D = str2num(B);
C = {A,D(1),D(2),D(3)};
%C is a 1x4 cell array that stores type F/m, dx, dy, dz
end