function [F,M] =  num_ext(cell_input)
%This function will read the line containing number of external forces and
%moments, split the line, and convert back to double

[F,M] = strtok(cell_input);
F = str2double(F); %convert from string to numeric
M = str2double(M);
end