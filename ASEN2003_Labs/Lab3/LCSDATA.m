%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to process the labview data by loading it in,
% subtract an integer number of full cycles , outputs the experimental
% angle, angular rate, and vertical velocity of the first 6 revolutions
%
% Created by: Lexie Marinelli and Johnathan Tucker
%
% Inputs: 
%           filename = name of the data file to be processed
%
%
% Outputs:
%           rad_exp = the experimental anglular position of the first 6
%                       revolutions [degrees]
%           w_exp = the experimental angular velocity of the first 6
%                   revolutions [degrees/s]
%           v_exp = the experimental vertical velocity of the first 6
%                   revolutions[cm/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rad_exp, w_exp, v_exp] = LCSDATA(filename)

% Open the file
file = fopen(filename);
% Get a cell array delimited by spaces
A1 = textscan(file, '%s', 'Delimiter', '\n');
% Get the first column of the cell array
A1 = A1{:};
% Parse the column using a regular expression to split it at spaces
A1 = regexp(A1, '\s+', 'split');
% Certically concatinate and save it to a cell array
A2 = vertcat(A1{:});
% Pre allocate final matrix
A = zeros(size(A1,1), 6);
% Loop through the cell array and convert it to a matrix
for i = 1:size(A1,1)
    A(i,1) = str2num(cell2mat(A2(i,1)));
    A(i,2) = str2num(cell2mat(A2(i,2)));
    A(i,3) = str2num(cell2mat(A2(i,3)));
    A(i,4) = str2num(cell2mat(A2(i,4)));
    A(i,5) = str2num(cell2mat(A2(i,5)));
    A(i,6) = str2num(cell2mat(A2(i,6)));

end
% Close the file
fclose(file);
% Save the matrix to variable lcs
lcs = A; 
%% Process the data
% Take the modulus of the number in the slide position column with the 
n = 0;
%columns: index   revolution_number   angular_pos   
for i = 1:length(lcs)
    theta_pos(i,1) = i;
    theta_pos(i,2) = mod(lcs(i,2),360);
end
% Create a counter variable that keeps track of the revolution number
counter = 0;
% Loop through the angular position matrix and add the rev number times 360
% to ensure the angular position is coninuous
for i = 1:length(lcs)-1
    m = (theta_pos(i,2) - theta_pos(i+1,2));
    if m > 100 
        indices(i,1) = i;
        counter = counter + 1;
    end
    theta_pos(i,2) = theta_pos(i,2) + 360*counter;
end
% Find the indices where the rotation is ended "at 360"
indices = find(indices);
% Create a vector of the angle values
b1 = theta_pos(indices,2);
% Subtract 360 to make the transition smooth between revolutions
b = b1 - 360;
% Overwrite the values in the angular position matrix
theta_pos(indices,2) = b;
% Output the required values
rad_exp = theta_pos(indices(1):indices(6),2);
w_exp = lcs(indices(1):indices(6),4);
v_exp = lcs(indices(1):indices(6),5) / 10;

end

