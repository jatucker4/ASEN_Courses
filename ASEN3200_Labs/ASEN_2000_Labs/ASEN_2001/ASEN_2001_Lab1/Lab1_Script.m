%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Clear the workspace
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Open the file and create a file ID
fid = fopen('Lab1_Input.txt');
%Create an incrementor variable
i = 1;
row_number = 0;
%While not at the end of file
while~feof(fid)
    %Put each line from the input file into the cell array input_data
    input_data{i,1} = fgetl(fid);
    %Increment the index variable
    i= i+1;
end
%Close the file for best practice
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Process and clean up the data so that each variable contains certain data
%from the cell array
%Load the second row into a variable split by whitespace
external_forces_moments = strsplit(input_data{2,1});
%Get the number of forces from external_forces_moments variable
num_external_forces = str2double(external_forces_moments{1,1});
%Get the number of moments from external_forces_moments variable
num_external_moments = str2double(external_forces_moments{1,2});

%Get the coordinates for each external force
for j = 1:num_external_forces
    coordinates_matrix{j,:} = strsplit(input_data{(4+j),1});
end
force_coordinates_matrix = str2double(vertcat(coordinates_matrix{:}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keep track of what row we are at
row_number = 4+num_external_forces+3;
%Get magnitude and direction of external forces
i = 1;
for j = row_number:(row_number+num_external_forces-1)
    force_mag_direction_matrix{i,:} = strsplit(input_data{j,1});
    i = i+1;
end
force_mag_direction_matrix = str2double(vertcat(force_mag_direction_matrix{:}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update row number counter
row_number = row_number + num_external_forces + 2;
i = 1;
for j = row_number : (row_number + num_external_moments - 1)
    moment_coordinates_matrix{i,:} = strsplit(input_data{j,1});
    i = i+1;
end
moment_coordinates_matrix = str2double(vertcat(moment_coordinates_matrix{:}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update row number counter
row_number = row_number + num_external_moments +2;
i = 1;
for j = row_number : (row_number + num_external_moments - 1)
    moment_mag_direction_matrix{i,:} = strsplit(input_data{j,1});
    i = i+1;
end
moment_mag_direction_matrix = str2double(vertcat(moment_mag_direction_matrix{:}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_number = row_number + num_external_moments + 2;
i = row_number;
j = 1;
while input_data{i,1} ~= "# type (F/M) and direction of reaction "
    support_coordinates_matrix{j,:} = strsplit(input_data{i,1});
    i = 1 + i;
    j = j + 1;
end
support_coordinates_matrix = str2double(vertcat(support_coordinates_matrix{:}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the row counter variable
[rows,columns] = size(support_coordinates_matrix);
row_number = row_number + rows + 2;

i = 1;
for j = row_number:length(input_data)
    type_direction_of_reaction{i,:} = strsplit(input_data{j, 1});
    i = i+1;
end
type_direction_of_reaction = vertcat(type_direction_of_reaction{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the Reaction force Magnitudes
b = force_mag_direction_matrix(:,1);
force_mag_direction_matrix(:,1) = [];
b = [b;moment_mag_direction_matrix(:,1)];
moment_mag_direction_matrix(:,1) = [];

force_mag_direction_matrix(1,:) = force_mag_direction_matrix(1,:)/norm(force_mag_direction_matrix(1,:));
force_mag_direction_matrix(2,:) = force_mag_direction_matrix(2,:)/norm(force_mag_direction_matrix(2,:));

moment_mag_direction_matrix(1,:) = moment_mag_direction_matrix(1,:)/norm(moment_mag_direction_matrix(1,:));

force_x_components = force_mag_direction_matrix(:,1);
force_y_components = force_mag_direction_matrix(:,2);
force_z_components = force_mag_direction_matrix(:,3);
moment_x_components = moment_mag_direction_matrix(:,1);
moment_y_components = moment_mag_direction_matrix(:,2);
moment_z_components = moment_mag_direction_matrix(:,3);
m = b;
c(1,1) = force_x_components(1,1)*m(1,1);
c(2,1) = force_x_components(2,1)*m(2,1);
b(1,1) = sum(c);

c(1,1) = force_y_components(1,1)*m(1,1);
c(2,1) = force_y_components(2,1)*m(2,1);
b(2,1) = sum(c);

c(1,1) = force_z_components(1,1)*m(1,1);
c(2,1) = force_z_components(2,1)*m(2,1);
b(3,1) = sum(c);

b(4,1) = moment_x_components*m(3,1);
b(5,1) = moment_y_components*m(3,1);
b(6,1) = moment_z_components*m(3,1);

type_direction_of_reaction = vertcat(type_direction_of_reaction);
type_direction_of_reaction = str2double(type_direction_of_reaction);
type_direction_of_reaction(:,1) = [];

for j = 1:length(type_direction_of_reaction)
    type_direction_of_reaction(j,:) = type_direction_of_reaction(j,:)/norm(type_direction_of_reaction(j,:));
end

type_direction_of_reaction = transpose(type_direction_of_reaction);

