function output_mat = sts_dataprocess(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was created to process the static test stand data in a
% general way that could be applied to any test file.
%
% Written by: Johnathan Tucker
%
% Inputs: 
%       filename: The name of the file
%
% Outputs:
%       output_mat: A cell array of the water mass used, pressure,
%       temperature, group, time of trial, thrust vector, ISP value, and
%       time of thrust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file and create a fileID object for closing
fileID = fopen(filename);
% Save everything on the file to a character array for language processing
char_data_array = fscanf(fileID,'%c');
% Conversion for sample rate
sample_rate = 1.652 * 1000;

% Get the amount of water tested
% Locate the string in the data file: NOTE strfind returns the index of the
% first letter in the search string hence the plus and minus 
mass_index_start = strfind(char_data_array,'Water')+7;
mass_index_end = strfind(char_data_array,'grams')-1;
% Convert the string to a double 
water_mass = str2double(string(char_data_array(mass_index_start:mass_index_end)));
% Be sure to convert units
water_mass = water_mass/1000;

% Get the pressure of the bottle rocket
% Locate the string in the data file: NOTE strfind returns the index of the
% first letter in the search string hence the plus and minus 
press_index_start = strfind(char_data_array,'Pressure')+10;
press_index_end = strfind(char_data_array,'(psi)')-1;
% Convert the string to a double 
press = str2double(string(char_data_array(press_index_start:press_index_end)));
% Be sure to convert units
press = press * 6894.76;

% Get the temperature
% Locate the string in the data file: NOTE strfind returns the index of the
% first letter in the search string hence the plus and minus 
temp_index_start = strfind(char_data_array,'Temperature')+13;
temp_index_end = strfind(char_data_array,'(C)')-1;
% Convert the string to a double
temp = str2double(string(char_data_array(temp_index_start:temp_index_end)));
% Be sure to convert units
temp = temp + 273.15;

% Get the group number
% Locate the string in the data file: NOTE strfind returns the index of the
% first letter in the search string hence the plus and minus 
group_index = strfind(char_data_array,'%4')-2;
% Convert the string to a double
group = str2double(string(char_data_array(group_index)));

% Get the time
% Locate the string in the data file: NOTE strfind returns the index of the
% first letter in the search string hence the plus and minus 
time_index_start = strfind(char_data_array,'2019')+4;
time_index_end = strfind(char_data_array,'AM')-2;
% Leave the time as a string
time = string(char_data_array(time_index_start:time_index_end));

% Be sure to close the file so we don't crash matlab
fclose(fileID);

% Get the thrust data from the file
% NOTE: Load pulls only the columns of data from the file and handles
% opening/closing itself
% I went with this because trying to search for numbers in a char array
% sounded like a slow and painful death
thrust_data_array = load(filename);
% Convert to Newtons
thrust = thrust_data_array(:,3)* 4.44822;
% Create time vector for the plot
time_vec = linspace(0,length(thrust),...
    length(thrust))/sample_rate;
% Plot the thrust vs time
plot(time_vec,thrust);
% Get chosen bounds
[x,y] = ginput(2);
clf
% Find smallest absolute distance between chosen points and thrust vector
[~,y_min_1] = min(abs(time_vec-x(1)));
[~,y_min_2] = min(abs(time_vec-x(2)));
% Create the thrust vector from the chosen bounds
thrust = thrust(y_min_1:y_min_2);
% Create a time vector
time_vec = linspace(0,length(thrust),length(thrust))/sample_rate;
% This was out of laziness so that I wouldn't have to change the name of
% one variable
thrust_int = thrust;

% Calculate ISP
ISP = trapz(time_vec,abs(thrust_int))/(water_mass*9.81);

% Output everything
output_mat = {water_mass,press,temp,group,time,...
    thrust,ISP,time_vec};

end

