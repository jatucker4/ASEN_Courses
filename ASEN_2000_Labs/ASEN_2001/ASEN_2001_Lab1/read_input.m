%This function will read in the ASCII file and return a cell array with each
%row as a new line in the .txt input file
function C = read_input(filename)

%open file
file = fopen(filename);

%initilize variable to count number of lines
i = 1;
while ~feof(file) %while the end of the file has not been reached, run the inside
    C{i,1} = fgetl(file);
    i = i+1; %increment i to input the proceeding line read into the row below
end

end