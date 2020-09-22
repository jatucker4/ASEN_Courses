function truss3d(input_filename)
%Primary Developer: Emma Markovich
%emma.markovich@colorado.edu
%Team Members: Johnathan Tucker, Brian Jackman
%johnathan.tucker@colorado.edu, brian.jackman@colorado.edu
%ASEN 2001 - Section 012
%Lab 1 - Computer Analysis of Structures

%% Purpose
%This code will input a static equilibrium problem filename and return the reaction forces of each support of a structure when 
%the forces/moments and their positions are specified by the user

%% (1) Set Up - Input Routine
%read in the static structure file and parse each new line as a cell in
%cell array C
filename = input_filename;
C = read_input(filename);

%define size of cell array that stores all data for reference in the rest of the script
[rmax, ~] = size(C);
%% (2) Set Up - Extracting and Organizing Data from Input File
% organize data into appropriate vectors/matrices for info given

%number of external forces and moments
[num_ext_forces,num_ext_moments] = num_ext(C{2,1});

%location of external forces
F_coordinates = zeros(num_ext_forces,3);%initialize F_coordinates as a matrix to store the coordinates of external forces

r = 5; %the row number of the cell with appropriate force data
 
for i = 1:num_ext_forces
    F_coordinates(i,:) = Coordinates(C{r,1});
    r = r+1; %iterating row counter
end

%magnitude and direction of forces
r = r+2; %defines the row first where this data is stored
%initialize matrix to store force mag and direction data
F_magnitude_direction = zeros(num_ext_forces,4);

for i = 1:num_ext_forces
    F_magnitude_direction(i,:) = Mag_Direc(C{r,1});
    r = r+1; %iterating row counter
end

%location of external moments
M_coordinates = zeros(num_ext_moments,3);%initialize F_coordinates as a matrix to store the coordinates of external forces

r = r+2; %define row that appropriate data begins on
for i = 1:num_ext_moments
    M_coordinates(i,:) = Coordinates(C{r,1});
    r = r+1; %iterating row counter
end

%Magnitude and Direction Moments
r = r+2; %defines the row first where this data is stored
%initialize matrix to store force mag and direction data
M_magnitude_direction = zeros(num_ext_moments,4);

for i = 1:num_ext_moments
    M_magnitude_direction(i,:) = Mag_Direc(C{r,1});
    r = r+1; %iterating row counter
end

%Location of Supports
r = r+2; %define row that appropriate data begins on
flag = 0;
i = 1;
while flag == 0  %define flag that will stop the loop when there are no more supports
    S_coordinates(i,:) = Coordinates(C{r,1});
    r = r+1; %iterating row counter
    i = i+1;  %iterate loop counter variable
    if contains(C{r,1},'#')
        flag = 1;
    end
end

%type and direction of reaction
%initializing counter/flag variables
r = r+2;
flag = 0;
i = 1;
while flag == 0 %define flag that will break the loop when there are no more lines of data
    reaction(i,:) = React(C{r,1});
    %iterating counting variables
    i = i+1;
    r=r+1;
    if r > rmax
        flag = 1;
    end
end

%clear unnecessary counter variables
clear i, clear r, clear flag;

%% (3) Computational Routine - Creating A, x, and b

%find/define size of necessary matrices
[rF,cF] = size(F_magnitude_direction);
[rM,cM] = size(M_magnitude_direction);

num_Fsupports = sum(count(reaction(:,1),'F'));
num_Msupports = sum(count(reaction(:,1),'M'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A - coefficient matrix to store presence of reaction forces (should be 6 x #supports)
%initializing A matrix to be correct size
A = zeros(6,num_Fsupports + num_Msupports);

    %convert reaction directions to unit vectors
    [rA,cA] = size(A);
    reaction_unit = zeros(cA,3);
    
    %normalize reaction directions to be unit vectors
    for i = 1:cA
        reaction_unit(i,:) = unitV(reaction(i,2:4));
    end
    
    %create 2 matrices that store sum Forces and sum Momenets (plus needed zero
    %padding)
    A_F = [reaction_unit(1:num_Fsupports,:)',zeros(3,num_Msupports)]; %sum Forces + zero padding
    
    %determining true reaction Mx,My,Mz by crossing radius (distance vector from point to respective axis)
    %with unit vector of moment at that support point
    m_fromF = zeros(num_Fsupports,3); %preallocating m_fromF as temporary variable for loop to store reaction moments cause by reaction forces

    for i = 1:num_Fsupports
        r = S_coordinates(i,:); %r is the radius from the origin
        u = reaction_unit(i,:); %u is unit reaction force
        M = cross(r,u); %M is the moments created by the unit vector of force
        m_fromF(i,:) = M;
    end
    
    m_fromM = reaction_unit(num_Fsupports+1:length(reaction_unit),:)'; %moments cause by reaction moments, calculated by projection of M onto unit vector axis
    A_M = [m_fromF',m_fromM]; % Moments due to reaction forces and moments
    
    %vertically concatenate F and M
    A = [A_F;A_M];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%x - column vector of magnitude of reaction forces and moments (should be #supports x 1)
x = zeros(length(reaction),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%b - column vector with elements that are the inverse of the totals for Fx, Fy,
%Fz, Mx, My, Mz from external forces and moments (should be 6x1)

    %convert external force/moment directions to unit vectors
    F_unit = zeros(rF,4);
    M_unit = zeros(rM,4);
    for i = 1:rF
        F_unit(i,1) = F_magnitude_direction(i,1);
        F_unit(i,2:4) = unitV(F_magnitude_direction(i,2:4));
    end
    
    for j = 1:rM
        M_unit(j,1) = M_magnitude_direction(j,1);
        M_unit(j,2:4) = unitV(M_magnitude_direction(j,2:4));
    end
    
    %user defined function create_b sums and negates
    %the external Fx, Fy, Fz, Mx, My, Mz strictly due to given components
    b = create_b(F_unit,M_unit,F_coordinates);
%% (4) Computational Routine - Solving Linear System of Equations for x
%A\b is a built in function that uses one of a plethora of options to solve
%the linear system, whether A is square or not

x = A\b;

%% (5) Output Routine - Write Results Cell Array
%use reaction directionaly unit vectors to print reaction components
results = cell(length(x)+1,4);

for j = 1:num_Fsupports
    RForce = x(j)*reaction_unit(j,:); %multiply magnitude by unit vector to find true component forces
    results(j,1:4) = {['F',num2str(j)],RForce(1),RForce(2),RForce(3)};      
end

for i = 1:num_Msupports %count number listing of the reaction moment
    j = j+1; %keeping count of element of x vector that we are pulling data from
    RMoment = x(j)*reaction_unit(j,:); %multiply magnitude by unit vector to find true component forces
    results(j,1:4) = {['M',num2str(i)],RMoment(1),RMoment(2),RMoment(3)};
end

%% (6) Output Routine - Assigning Names to Support Reactions
%pop-up a menu for the user to enter 
choice = menu('Would you like to enter reaction force/moment names?','Yes','No');

%Prompt user to enter names for the reaction forces
if choice == 1  %if user chooses "yes"
    prompt = {'Enter comma-seperated reaction force/moment name (in order of input):'};
    title = 'Reaction Force Labels';
    dimensions = [1 73];
    example_input = {'Ax,Ay,Az,Bx,Bz,Cz'};
    answer = inputdlg(prompt,title,dimensions,example_input); %bring up dialogue box for user input
    %convert answer from cell array to string
    answer = string(answer);
    %separate user-input string into 6 labels for forces
    Reaction_Labels = strsplit(answer,',');
    
    for i = 1:6
        results(i,1) = {Reaction_Labels(i)};
    end
end

%% (7) Output Routine - Save Results to ACSII file
[rResults,cResults] = size(results);
%define filename and open that file for writing - note: will write over anything already inside of file
filename_results = strcat('ReactionResults_',filename); %different for each different input file
fileResults = fopen(filename_results,'w');

%print input file name on top of results 
fprintf('Input File: %s\r\n',filename);
%begin text with column headers
fprintf(fileResults,'%1s %7s %8s %7s \r\n', 'Reaction','x','y','z');
%use for loop to print data in columns aligning with the headers
for i = 1:num_Fsupports
    fprintf(fileResults,'%3s %14.3f %8.3f %8.3f\r\n',results{i,1},results{i,2:cResults});
end
for i = 1:num_Msupports
    fprintf(fileResults,'%3s %14.3f %8.3f %8.3f\r\n',results{i+num_Fsupports,1},results{i+num_Fsupports,2:cResults});
end
%view file in command window
type(filename_results);
%close file
fclose(fileResults);

end