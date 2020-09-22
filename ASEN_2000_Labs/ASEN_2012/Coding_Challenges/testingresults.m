%% Housekeeping: Clear the variables and close all open plots
clear all;
close all;
clc;
%
%% Initialization: Define the Functions, Load Data, and Set Known Values
% 
% 1. The function f(x) is specified in the file func.txt. Open the file
% using fopen() to read the function handle.
filename = 'func.txt';  
fh = fopen(filename);       % Open the file, return a handle
% These two lines create a function handle from the contents of the file.
% Setting up the function this way is atypical, but it could be useful 
% to test many arbitrary functions in a single file.
string_exp = fgetl(fh);     % Read the single line
f_x = eval(string_exp);     % Evaluate the single line to get a handle

% 2. Specify an output file to save results to
output_filename = 'bisection-results.txt';
output_file = fopen(output_filename, 'w');

% 3. Set the input variables
% Interval a,b; convergence tolerance; maximum number of iterations
a = 0.0;
b = 6.0;
tol = 10^-6;
max_iter = 1000;

%% Begin The Algorithm
% 4. Begin iterative procedure.
for i = 1:1:max_iter

    % 5. Calculate the mid-point, p
    p = (a+b)/2;
    
    % 6. If f(p) has the same sign as f(b) set b=p, otherwise, set a=p
    if f_x(p)*f_x(b) > 0 
        b = p;
    else
        a = p;
    end

    
    % 7. If the absolute value of (b-a) is less than the convergence 
    % tolerance, print the solution and end the program
    if abs(a-b) <= tol
        fprintf(output_file, 'The converged root of the function is %f\n',p);
        fprintf(output_file, 'The function value at the converged root is %f\n',f_x(p));
        break
    end   

end

% 8. If the number of iterations is equal to the maximum number of iterations,
% print maximum number of iterations exceeded? and end the program
if i==max_iter
    fprintf(output_file, 'Maximum number of iterations, %d, reached\n',i)
end