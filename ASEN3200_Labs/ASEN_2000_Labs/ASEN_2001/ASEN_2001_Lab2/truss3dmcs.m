function truss3dmcs(inputfile)
%
% Stochastic analysis of 2-D statically determinate truss by
% Monte Carlo Simulation. Only positions and strength of joints 
% treated as random variables
%
% Assumption: variation of joint strength and positions described 
%             via Gaussian distributions
% 
%             joint strength : mean = 4.8
%                              coefficient of variation = 0.4
%             joint position : 
%                              coefficient of variation = 0.01
%                              (defined wrt to maximum dimension of truss)
%
%             number of samples is set to 1e5
%
% Input:  inputfile  - name of input file
%
% Author: Kurt Maute for ASEN 2001, Oct 13 2012

% parameters
% mean of joint strength 4.8 N
jstrmean   = 4.8;
% coefficient of variation (sigma/u) of joint strength = 0.4/4.8 N
jstrcov    = 0.08;
% coefficient of variation of joint position percent of length of truss(ext)
jposcov    = 0.01;
% number of samples
numsamples = 1e5;   
desired = input('What is the desired probability of failure?\n');
% Calculate the factor of safety given a chosen desired probability for the
% design
FOS = factor_of_safety(desired);

% read input file
[joints,connectivity,reacjoints,reacvecs,...
    loadjoints,loadvecs]=readinput(inputfile);

% determine extension of truss
ext_x=max(joints(:,1))-min(joints(:,1));   % extension in x-direction
ext_y=max(joints(:,2))-min(joints(:,2));   % extension in y-direction
ext_z=max(joints(:,3))-min(joints(:,3));
ext=max([ext_x,ext_y,ext_z]);

% loop overall samples
numjoints=size(joints,1);       % number of joints
maxforces=zeros(numsamples,1);  % maximum bar forces for all samples
maxreact=zeros(numsamples,1);   % maximum support reactions for all samples
failure=zeros(numsamples,1);    % failure of truss

for is=1:numsamples 
    
    % generate random joint strength limit
    varstrength = (jstrcov*jstrmean)*randn(1,1);
    
    jstrength = jstrmean + varstrength;
    jstrength = jstrength/1.44;
    % generate random samples
    varjoints = (jposcov*ext)*randn(numjoints,3);
    
    % perturb joint positions
    randjoints = joints + varjoints;
    
    % compute forces in bars and reactions
    [barforces,reacforces] = forceanalysis(randjoints,...
        connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
    
    % determine maximum force magnitude in bars and supports
    maxforces(is) = max(abs(barforces));
    maxreact(is)  = max(abs(reacforces));
    
    % determine whether truss failed
    failure(is) = maxforces(is) > jstrength || maxreact(is) > jstrength;
end

figure(1);
subplot(1,2,1);
histogram(maxforces,30);
title('Histogram of maximum bar forces');
xlabel('Magnitude of bar forces');
ylabel('Frequency');

subplot(1,2,2);
histogram(maxreact,30);
title('Histogram of maximum support reactions');
xlabel('Magnitude of reaction forces');
ylabel('Frequency');

fprintf('\nFailure probability : %e With a factor of safety: %f\n\n',...
    sum(failure)/numsamples,FOS);
end
