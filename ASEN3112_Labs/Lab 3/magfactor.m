function [t, mag] = magfactor(tdat, d0, dx)
%% Function magfactor
%
% Purpose - to calculate the magnification factor as a function of time for
% datasets
%
% Method - there are very distinct points in the data set where the
% magnification factor can be taken for each frequency, and this is done by
% finding the local minima and maxima of the sample displacement dx divided
% by the input (shaker) displacement d0, and finding the value of d0/dx in
% the middle of these peaks. The peaks correlate strongly to a 0 dx value,
% so ideally the value of d0/dx in between peak indices should correspond
% to the magnification factor for the given cycle. The peaks are
% implemented using findpeaks, with local minima found by using findpeaks
% for negative data
%
% Inputs:
% tdat  - time data 
% d0    - input displacement data
% dx    - test displacement data
%
% Outputs:
% t     - time data corresponding to the magnification factor values
% mag   - magnification factor

d = dx./d0;

[~, imax] = findpeaks(d);
[~, imin] = findpeaks(-d);

% Trim peak data

limax = length(imax);
limin = length(imin);

if (limax > limin)
    imax = imax(1:limin);
else 
    imin = imin(1:limax);
end

% Fold and sort data

ipeaks = sort([imax; imin]);

imean = 0.5*(ipeaks(1:end-1) + ipeaks(2:end));

t = interp1(1:length(tdat), tdat, imean);
mag = interp1(1:length(d), d, imean);

[mag, rmi] = rmoutliers(mag, 'movmedian', 3);
t = t(~rmi);

[mag, rmi] = rmoutliers(mag, 'movmedian', 5);
t = t(~rmi);
% 
% [mag, rmi] = rmoutliers(mag, 'movmedian', 5);
% t = t(~rmi);
% 
% [mag, rmi] = rmoutliers(mag, 'movmedian', 5);
% t = t(~rmi);

end