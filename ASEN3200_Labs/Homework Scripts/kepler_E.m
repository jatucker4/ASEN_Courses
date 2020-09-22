% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
function E = kepler_E(e, M)
% wwwwwwwwwwwwwwwwwwwwwwwwwww
%{
This function uses Newton’s method to solve Kepler’s
equation E - e*sin(E) ¼ M for the eccentric anomaly,
given the eccentricity and the mean anomaly.
E - eccentric anomaly (radians)
e - eccentricity, passed from the calling program
M - mean anomaly (radians), passed from the calling program
pi - 3.1415926...
User m-functions required: none
%}
% ----------------------------------------------
%...Set an error tolerance:

% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww