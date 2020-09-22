function [r] = ellipticalOrbit(e, a, theta)
%{
    This function calculates the distance between an orbiting object and
    its central body when it is at the angle theta up from the local
    horizontal.
    Inputs:
        e     = ellipticity of orbit
        a     = semi-major axis of orbit
        theta = angle between orbiting object and the local horizontal
    Outputs:
        r     = distance between orbiting object and central body
%}
    r = (a)/(1-(e*cos(theta)));
    
end