function [I_b,I_f] = beammoment(width)
h_f = 0.75 / 39.37;
h_b = (1/32) / 39.37;
I_f = (1/12)*width*(h_f)^3;
I_b = 2 * ((1/12) * width * (h_b)^3 + (width * h_b) * (0.5 * (h_f + h_b))^2);
end