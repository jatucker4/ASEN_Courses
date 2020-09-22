function FOS = factor_of_safety(desired)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is to return a factor of safety given a
% desired probability
%
% Input: Desired probability
% Output: Factor of Safety
%
% Author: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the desired force
f_dsr = icdf('normal', desired, 4.8, 0.4);
%Using the desired force calculate a factor n
n = (4.8-f_dsr)/0.4;
%Then using n find the factor of safety
FOS = 1/(1-(n*(0.4/4.8)));
end

