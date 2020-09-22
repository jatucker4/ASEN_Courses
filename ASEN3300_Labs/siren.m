function S = siren(A,O,F)
%PRELAB_9 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(A)
    if F(i) || (A(i) && O(i))
        S(i) = 1;
    else
        S(i) = 0;
    end
end
end

