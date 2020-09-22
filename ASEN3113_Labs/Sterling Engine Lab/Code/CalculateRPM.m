m = csvread('Flywheel-1.csv');

rotationIndex = 0;
isNegative = false;
start = true;
numHalfRot = 0;
for i = 1:length(m)
    if isNegative == false && m(i,2) < 0 && start == true;
        isNegative = true;
        start = false;
    elseif isNegative == false && m(i,2) < 0
        numHalfRot = numHalfRot + 1;
        rotationIndex = i;
    elseif isNegative == true && m(i,2) > 0
        numHalfRot = numHalfRot + 1;
        rotationIndex = i;
    end
end

numRot = numHalfRot*2;

rpm = numRot/m(rotationIndex,1);