clear
clc

for whichDataFile = 1:3
    
    if whichDataFile == 1
        data = tableread('lab1_8.txt');
        solidWorksData = csvread('90rpmAng.csv');
        displacementData = csvread('90rpmSmall.csv');
    elseif whichDataFile == 2
        data = textscan('lab1_10.txt');
        solidWorksData = csvread('113rpmAng.csv');
        displacementData = csvread('113rpmSmall.csv');
    else
        data = textscan('lab1_12.txt');
        solidWorksData = csvread('142rpmAng.csv');
        displacementData = csvread('142rpmSmall.csv');
    end
    
    %% Find RPMs from the experimental data
    onesandzeros = table2array(data(:,8));
    times = table2array(data(:,1));
    timeOfRotation = [];
    bobap = 0;
    newZeroes = true;
    index = 1;
    for j = 1:length(onesandzeros)
       if onesandzeros(j) > 0 && newZeroes == true
          bobap = bobap + 1;
          timeOfRotation(index) = times(j);
          index = index + 1;
          newZeroes = false;
       end
       if onesandzeros(j) < 1
           newZeroes = true;
       end
    end

    differences = [];
    for op = 1:length(timeOfRotation)
        if op > 1
            differences(op-1) = timeOfRotation(op) - timeOfRotation(op - 1);
        end
    end

    avgOfDifferences = mean(differences);

    RPMs = (1/avgOfDifferences)*60;
    
    %% Find the RPM from the Solidworks Data
    rotationIndex = 0;
    isNegative = false;
    start = true;
    numHalfRot = 0;
    foundNinety = false;
    minVolumeIndex1 = 0;
    bobap = 0;
    rotationTimes = [];
    index = 1;
    newRotation = true;
    found90 = false;
    for i = 1:length(solidWorksData)
       if solidWorksData(i,2) > -2 && solidWorksData(i,2) < 2 && newRotation == true && i > 5
          bobap = bobap + 1;
          rotationTimes(index) = solidWorksData(i,1);
          index = index + 1;
          newRotation = false;
       elseif newRotation == false
           newRotation = true;
       end
       if solidWorksData(i,2) > -92 && solidWorksData(i,2) < -88 && found90 == false
           timeWhen90 = solidWorksData(i,1);
           found90 = true;
       end
    end

    differences = [];
    for op = 1:length(rotationTimes)
        if op > 1
            differences(op-1) = rotationTimes(op) - rotationTimes(op - 1);
        end
    end

    avgOfDifferences = mean(differences);

    solidworksRPM = (1/rotationTimes(1))*60;

    %% Find the displacement of the piston from the solidworks data
    j = 0;
    k = 0;
    maxes = [];
    mins = [];
    displacements = displacementData(:,2);
    for i = 1:length(displacements)
        if i > 2 && i < length(displacements)
            if displacements(i-1) < displacements(i) && displacements(i+1) < displacements(i)
                j = j + 1;
                maxes(j) = displacements(i);
            elseif displacements(i-1) > displacements(i) && displacements(i+1) > displacements(i)
                k = k + 1;
                mins(k) = displacements(i);
            end
        end
    end
    maxDisplacement = max(maxes) - min(mins);   %mm
    volumeChange = maxDisplacement*(1/10)*pi*(1.5/2)^2; %cc

    foundMin = false;
    minDispIndex = 0;
    for w = 1:length(displacements)
        if w > 1 && w < length(displacements) && displacements(w) < displacements(w-1) && displacements(w) < displacements(w+1) && foundMin == false
            minDispIndex = w;
            foundMin = true;
        end
    end

    minWorkingVolume = 342; %cc


    %% Find Pressure from the experimental data
    pressures = table2array(data(:,2));
    times = table2array(data(:,1));
    multiplier = 1;
    newPress = [];
    for w = 1:length(pressures)
        if w == 50*multiplier
            newPress(multiplier) = pressures(w);
            multiplier = multiplier + 1;
        end
    end

    foundMax = false;
    maxPressureIndex = 0;
    for w = 1:length(newPress)
        if w > 1 && w < length(newPress) && newPress(w) > newPress(w-1) && newPress(w) > newPress(w+1) && foundMax == false
            maxPressureIndex = w;
            foundMax = true;
        end
    end

    maxPressureTime = times(maxPressureIndex-1);


    %% Make the experimental time scale match up with the solidworks time scale and plot

    displacementsCutoff = displacements(minDispIndex:length(displacements));

    normalizeZeroDisp = displacementsCutoff - min(mins);

    volumes = normalizeZeroDisp*(1/10)*pi*(1.5/2)^2; %cc

    volumeTimes = displacementData(:,1);
    pressureTimes = table2array(data(:,1));

    maxTimeV = volumeTimes(length(volumeTimes));
    maxTimeP = pressureTimes(length(pressureTimes));

    endPressIndex = 0;
    for j = 1:length(pressureTimes)
        if maxTimeV > maxTimeP
            if volumeTimes(j) > maxTimeP
                endPressIndex = j;
                break
            end
        else
            if pressureTimes(j) > maxTimeV
                endPressIndex = j;
                break
            end
        end
    end

    pressuresFinal = pressures(1:endPressIndex);

    multiplier2 = length(pressuresFinal)/length(volumes);

    index = 1;
    pressuresAdjust = [];
    for w = 1:length(pressuresFinal)
        if w == fix(index*multiplier2)
            pressuresAdjust(index) = pressuresFinal(w);
            index = index + 1;
        end
    end

    figure
    plot(volumes,pressuresAdjust);

    oneCycleLength = 0;
    for solidWorksData = 1:length(volumes)
        if solidWorksData > 5 && volumes(solidWorksData) == 0
            oneCycleLength = solidWorksData;
            break
        end
    end

    figure
    plot(volumes(1:oneCycleLength),pressuresAdjust(1:oneCycleLength));

end
