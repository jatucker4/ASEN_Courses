clear
clc

for whichDataFile = 1:3
    
    if whichDataFile == 1
        data = Bs('lab1_8.txt');
        solidWorksData = csvread('90rpmAng.csv',2,0);
        displacementData = csvread('90rpmSmall.csv',2,0);
    elseif whichDataFile == 2
        data = Bs('lab1_10.txt');
        solidWorksData = csvread('113rpmAng.csv',2,0);
        displacementData = csvread('113rpmSmall.csv',2,0);
    else
        data = Bs('lab1_12.txt');
        solidWorksData = csvread('142rpmAng.csv',2,0);
        displacementData = csvread('142rpmSmall.csv',2,0);
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
    deg_per_second = RPMs *(180/30);
    time_180 = 180/deg_per_second;
    time_270 = 270/deg_per_second;
    time_90 = 90/deg_per_second;
    time_0 = 0;
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
    avg_current = mean(table2array(data(:,7)));
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
    
    if whichDataFile ~= 3
        pressuresFinal = pressures(1:endPressIndex);
        multiplier2 = length(pressuresFinal)/length(volumes);
    else
%         pressuresFinal = volumes(1:endPressIndex);
        pressuresFinal = pressures(1:endPressIndex);
%         multiplier2 = length(pressuresFinal)/length(pressuresFinal);
        multiplier2 = 23124/length(volumes);
    end
    
    index = 1;
    pressuresAdjust = [];
    for w = 1:length(pressuresFinal)
        if w == fix(index*multiplier2)
            pressuresAdjust(index) = pressuresFinal(w);
            index = index + 1;
        end
    end
    temp1 = find(volumes<1*10^-3,2);
    volumes = volumes(1:temp1(2)) * 1e-6;
%     index_180 = find(volumeTimes>time_180,1) - 1;
%     index_270 = find(volumeTimes>time_270,1) - 1;
%     index_90 = find(volumeTimes>time_90,1) - 1;
%     index_0 = 1;
    pressuresAdjust = pressuresAdjust(1:temp1(2));
    if whichDataFile == 3
        pressuresAdjust = (-pressuresAdjust(1:temp1(2)));
    end
    pressuresAdjust = pressuresAdjust' * 6894.7600026396;
    index_180 = find(pressuresAdjust>360,1);
    index_270 = find(pressuresAdjust<-120,1);
    index_90 = find(pressuresAdjust>100,1);
    index_0 = find(pressuresAdjust<-350,1);
    figure
    plot(volumes,pressuresAdjust);
    ylabel("$Pressure\:[Pa]$",'Interpreter','latex','FontSize',26)
    xlabel("$Volume\:[m^3]$",'Interpreter','latex','FontSize',26)
    title("$Pressure\:Versus\:Volume\:Graph$",'Interpreter',...
        'latex','FontSize',26)
%     [k,work1] = boundary(volumes,pressuresAdjust);
    if index_270 < index_180
        work_out = polyarea(volumes(index_270:index_180),pressuresAdjust(index_270:index_180));
    else 
        work_out = polyarea(volumes(index_180:index_270),pressuresAdjust(index_180:index_270));
    end
    if index_0 > index_90
        work_in = polyarea(volumes(index_90:index_0),pressuresAdjust(index_90:index_0));
    else
        work_in = polyarea(volumes(index_0:index_90),pressuresAdjust(index_0:index_90));        
    end
    work2 = polyarea(volumes,pressuresAdjust);
    fprintf("The work done in trial %d is %f\n",whichDataFile,work2);
    fprintf("The work out of the system is %f for trial %d\n",work_out,whichDataFile);
    fprintf("The work in to the system is %f for trial %d\n",work_in,whichDataFile);
%     oneCycleLength = 0;
%     for solidWorksData = 1:length(volumes)
%         if solidWorksData > 5 && volumes(solidWorksData) == 0
%             oneCycleLength = solidWorksData;
%             break
%         end
%     end
% 
%     figure
%     plot(volumes(1:oneCycleLength),pressuresAdjust(1:oneCycleLength));

end
