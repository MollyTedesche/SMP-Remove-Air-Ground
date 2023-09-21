%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Script Module for Systematically Removing Top Air Space and Bottom Ground Space
%%%% from SnowMicroPenetrometer Raw Force Profiles 
%%%% using statistical approach based on Satyawali et al. (2009) paper:
%%%% "Preliminary characterization of Alpine snow using SnowMicroPen"

function x = RemoveTopBottom_SMP(zF,force)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART (1) - Algorithm to remove top of SMP Profile that indicates instrument moving through  
    % air space before reaching the snow surface (top 10 mm) - from Satyawali et al., 2009

    % Rule 1 - first 10 mm of the SMP profile is marked as air because cone tip of the SMP travels 
    % about this length in air before hitting the snow surface. Then estimate MPF, SD, & CV values 
    % averaged over 10 mm profile length (Mean Force, Standard Deviation, Coefficient of Variation)

    profileSMP = [zF,force];
    
    % remove first force msmnt at zero - sometimes erroneously large and could change stats
    profileSMP(1,:) = [];
    
    % Find top 10 mm of SMP profile
    topAirBin = profileSMP(profileSMP(:,1)>=0 & profileSMP(:,1)<11, :);

    % Find statistics of top 10 mm of air
    %disp('Mean Force, first 10 mm')
    MPF = mean(topAirBin(:,2)); 
    %disp('STD, first 10 mm')
    SD = std(topAirBin(:,2)); 
    %disp('Coef of Var, first 10 mm')
    CV = SD/MPF;


    % Rule 2 - from 11 mm onward till the end of SMP profile, if statistics (MPF,SD,CV) of msmnt 
    % are within the range of ±10% of the first 10 mm, it's marked as air. The first set of 10, 1mm 
    % consecutive non-air measurements is set as snow surface (snow-air interface is 10 mm thick)

    % Mark top 10 mm as air and evaluate the rest of the SMP profile
    statusMatrixTop = repmat("air",height(topAirBin),1);

    SMPprof11 = profileSMP(profileSMP(:,1)>=11, :);
    
    % Find ±10% ranges for MPF,SD,CV
    lowMPF = MPF-(MPF*0.1);
    hiMPF = MPF+(MPF*0.1);
    lowSD = SD-(SD*0.1);
    hiSD = SD+(SD*0.1);
    lowCV = CV-(CV*0.1);
    hiCV = CV+(CV*0.1);

    % Find statistics for each 1 mm from 11 mm till end of profile and mark as "air" or "non-air"

    % Split profile into 1 mm sections and get statistics
    statusMatrix = [];
    
    %%%%% Preallocating the matrices doesn't seem to work here
    mmBottom = max(SMPprof11(:,1));

    MPFstatus = ("");
    SDstatus = ("");
    CVstatus = ("");
    
    for i = 11:(mmBottom)
        oneMM = SMPprof11(SMPprof11(:,1)>=i & SMPprof11(:,1)<(i+1) , :);
        mmMPF = mean(oneMM(:,2));
        mmSD = std(oneMM(:,2));
        mmCV = mmSD/mmMPF;

        if (lowMPF <= mmMPF) 
            if (hiMPF >= mmMPF)
                MPFstatus = "air";
            else 
                MPFstatus = "non-air";
            end
        end 

        if (lowSD <= mmSD)
            if (hiSD >= mmSD)
                SDstatus = "air";
            else 
                SDstatus = "non-air";
            end
        end 

        if (lowCV <= mmCV) 
            if (hiCV >= mmCV)
                CVstatus = "air";
            else 
                CVstatus = "non-air";
            end 
        end

        if (MPFstatus ~= "air") || (SDstatus ~= "air") || (CVstatus ~= "air")
            status = "non-air";
        else 
            status = "air";
        end    

        airColumn = repmat(status,height(oneMM),1);

        %%%%%%% Tried not to resize matrix; tried preallocating, but this doesn't work: 
        % statusRangeStart = ((i-11)*height(oneMM)) + 1;
        % statusRangeEnd = (i-10)*height(oneMM);
        % statusMatrix(statusRangeStart:statusRangeEnd,1) = airColumn;

        statusMatrix = [statusMatrix;airColumn];
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART (2) - Algorithm to remove bottom of SMP Profile that indicates below snow-ground
    % interface (also from Satyawali et al., 2009)

    % Note that in Satyawali et al., 2009 - the densest/hardest snow types
    % of "Melt-Freeze" and "Hard Crust" have averages of MPF = 6.5N, SD = 1.1N, and CV = 0.24; 
    % and maximums of MPF = 23.5N, SD = 5.7N, and CV = 2.0  

    % Rule 1 - Mark the last/bottom 10 mm of the SMP profile as ground, then estimate MPF, SD, & CV  
    % values averaged over last 10 mm profile length (Mean Force, Standard Dev, Coeff of Variation)

    % Find bottom 10 mm of SMP profile
    grndTop = mmBottom - 10;
    bottomGrndBin = profileSMP(profileSMP(:,1)>=grndTop & profileSMP(:,1)<=mmBottom, :);

    % Find statistics of bottom 10 mm of ground
    grndMPF = mean(bottomGrndBin(:,2)); 
    grndSD = std(bottomGrndBin(:,2)); 
    grndCV = grndSD/grndMPF;

    % Mark bottom 10 mm as ground
    statusMatrixBottom = repmat("ground",height(bottomGrndBin),1);


    % Rule 2 - From 11 mm from the bottom upward, if statistics (MPF,SD,CV) over each ascending 
    % 1 mm length are within the range of ±10% of the last 10mm, it's marked as ground. The first
    % set of 10, 1 mm ascending consecutive non-ground measurements (starting from the
    % bottom up) is set as bottom of the snowpack (snow-ground interface as 10mm thick)

    SMPprof0 = profileSMP(profileSMP(:,1)<grndTop, :);
    
    % Find ±10% ranges for ground MPF,SD,CV
    lowMPFgrnd = grndMPF-(grndMPF*0.1);
    hiMPFgrnd = grndMPF+(grndMPF*0.1);
    lowSDgrnd = grndSD-(grndSD*0.1);
    hiSDgrnd = grndSD+(grndSD*0.1);
    lowCVgrnd = grndCV-(grndCV*0.1);
    hiCVgrnd = grndCV+(grndCV*0.1);

    % Find statistics for each 1 mm from SMP profile top till top  
    % of pre-determined ground layer and mark as "ground" or "non-ground"

    % Split profile into 1 mm sections and get statistics
    statusMatrixGrnd = [];

    MPFstatusGrnd = ("");
    SDstatusGrnd = ("");
    CVstatusGrnd = ("");

    for i = 1:(grndTop+1)
        oneMMgrnd = SMPprof0(SMPprof0(:,1)>(i-1) & SMPprof0(:,1)<=i, :);
        mmMPFgrnd = mean(oneMMgrnd(:,2));
        mmSDgrnd = std(oneMMgrnd(:,2));
        mmCVgrnd = mmSDgrnd/mmMPFgrnd;

        if (lowMPFgrnd <= mmMPFgrnd) 
            if (hiMPFgrnd >= mmMPFgrnd)
                MPFstatusGrnd = "ground";
            else 
                MPFstatusGrnd = "non-ground";
            end
        end 

        if (lowSDgrnd <= mmSDgrnd)
            if (hiSDgrnd >= mmSDgrnd)
                SDstatusGrnd = "ground";
            else 
                SDstatusGrnd = "non-ground";
            end
        end 

        if (lowCVgrnd <= mmCVgrnd) 
            if (hiCVgrnd >= mmCVgrnd)
                CVstatusGrnd = "ground";
            else 
                CVstatusGrnd = "non-ground";
            end 
        end

        if (MPFstatusGrnd ~= "ground") || (SDstatusGrnd ~= "ground") || (CVstatusGrnd ~= "ground")
            statusGrnd = "non-ground";
        else 
            statusGrnd = "ground";
        end    

        grndColumn = repmat(statusGrnd,height(oneMMgrnd),1);
        
        %%%%%%% Tried not to resize matrix; tried preallocating, but this doesn't work: 
        % statusRangeStart = ((i-1)*height(oneMM)) + 1;
        % statusRangeEnd = i*height(oneMM);
        % statusMatrixGrnd(statusRangeStart:statusRangeEnd,1) = grndColumn;
        
        statusMatrixGrnd = [statusMatrixGrnd;grndColumn];
    end
    
    x.statusAir = [statusMatrixTop;statusMatrix];
    x.statusGrnd = [statusMatrixGrnd;statusMatrixBottom]; 
    x.zFstatus = profileSMP(:,1);
    x.forceStatus = profileSMP(:,2);  
    
    %%%%%%%%%%%%%% FIND SNOW SURFACE %%%%%%%%%%%%%%
    %%% Find first set of 10mm consecutive non-air msrmnts (1mm ~ 240 lines; 10 mm ~ 2400 lines)
    %%% index of first line of 10 mm set of measurements is snow surface
    binaryAir = contains(x.statusAir,"non-air");
    binaryAirVector = reshape(binaryAir,1,[]);
    trans = find(diff([0,binaryAirVector,0]==1)); % Transitions 0->1 and 1->0 indices
    p = trans(1:2:end-1);  % Start indices
    y = trans(2:2:end)-p;  % Consecutive ones count
    tenMM = find(y>=2400);
    snowSurfIndx = p(1,tenMM(1,1));
    x.snowSurfVal = x.zFstatus(snowSurfIndx,1);
    
    %%%%%%%%%%%%%% FIND GROUND SURFACE %%%%%%%%%%%%%%
    %%% Find last set of 10mm consecutive non-ground msrmnts (1mm ~ 240 lines; 10 mm ~ 2400 lines)
    %%% index of last line above 10 mm set of ground measurements is snow bottom
    binaryGrnd = contains(x.statusGrnd,"non-ground");
    binaryGrndVector = reshape(binaryGrnd,1,[]);
    transGrnd = find(diff([1,binaryGrndVector,1]==0)); % Transitions 1->0 and 0->1 indices
    a = transGrnd(1:2:end-1);  % Start indices
    c = transGrnd(2:2:end)-a;  % Consecutive ones count
    tenMMGrnd = find(c>=2400);

    % Only use default "ground" layer (statusMatrixBottom), if there are no additional 10mm 
    % or thicker ground sections found statistically above the manually added (default) layer 
    numTenMM = numel(tenMMGrnd);
    if numTenMM == 1
        snowBtmIndx = (a(1,tenMMGrnd(1,1)))-1;
    else
        revTenMM = flip(tenMMGrnd); % reverse order of 10mm non-ground indices to find last set
        snowBtmIndx = (a(1,revTenMM(1,2)))-1;
    end
    x.snowBtmVal = x.zFstatus(snowBtmIndx,1);

    %%%%% Crop the depth and force vectors in the SMP profile
    x.zFcrop = profileSMP(snowSurfIndx:snowBtmIndx,1);  
    x.forceCrop = profileSMP(snowSurfIndx:snowBtmIndx,2);  

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Optional: Find Max Normal Forces to calculate density later during SMP Inversion 
    forceNorm = x.forceCrop*2;
    x.forceMax = NaN(height(forceNorm),1);

    for i = 2:((height(forceNorm))-1)
        if forceNorm(i) > forceNorm(i-1) 
            if forceNorm(i) > forceNorm(i+1)
                maxForce = forceNorm(i);
            else 
                maxForce = NaN;
            end
            x.forceMax(i,:) = maxForce;
        end
    end
end
