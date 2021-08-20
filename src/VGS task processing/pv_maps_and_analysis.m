% Heat Map Difference Plots
% Loading pre and post data, and then creating a difference in peak
% velocity map on to the visual field in cartesian coordinates and also
% onto the visual field in polar coordinates and also on to the SC Map to
% estimate the spread of muscimol
clear all
close all

%-----------MODIFY MAIN SETTINGS HERE-------------------------------------%
xx = 100; %receptive field (RF) in horizontal degrees in cartesian coordinates
yy = -60; %RF in vertical degrees in cartesian coordinates
target_in_RF = 1; %1 if one of the visually guided targets was in the RF, 0 if not
%-------------------------------------------------------------------------%
%% Loading File 1 (pre) and File 2 (post)
%should have the following data for each valid trial in the visually guided
%saccade task: target location, reaction time and/or peak velocity
addpath(genpath('C:\Users\eliza\Desktop\Matlab Scripts\'))

%receptive field for plotting on visual field
plotreceptivefieldx = xx + 201;
plotreceptivefieldy = abs(yy - 201);

%file loading parameters
continueAction_pre = 1;
continueAction_post = 1;
files_inputted_pre = 1;
files_inputted_post = 1;

%For Peak Velocity (PV) Difference Map...

%Load File 1 (Pre-Injection Data)
while continueAction_pre
    display('-----Open pre files------');
    [fileNamepre, pathnamepre] = getFile('*.nev');
    filenamepre = fileNamepre(1,1:end-4);
    filePathpre = ['Z:\Lab Staff Folders\Liz\Matlab Scripts\Latency_Velocity_Map\datafiles\' filenamepre '_filtered.mat'];
    % filePathpre = [pwd 'datafiles\' filenamepre '_filtered.mat'];
    datapre{files_inputted_pre} = load(filePathpre);
    if files_inputted_pre == 1
        Trial_3_pre = datapre{files_inputted_pre}.Trial_3;
    else
        Trial_3_pre = [Trial_3_pre, datapre{files_inputted_pre}.Trial_3];
    end
    display('Last opened file was');
    fileNamepre
    prompt_0 = ('Open another pre session? (y/n)');
    add_file = input(prompt_0, 's');
    if add_file == 'Y' || add_file == 'y'
        continueAction_pre = 1;
        files_inputted_pre = files_inputted_pre +1;
    elseif add_file == 'N' || add_file == 'n'
        continueAction_pre = 0;
    end
    
end


%Load File 2 (Post-Injection Data)
while continueAction_post
    display('-----Open post files------');
    [fileNamepost, pathnamepost] = getFile('*.nev');
    filenamepost = fileNamepost(1,1:end-4);
    filePathpost = ['Z:\Lab Staff Folders\Liz\Matlab Scripts\Latency_Velocity_Map\datafiles\' filenamepost '_filtered.mat'];
    % filePathpre = [pwd 'datafiles\' filenamepost '_filtered.mat'];
    datapost{files_inputted_post} = load(filePathpost);
    if files_inputted_post == 1
        Trial_3_post = datapost{files_inputted_post}.Trial_3;
    else
        Trial_3_post = [Trial_3_post, datapost{files_inputted_post}.Trial_3];
    end
    display('Last opened file was');
    filenamepost
    prompt_0 = ('Open another post session? (y/n)');
    add_file = input(prompt_0, 's');
    if add_file == 'Y' || add_file == 'y'
        continueAction_post = 1;
        files_inputted_post = files_inputted_post +1;
    elseif add_file == 'N' || add_file == 'n'
        continueAction_post = 0;
    end
    
    
end
filename = filenamepre(1:8);

% Averaging the peak velocity with the same target location values and
% getting a count of how many trials for each target

%First, create a list of all the targets in the datafile for pre
rfindatapre = struct('Tx', [500], 'Ty', [500], 'AveragePV', [500], 'TrialCount', [500]);              %(500 is an arbitrary number as a filler to make the code work)
p = 1;
for r = 1:numel(Trial_3_pre)        %for each trial
    currentx = Trial_3_pre(r).RFx;
    currenty = Trial_3_pre(r).RFy;
    duplicate = 0;               %set duplicate to 0 for every new trial's RF
    for f = 1:numel([rfindatapre.Tx])    %go through the rf in data list to check to see if the RF data value is already listed
        if rfindatapre(f).Tx ~= 500     %if it's not the first index
            if currentx == rfindatapre(f).Tx && currenty == rfindatapre(f).Ty
                duplicate = duplicate + 1;
            else
            end
        end
    end
    
    if duplicate == 0     %if this RF hasn't been listed yet, then list it in the rfindatapre
        rfindatapre(p).Tx = currentx;
        rfindatapre(p).Ty = currenty;
        p = p +1;   %go to the next empty row to fill in
    else
    end
    
    
end

checkpoint = 1;
%Then, take the average and median of the RF data points (rt and pv) with the same RF's
for f = 1:numel([rfindatapre.Tx])    %going through the list of all the unique RF's
    currentx = rfindatapre(f).Tx;
    currenty = rfindatapre(f).Ty;
    samerfindx = find([Trial_3_pre.RFx] == currentx & [Trial_3_pre.RFy] == currenty);   %find which trials have the same RF's
    samerfpv = [Trial_3_pre(samerfindx).PeakVelocity];                                                    %get a list of the PV's in the trials with the same RF's
    samerfpv(find(isnan(samerfpv))) = [];
    avgpvpre = mean(samerfpv);                                                          %get average of the pv
    medpvpre = median(samerfpv);
    rfindatapre(f).AveragePV = avgpvpre; %save it in the structure
    rfindatapre(f).StdPV = std(samerfpv);
    rfindatapre(f).MedianPV = medpvpre;
    samerfrt = [Trial_3_pre(samerfindx).RT];
    samerfrt(find(isnan(samerfrt))) = [];
    avgrtpre = mean(samerfrt);
    rfindatapre(f).AverageRT = avgrtpre;
    rfindatapre(f).StdRT = std(samerfrt);
    rfindatapre(f).TrialCount = numel(samerfindx);
end
done =1 ;


%First, create a list of all the RF's in the datafile for post
rfindatapost = struct('Tx', [500], 'Ty', [500], 'AveragePV', [500], 'TrialCount', [500]);
p = 1;
for r = 1:numel(Trial_3_post)                                               %for each trial
    currentx = Trial_3_post(r).RFx;                                         %get the Tx for that trial
    currenty = Trial_3_post(r).RFy;                                         %get the Ty for that trial
    duplicate = 0;               %set duplicate to 0 for every new trial's target
    for f = 1:numel([rfindatapost.Tx])    %go through the targets in data list to check to see if the target data value is already listed
        if rfindatapost(f).Tx ~= 500     %if it's not the first index
            if currentx == rfindatapost(f).Tx && currenty == rfindatapost(f).Ty
                duplicate = duplicate + 1;
            else
            end
        end
    end
    
    if duplicate == 0     %if this target hasn't been listed yet, then list it in the rfindatapost
        rfindatapost(p).Tx = currentx;
        rfindatapost(p).Ty = currenty;
        p = p +1;   %go to the next empty row to fill in
    else
    end
end
%Then, take the average of the peak velocity (PV) and the reaction time (RT) of all
%the trials with the same target
for f = 1:numel([rfindatapost.Tx])    %going through the list of all the unique targets
    currentx = rfindatapost(f).Tx;
    currenty = rfindatapost(f).Ty;
    samerfindx = find([Trial_3_post(1:end).RFx] == currentx & [Trial_3_post(1:end).RFy] == currenty);   %find which trials have the same targets
    samerfpv = [Trial_3_post(samerfindx).PeakVelocity];                                                    %get a list of the PV's in the trials with the same targets
    samerfpv(find(isnan(samerfpv))) = [];
    avgpvpost = mean(samerfpv);                                                          %get average of the pv
    medpvpost = median(samerfpv);
    rfindatapost(f).AveragePV = avgpvpost;  %save it in the structure
    rfindatapost(f).StdPV = std(samerfpv);
    rfindatapost(f).MedianPV = medpvpost;
    samerfrt = [Trial_3_post(samerfindx).RT];
    samerfrt(find(isnan(samerfrt))) = [];
    avgrtpost = mean(samerfrt);
    rfindatapost(f).AverageRT = avgrtpost;
    rfindatapost(f).StdRT = std(samerfrt);
    rfindatapost(f).TrialCount = numel(samerfindx);
end

%first checking if there are missing RF data points in pre or post data
%files, or if an RF datapoint has less than 5 trials
if numel([rfindatapost.Tx]) < 32 || numel([rfindatapost.Tx]) < 32 || sum([rfindatapost.TrialCount] < 5) >= 1 || sum([rfindatapre.TrialCount] < 5) >= 1
    prompt_4 = ['\n Error: Post and Pre injection data dont have the same number of RF data points OR theres less than 5 trials in a data point. Heres a list of them...'...
        '\n For Pre-injection...'...
        '\n Num  ', num2str(1:numel([rfindatapre.Tx])),...
        '\n Tx  ', num2str([rfindatapre.Tx]),...
        '\n Ty  ', num2str([rfindatapre.Ty]), ...
        '\n Count', num2str([rfindatapre.TrialCount]), ...
        '\n'...
        '\n For Post-injection...'...
        '\n Num  ', num2str(1:numel([rfindatapost.Tx])),...
        '\n Tx  ', num2str([rfindatapost.Tx]),...
        '\n Ty  ', num2str([rfindatapost.Ty]), ...
        '\n Count', num2str([rfindatapost.TrialCount]), ...
        '\n'...
        'Continue? y/n'];
    input_prompt_4 = input(prompt_4, 's');
    if input_prompt_4 == 'n'
        dbstop in 'peak_velocity_maps_plotter' at '204';
    else
    end
end

% Substracting average PV first, and then normalizing, then interpolating
DiffPVStruct = struct('Tx', [500], 'Ty', [500], 'DiffPV', [500], 'NormalizedDiffPV', [500], 'DiffRT', [500], 'Normalized_RT', [500]);

onlyoneTxpre = [];
onlyoneTypre = [];
%Subtracting average PV
for f = 1:numel([rfindatapre.Tx])        %going through all of the list of targets in pre and matching it's index/data to post
    currentprex = rfindatapre(f).Tx;
    currentprey = rfindatapre(f).Ty;
    matchingrfindx = find([rfindatapost.Tx] == currentprex & [rfindatapost.Ty] == currentprey);   %get the index in post data struct where it has the same target
    if isempty(matchingrfindx)
        onlyoneTxpre = [onlyoneTxpre, currentprex];
        onlyoneTypre = [onlyoneTypre, currentprey];
    else
        diffPV = rfindatapost(matchingrfindx).AveragePV - rfindatapre(f).AveragePV;   %subtracting post - pre for PV average
        diffPVmed = rfindatapost(matchingrfindx).MedianPV - rfindatapre(f).MedianPV;
        diffRT = rfindatapost(matchingrfindx).AverageRT - rfindatapre(f).AverageRT;
        DiffPVStruct.Tx(f) = currentprex; %saving the target locations in the new structure along with the PV data
        DiffPVStruct.Ty(f) = currentprey;
        DiffPVStruct.DiffPV(f) = diffPV;
        DiffPVStruct.NormalizedDiffPV(f) = (diffPV/rfindatapre(f).AveragePV)*100; %normalizing it by (post - pre)/(pre) PV. The difference normalized to the pre-injection PV value
        DiffPVStruct.DiffPVMed(f) = diffPVmed;
        DiffPVStruct.DiffRT(f) = diffRT;
    end
end

%checking to see if there's any missing data for targets in pre or post
onlyoneTxpost = [];
onlyoneTypost = [];
for f = 1:numel([rfindatapost.Tx])        %going through all of the list of targets in pre and matching it's index/data to post
    currentpostx = rfindatapost(f).Tx;
    currentposty = rfindatapost(f).Ty;
    matchingrfindx = find([rfindatapre.Tx] == currentpostx & [rfindatapre.Ty] == currentposty);   %get the index in post data struct where it has the same target
    if isempty(matchingrfindx)
        onlyoneTxpost = [onlyoneTxpost, currentpostx];
        onlyoneTypost = [onlyoneTypost, currentposty];
    end
end
remove_indexpre = [];
remove_indexpost = [];
if numel(onlyoneTxpre) > 0 || numel(onlyoneTypre) > 0 || numel(onlyoneTxpost) > 0 || numel(onlyoneTypost) > 0
    prompt_4 = ['\n Error: Theres missing target datapoints. Heres a list of them...'...
        '\n For Pre-injection...'...
        '\n Num  ', num2str(1:numel(onlyoneTxpost)),...
        '\n Tx  ', num2str(onlyoneTxpost),...
        '\n Ty  ', num2str(onlyoneTypost), ...
        '\n'...
        '\n For Post-injection...'...
        '\n Num  ', num2str(1:numel(onlyoneTxpre)),...
        '\n Tx  ', num2str(onlyoneTxpre),...
        '\n Ty  ', num2str(onlyoneTypre), ...
        '\n'...
        'Continue and remove these points from both data sets? y/n'];
    input_prompt_4 = input(prompt_4, 's');
    if input_prompt_4 == 'n'
        dbstop in 'peak_velocity_maps_plotter' at '259';
        stop = 1;
    elseif input_prompt_4 == 'y'
        for i = 1:numel(onlyoneTxpre)
            currentprex = onlyoneTxpre(i);
            currentprey = onlyoneTypre(i);
            remove_indexpre = [remove_indexpre, find([rfindatapre.Tx] == currentprex & [rfindatapre.Ty] == currentprey)];
        end
        for i = 1:numel(onlyoneTxpost)
            currentpostx =onlyoneTxpost(i);
            currentposty =onlyoneTypost(i);
            remove_indexpost = [remove_indexpost, find([rfindatapost.Tx] == currentpostx & [rfindatapost.Ty] == currentposty)];
        end
        rfindatapre(remove_indexpre) = [];
        rfindatapost(remove_indexpost) = [];
        
    end
end


% Plotting Peak velocity map onto the visual field in cartesian coordinates

Tx=DiffPVStruct.Tx+201;    %RF values for indexing 401x401 matrix --> 20 x 20 degrees visual field
Ty=abs(DiffPVStruct.Ty-201);

Txpre = [rfindatapre.Tx]+201;
Typre = abs([rfindatapre.Ty]-201);

Txpost = [rfindatapost.Tx] +201;
Typost = abs([rfindatapost.Ty]-201);

%Linear Interpolation of peak velocity data
[Xq,Yq] = meshgrid(1:401,1:401);

FDiffPV_norm = scatteredInterpolant(Tx',Ty', [DiffPVStruct.NormalizedDiffPV]');
interp_DiffPVmatrix_norm = FDiffPV_norm(Xq,Yq);

%Smoothing
smooth_DiffPVmatrix1_norm = smooth2a(interp_DiffPVmatrix_norm, 20, 20);

close all
figure(1);
imagesc(smooth_DiffPVmatrix1_norm)
set(gca, 'XTick', [ 0 50 100 150 200 250 300 350 400  ], 'XTickLabel', {'-20', '-15', '-10', '-5', '0', '5', '10', '15', '20'});
set(gca, 'YTick', [ 0 50 100 150 200 250 300 350 400  ], 'YTickLabel', {'20', '15', '10', '5', '0', '-5', '-10', '-15', '-20'});
hold on
scatter(Tx,Ty, 20, 'filled', 'k');
scatter(plotreceptivefieldx, plotreceptivefieldy, 150, 'x', 'w', 'linewidth', 2);

title(['Difference in Peak Velocity -- ',filenamepre(1:8)]);
ylabel('Vertical Degrees');
xlabel('Horiztonal Degrees');
c = colorbar;
caxis([-60 60]);  %set the colorbar scale min and max values

%% Getting the information to plot the peak velocity map onto the visual field in polar coordinates and onto the SC Map
datastructs = {'rfindatapre', 'rfindatapost'};

for i = 1:numel(rfindatapre)
    rfindatapre(i).Tx = rfindatapre(i).Tx/10; %divide by 10 because REX degree angles is 10x greater than regular degree
    rfindatapre(i).Ty = rfindatapre(i).Ty/10;
    rfindatapre(i).angle = (atan2(rfindatapre(i).Ty, rfindatapre(i).Tx)) * 180/pi;    % get the polar angle in degrees
    if rfindatapre(i).angle < 0
        rfindatapre(i).angle = 360 + rfindatapre(i).angle;
    end
    rfindatapre(i).amplitude = sqrt((rfindatapre(i).Tx^2) + (rfindatapre(i).Ty^2));   %amplitude for the saccade vector
    
    %do same for post data
    rfindatapost(i).Tx = rfindatapost(i).Tx/10; %divide by 10 because REX degree angles is 10x greater than regular degree
    rfindatapost(i).Ty = rfindatapost(i).Ty/10;
    rfindatapost(i).angle = (atan2(rfindatapost(i).Ty, rfindatapost(i).Tx)) * 180/pi;    % get the polar angle in degrees
    if rfindatapost(i).angle < 0
        rfindatapost(i).angle = 360 + rfindatapost(i).angle;
    end
    rfindatapost(i).amplitude = sqrt((rfindatapost(i).Tx^2) + (rfindatapost(i).Ty^2));   %amplitude for the saccade vector
end

% the script used to plot the peak velocity difference onto the visual
% field with polar coordinates and SC Map
SCMap(xx, yy, rfindatapre, rfindatapost,filename)


%% Get the  diffPV value from the interpolation at the RF -> get the diffPV interpolated value -> diffPV/PVpre
% For experiments that did not have the RF in the targets; before 10/1/19

if ~target_in_RF
    [Xq,Yq] = meshgrid(1:401,1:401);
    
    %For PV - mean
    FPVpre = scatteredInterpolant(RFxpre',RFypre', rfindatapre.AveragePV', 'linear', 'none');
    interp_PVprematrix = FPVpre(Xq,Yq);                                         %interpolated values for the pre-PV matrix
    PVRFpre = interp_PVprematrix(plotreceptivefieldy, plotreceptivefieldx);
    
    FPVpost = scatteredInterpolant(RFxpost',RFypost', rfindatapost.AveragePV','linear', 'none');
    interp_PVpostmatrix = FPVpost(Xq,Yq);                                         %interpolated values for the post-PV matrix
    PVRFpost = interp_PVpostmatrix(plotreceptivefieldy, plotreceptivefieldx);
    
    %For PV - median
    FPVpremed = scatteredInterpolant(RFxpre',RFypre', rfindatapre.MedianPV');
    interp_PVpremedmatrix = FPVpremed(Xq,Yq);                                         %interpolated values for the pre-PV matrix
    PVRFpremed = interp_PVpremedmatrix(plotreceptivefieldy, plotreceptivefieldx);
    
    FPVpostmed = scatteredInterpolant(RFxpost',RFypost', rfindatapost.MedianPV');
    interp_PVpostmedmatrix = FPVpostmed(Xq,Yq);                                         %interpolated values for the post-PV matrix
    PVRFpostmed = interp_PVpostmedmatrix(plotreceptivefieldy, plotreceptivefieldx);
    
    %For RT
    FRTpre = scatteredInterpolant(RFxpre',RFypre', rfindatapre.AverageRT');
    interp_RTprematrix = FRTpre(Xq,Yq);                                         %interpolated values for the pre-RT matrix
    RTRFpre = interp_RTprematrix(plotreceptivefieldy, plotreceptivefieldx);
    
    FRTpost = scatteredInterpolant(RFxpost',RFypost', rfindatapost.AverageRT');
    interp_RTpostmatrix = FRTpost(Xq,Yq);                                         %interpolated values for the post-RT matrix
    RTRFpost = interp_RTpostmatrix(plotreceptivefieldy, plotreceptivefieldx);
    
    % for the contralateral RF, PV's and RT's
    if xx > 0
        plotreceptivefieldx_contra = plotreceptivefieldx-(xx*2);
    elseif xx < 0
        plotreceptivefieldx_contra = plotreceptivefieldx+abs((xx*2));
    end
    plotreceptivefieldy_contra = plotreceptivefieldy;
    %     mean PV contra
    PVoutRFpre = interp_PVprematrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    PVoutRFpost = interp_PVpostmatrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    %     median PV contra
    PVoutRFpremed = interp_PVpremedmatrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    PVoutRFpostmed = interp_PVpostmedmatrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    %     mean RT contra
    RToutRFpre = interp_RTprematrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    RToutRFpost = interp_RTpostmatrix(plotreceptivefieldx_contra, plotreceptivefieldy_contra);
    
    
    DiffPVRF = interp_DiffPVmatrix_nonorm(plotreceptivefieldy,plotreceptivefieldx);    %interpolated (Pre - post PV) -- using this one for analysis
    DiffPVPercent_nonorm = DiffPVRF/PVRFpre;
    PVEffect_nonorm = DiffPVPercent_nonorm*100;
    
    DiffPVRF1 = PVRFpost - PVRFpre;    %interpolated Pre PV - interpolated Post PV
    DiffPVPercent_nonorm = DiffPVRF1/PVRFpre;
    PVEffect_nonorm1 = DiffPVPercent_nonorm*100;
    
    DiffRTRF = interp_DiffRTmatrix_nonorm(plotreceptivefieldy,plotreceptivefieldx);    %interpolated (Pre - post RT) -- using this one for analysis
    DiffRTPercent_nonorm = DiffRTRF/RTRFpre;
    RTEffect_nonorm = DiffRTPercent_nonorm*100;
    
    DiffRTRF1 = RTRFpost - RTRFpre;    %interpolated Pre RT - interpolated Post RT
    DiffRTPercent_nonorm = DiffRTRF1/RTRFpre;
    RTEffect_nonorm1 = DiffRTPercent_nonorm*100;
    
    %store data in the LMData structure; LMData structure will be later
    %processed all together for all experiments
    LMData = struct();
    LMData.DiffPVRF = DiffPVRF;
    LMData.DiffPVRF1 = DiffPVRF1;
    LMData.DiffRTRF = DiffRTRF;
    LMData.DiffRTRF1 = DiffRTRF1;
    LMData.DiffPercentPV = PVEffect_nonorm;
    LMData.DiffPercentPV1 = PVEffect_nonorm1;
    LMData.DiffPercentRT = RTEffect_nonorm;
    LMData.DiffPercentRT1 = RTEffect_nonorm1;
    
    LMData.MeanPVInRF{1,1} = PVRFpre;
    LMData.MeanPVInRF{2,1} = PVRFpost;
    LMData.MeanPVInRF{3,1} = nan;
    LMData.MedianPVInRF{1,1} = PVRFpremed;
    LMData.MedianPVInRF{2,1} = PVRFpostmed;
    LMData.MeanPVInRF{3,1} = nan;
    LMData.MeanPVOutRF{1,1} = PVoutRFpre;
    LMData.MeanPVOutRF{2,1} = PVoutRFpost;
    LMData.MeanPVOutRF{3,1} = nan;
    LMData.MedianPVOutRF{1,1} = PVoutRFpremed;
    LMData.MedianPVOutRF{2,1} = PVoutRFpostmed;
    LMData.MeanPVOutRF{3,1} = nan;
    
    LMData.MeanRTInRF{1,1} = RTRFpre;
    LMData.MeanRTInRF{2,1} = RTRFpost;
    LMData.MeanRTInRF{3,1} = nan;
    LMData.MeanRTOutRF{1,1} = RToutRFpre;
    LMData.MeanRTOutRF{2,1} = RToutRFpost;
    LMData.MeanRTOutRF{3,1} = nan;
    
    LMData.Trial_3{1,1} = Trial_3_pre;
    LMData.Trial_3{1,2} = Trial_3_post;
    LMData.Trial_3{1,3} = nan;
    LMData.RFx = xx;
    LMData.RFy = yy;
elseif target_in_RF
    %% Get the % diff PV direct datapoint (only for experiments 10/1/19 and beyond)
    % For experiment that did have the RF in the targets; 10/1/19 and beyond
    
    indxRFdiffPV = find(DiffPVStruct.RFx == xx & DiffPVStruct.RFy == yy);
    
    DiffPVRF = DiffPVStruct.DiffPV(indxRFdiffPV);
    DiffRTRF = DiffPVStruct.DiffRT(indxRFdiffPV);
    
    indxRFprePV = find(rfindatapre.RFx == xx & rfindatapre.RFy == yy);
    prePVRF = rfindatapre.AveragePV(indxRFprePV);
    prePVmedRF = rfindatapre.MedianPV(indxRFprePV);
    preRTRF = rfindatapre.AverageRT(indxRFprePV);
    
    PVEffect_nonorm = (DiffPVRF/prePVRF)*100;
    RTEffect_nonorm = (DiffRTRF/preRTRF)*100;
    
    %store data in the LMData structure; LMData structure will be later
    %processed all together for all experiments
    LMData = struct();
    LMData.DiffPVRF = DiffPVRF;
    LMData.DiffRTRF = DiffRTRF;
    LMData.DiffPercentPV = PVEffect_nonorm;
    LMData.DiffPercentRT = RTEffect_nonorm;
    
    LMData.MeanPVInRF{1,1} = rfindatapre.AveragePV(find(rfindatapre.RFx == xx & rfindatapre.RFy == yy));
    LMData.MeanPVInRF{2,1} = rfindatapost.AveragePV(find(rfindatapost.RFx == xx & rfindatapost.RFy == yy));
    LMData.MedianPVInRF{1,1} = rfindatapre.MedianPV(find(rfindatapre.RFx == xx & rfindatapre.RFy == yy));
    LMData.MedianPVInRF{2,1} = rfindatapost.MedianPV(find(rfindatapost.RFx == xx & rfindatapost.RFy == yy));
    
    LMData.MeanPVOutRF{1,1} = rfindatapre.AveragePV(find(rfindatapre.RFx == -xx & rfindatapre.RFy == yy));
    LMData.MeanPVOutRF{2,1} = rfindatapost.AveragePV(find(rfindatapost.RFx == -xx & rfindatapost.RFy == yy));
    LMData.MedianPVOutRF{1,1} = rfindatapre.MedianPV(find(rfindatapre.RFx == -xx & rfindatapre.RFy == yy));
    LMData.MedianPVOutRF{2,1} = rfindatapost.MedianPV(find(rfindatapost.RFx == -xx & rfindatapost.RFy == yy));
    
    LMData.MeanRTInRF{1,1} = rfindatapre.AverageRT(find(rfindatapre.RFx == xx & rfindatapre.RFy == yy));
    LMData.MeanRTInRF{2,1} = rfindatapost.AverageRT(find(rfindatapost.RFx == xx & rfindatapost.RFy == yy));
    
    LMData.MeanRTOutRF{1,1} = rfindatapre.AverageRT(find(rfindatapre.RFx == -xx & rfindatapre.RFy == yy));
    LMData.MeanRTOutRF{2,1} = rfindatapost.AverageRT(find(rfindatapost.RFx == -xx & rfindatapost.RFy == yy));
    
    LMData.Trial_3{1,1} = Trial_3_pre;
    LMData.Trial_3{2,1} = Trial_3_post;
    LMData.RFx = xx;
    LMData.RFy = yy;
end
%% Saving the data

if input('Save LM Data? y/n', 's') == 'y'
    savingFilename = input('Name the saved variables with initial of monkey, date, and "LMData": ', 's'); % Name of file
    savingPath = [pwd]; % Location to save the file in
    save([savingPath '/' savingFilename],'LMData'); % Save the file
else
end
