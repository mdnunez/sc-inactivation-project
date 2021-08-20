clc;
clear all;
close all;

%% Load datafile 
%Getting filename and datafile directory and loading NSx (Trial) and NEV files (Trial_1)
[fileName, pathname] = getFile('*.nev');
filename = fileName(1,1:end-4);

file_NSx = [filename '.ns5'];

file_in = strcat(pathname,file_NSx);
fileNEV_in = strcat(pathname,fileName);

Trial = openNSx(file_in);
Trial_1 = openNEV(fileNEV_in);

%% Create Data Structure (Trial_3)
Trial_3 = struct('eCodes', [], 'correct', [],'timeindexes', [], 'times', [], 'RT', [], 'RT_1503', [], 'PeakVelocity', [], 'RFx', [], 'RFy', [], 'first_saccade_time', [], 'first_saccade_index', []);
% Retrive Event Code Structure Array from unpacked file_in
starttrial_indices = find(Trial_1.Data.SerialDigitalIO.UnparsedData == 1001);

for s = 1:numel(starttrial_indices)
    if s ~= numel(starttrial_indices)
        Trial_3(s).eCodes = double([Trial_1.Data.SerialDigitalIO.UnparsedData(starttrial_indices(s):starttrial_indices(s+1)-1)]);
        Trial_3(s).timeindexes = double([Trial_1.Data.SerialDigitalIO.TimeStamp(starttrial_indices(s):starttrial_indices(s+1)-1)]);
        Trial_3(s).times = double([Trial_1.Data.SerialDigitalIO.TimeStampSec(starttrial_indices(s):starttrial_indices(s+1)-1)]);
    else
        Trial_3(s).eCodes = double([Trial_1.Data.SerialDigitalIO.UnparsedData(starttrial_indices(s):end)]);
        Trial_3(s).timeindexes = double([Trial_1.Data.SerialDigitalIO.TimeStamp(starttrial_indices(s):end)]);
        Trial_3(s).times = double([Trial_1.Data.SerialDigitalIO.TimeStampSec(starttrial_indices(s):end)]);
    end
end
remove = [];    %create a trash bin

%%
%Filter out invalid, incorrect trials or trials with no RF codes
for i=1:numel(Trial_3)     %go through all trials
    
    codes = Trial_3(i).eCodes; % Event Codes for the current trial
    [indx_RFx] = find(codes >= 6200 & codes <= 6800,1); %find the indices of the RF ecodes
    [indx_RFy] = find(codes >= 4200 & codes <= 4800,1);
    idx_3 = [indx_RFx, indx_RFy];
    
    if size(idx_3,2) ~= 2    %if the number of columns (number of values) is less than 2 for the indices for RF's
        remove = [remove, i]; % Add to remove list
    elseif size(find(codes == 1001),2) ~= 1    %make sure it has exactly 1 start 1001 code
        remove = [remove, i];
    elseif size(find(codes == 3000),2) ~= 1      %make sure it has exactly 1 fpoff 3000 code
        remove = [remove, i];
    elseif size(find(codes == 5050),2) ~= 1     %make sure it has exactly  1 reward 5050 code
        remove = [remove, i];
    else
        Trial_3(i).correct = 1;         %if it passes through these filters, then label as correct in the data structure
    end
end

clear x y
Trial_3(remove) = [];       %now Trial_3 has a list of all the valid, correct trials

%% detect when sac was made, get peak velocity, and update Trial_3 struc
[Trial_4,allh,allv,saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold, accelerationThreshold] = Saccade_Detector_VGST(Trial, Trial_3);

N = numel(Trial_4); %new number of elements in filtered trial structure from fn_pos_vel_accel_post_sac_EJ

RFxyM = [];
RFx = [];
RFy = [];
remove_i = [];
%Getting all the RF's and RT's (from ecode and SD) in the structure
for i=1:N   %for all the trials in valid trial struc
    
    codes = Trial_4(i).eCodes;  %get current trial codes
    Times = Trial_4(i).times;   %get current trial times for the codes
    
    idx = find(codes == 3000);    %get index for FPOFF ecode
    idx_2 = find(codes == 1503); %get index for SACCD ecode
    idx_RFx = find(codes >= 6200 & codes <= 6800,1,'first'); %get indices of RF ecodes
    idx_RFy = find(codes >= 4200 & codes <= 4800,1,'first');
    
    % RF coordinates
    if ~isempty(codes(idx_RFx)) && ~isempty(codes(idx_RFy)) && ~isempty(codes(idx))
            RFxy1(1,1) = codes(idx_RFx)-6500;  %x position
            RFxy1(1,2) = codes(idx_RFy)-4500;  %y position
        Trial_4(i).RFx = RFxy1(1,1);
        Trial_4(i).RFy = RFxy1(1,2);
        
        % Fixation point off and saccade times
        Trial_4(i).FpOffTime = Times(1,idx);  %list of all the FPOFF abosolute times in sec
        Trial_4(i).SACCDTime = Times(1,idx_2); %List of all the SACCD absolute times in sec -- sac inititation time derived from ecodes, will not use this for analysis
        
        % Storing times in Trial structure
        Trial_4(i).RT_1503 = Trial_4(i).SACCDTime -  Trial_4(i).FpOffTime;  %RT_1503 is RT derived from event codes -- sac inititation time derived from ecodes, will not use this for analysis
        Trial_4(i).RT =  Trial_4(i).first_saccade_time -  Trial_4(i).FpOffTime;   %RT is RT derived from saccade detector saccade times
        
        %filter out trials when there's no FP off time
        if Trial_4(i).FpOffTime == 0
            remove_i = [remove_i, i];
        end
        %making a list of RF's
        RFx = [RFx; Trial_4(i).RFx];
        RFy = [RFy; Trial_4(i).RFy];
    else
        remove_i = [remove_i, i];
    end
    
end
%put into scale for the plot for peak velocity
RFx = RFx + 251;
RFy = abs(RFy - 251);
RFxyM = [RFx, RFy];

Trial_4(remove_i) = [];
%%
%Cleaning up data -- to manually double check the saccades if you want and
%get rid of trials in which the monkey made weird saccades 
N = numel(Trial_4);
remove_3 = [];

[Trial_5, remove_3] = SacDetectManualCheck_VGST(Trial, Trial_4, 1, 0, saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold,accelerationThreshold,allh, allv);

Trial_5(unique(remove_3)) = [];

%Save this session's data for further processing of aggreggate data later 
input_save = 'z';
while input_save ~= 'y' && input_save ~= 'n'
    input_save = input('Save Trials Structure? y/n', 's');
    if input_save == 'y'
        savingFilename = [filename '_LMdata']; % Name of file
        savingPath = [pwd]; % Location to save the file in
        save([savingPath '/' savingFilename], 'Trials'); % Save the file
    else
    end
end

