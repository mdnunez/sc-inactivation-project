function [saccAwayRF, saccToRF, Trial_1, TrialTargInRF,TrialTargOutRF,PropTrialInRF,PropTrialOutRF, percCin, percCopp,filename] = selectiontask_rawdataprocessor(a,b,files_inputted,NEV, REX)

if NEV == 1 %unpack NEV file if data was recorded with NEV
    [fileName, pathname] = getFile('*.nev');
    filename = fileName(1,1:end-4);
    fileNEV_in = strcat(pathname,fileName);
    [Trial] = openNEV(fileNEV_in);
    
    file_NSx = [filename '.ns5'];
    file_in = strcat(pathname,file_NSx);
    Trial_nsx = openNSx(file_in, 'report');
    
    Trial_1 = struct('eCodes', [], 'timeindexes', [], 'times', [], 'RT', [], 'PeakVelocity', [], 'first_saccade_time', [], 'first_saccade_index', []);
    
    % Retrive Event Code Structure Array from unpacked file 
    starttrial_indices = find(Trial.Data.SerialDigitalIO.UnparsedData == 1001);
    
    for s = 1:numel(starttrial_indices)
        if s ~= numel(starttrial_indices)
            Trial_1(s).eCodes = double([Trial.Data.SerialDigitalIO.UnparsedData(starttrial_indices(s):starttrial_indices(s+1)-1)]);
            Trial_1(s).timeindexes = double([Trial.Data.SerialDigitalIO.TimeStamp(starttrial_indices(s):starttrial_indices(s+1)-1)]);
            Trial_1(s).times = double([Trial.Data.SerialDigitalIO.TimeStampSec(starttrial_indices(s):starttrial_indices(s+1)-1)]);
        else
            Trial_1(s).eCodes = double([Trial.Data.SerialDigitalIO.UnparsedData(starttrial_indices(s):end)]);
            Trial_1(s).timeindexes = double([Trial.Data.SerialDigitalIO.TimeStamp(starttrial_indices(s):end)]);
            Trial_1(s).times = double([Trial.Data.SerialDigitalIO.TimeStampSec(starttrial_indices(s):end)]);
        end
    end
    
elseif REX == 1 %unpack REX file if data was recorded with REX
    [fileName pathname] = getFile('*.*');
    filename = fileName(1,1:end-4);
    file_in = strcat(pathname,fileName);
    Trial = rex_unpack(file_in);
    % Create Data Structure (Trial_3)
    Trial_1 = struct('eCodes', [], 'timeindexes', [], 'times', [], 'RT', [], 'PeakVelocity', [], 'first_saccade_time', [], 'first_saccade_index', []);
    for i = 1:numel(Trial)
        Trial_1(i).eCodes = Trial(i).codes;
        Trial_1(i).timeindexes = Trial(i).t;
        Trial_1(i).times = Trial(i).t./2000;
    end
end

rf = [a,b]; %recpetive field 
m = numel(Trial_1); 
remove = [];
%see whether target is inRF or outRF and whether trial was correct or not 
for i=1:m
    if ~isempty(find(Trial_1(i).eCodes>=4200 & Trial_1(i).eCodes<=4800, 1))
        tg = find(Trial_1(i).eCodes>=4200 & Trial_1(i).eCodes<=4800); %ecodes for the choice target locations in x and y 
        if Trial_1(i).eCodes(tg(1)) == 4500 + rf(1)
            Trial_1(i).TGinRF = 1; %target is inRF
        elseif Trial_1(i).eCodes(tg(1)) == 4500 - rf(1)
            Trial_1(i).TGinRF = 0; %target is outRF 
        else
            Trial_1(i).TGinRF = NaN;
            remove = [remove, i];
        end
    end
    
    if ~isempty(find(Trial_1(i).eCodes ==5000))
        Trial_1(i).outcome=1;   %correct
    elseif ~isempty(find(Trial_1(i).eCodes==5006))
        Trial_1(i).outcome=0;   %incorrect
    end
    
end

NumTrials = numel(Trial_1);
remove_1 = [];
%extracting out event times 
for i = 1:NumTrials
    codes = Trial_1(i).eCodes;  %get current trial codes
    Times = Trial_1(i).times;   %get current trial times for the codes
    idx_RFx = find(codes >= 7200 & codes <= 7800); %get indices of RF ecodes
    idx_RFy = find(codes >= 8200 & codes <= 8800);
    idx_FPOFF = find(codes == 3000, 1, 'first'); %get index for FPOFF ecode
    idx_TGON = find(codes == 2000, 1, 'first');
    % Fixation point off and saccade times
    Trial_1(i).FPOff = Times(1,idx_FPOFF);
    Trial_1(i).TGOn = Times(1,idx_TGON);
    if isempty(Trial_1(i).FPOff)
        Trial_1(i).FPOff = nan;
    end
    if isempty(Trial_1(i).TGOn)
        Trial_1(i).TGON = nan;
    end
    Trial_1(i).DelayPeriod = Trial_1(i).FPOff - Trial_1(i).TGOn;
    
    if Trial_1(i).DelayPeriod < 0.4  %remove the trial if the delay period was less the 400ms
        remove_1 = [remove_1,i];
    end
end

if ~isempty(remove_1)
    Trial_1(unique(remove_1)) = [];
end
%% Saccade detector
NumTrials = numel(Trial_1);

[Trial_1,allh,allv,saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold, accelerationThreshold] = Saccade_Detector_STO(Trial_nsx, Trial_1, files_inputted);
%%
%Getting all the RT's (from ecode and saccade detector) in the structure
remove = [];
N = numel(Trial_1);
for i=1:N   %for all the trials in valid trial struc
    
    codes = Trial_1(i).eCodes;  %get current trial codes
    Times = Trial_1(i).times;   %get current trial times for the codes
    
    idx_RFx = find(codes >= 7200 & codes <= 7800); %get indices of RF ecodes
    idx_RFy = find(codes >= 8200 & codes <= 8800);
    idx_FPOFF = find(codes == 3000, 1, 'first'); %get index for FPOFF ecode
    
    % Fixation point off and saccade times
    Trial_1(i).FPOff = Times(1,idx_FPOFF);
    if isempty(Trial_1(i).FPOff)
        Trial_1(i).FPOff = nan;
    end
    
    %Storing times in Trial structure
        Trial_1(i).RT = Trial_1(i).first_saccade_time - Trial_1(i).FPOff;
end

if ~isempty(remove)
    Trial_1(unique(remove)) = [];
end
%% Manually check saccades if you want 
remove_3 = [];
[Trial_2, remove_3] = SacDetectManualCheck_STO(Trial_nsx, Trial_1, saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold,accelerationThreshold,allh, allv);

remove_3 = unique(remove_3);
for r = 1:numel(remove_3)
    Trial_2(remove_3(r)).PeakVelocity = nan;
    Trial_2(remove_3(r)).RT = nan;
end
numel(remove_3)
%%
NumTrials = numel(Trial_2);
m = NumTrials;
TinC=NaN(1,m); %correct trial when target is inRF
ToppC=NaN(1,m); %correct trial when target is outRF
TinDS=NaN(1,m); %incorrect trial when target is inRF
ToppDS=NaN(1,m); %incorrect trial when target is outRF

for i=1:m
    if Trial_2(i).TGinRF == 1 &  Trial_2(i).outcome==1     %got it correct, target is in RF, chose inRF
        TinC(i)=(1);
    end
    if Trial_2(i).TGinRF == 0 &  Trial_2(i).outcome==1     %got it correct, target is out RF, chose outRF
        ToppC(i)=(1);
    end
    if Trial_2(i).TGinRF == 1 &  Trial_2(i).outcome==0     %got it incorrect, target is inRF, chose outRF
        TinDS(i)=1;
    end
    if Trial_2(i).TGinRF == 0 &  Trial_2(i).outcome==0     %got it incorrect, target is outRF, chose inRF
        ToppDS(i)=1;
    end
end

nTinC= length(find(TinC==1));
nToppC= length(find(ToppC==1));
nTinDS= length(find(TinDS==1));
nToppDS= length(find(ToppDS==1));

percCin = nTinC/(nTinC+nTinDS); %accuracy toRF (inRF)
percCopp = nToppC/(nToppC+nToppDS); %accuracy awayRF (outRF)

saccToRF = (nTinC+nToppDS)/(nTinC+nToppDS+nTinDS+nToppC); %percent of saccades toRF
saccAwayRF = (nTinDS+nToppC)/(nTinC+nToppDS+nTinDS+nToppC); %percent of saccades awayRF

%Getting info for target trial distribution
TrialTargInRF = 0;
TrialTargOutRF = 0;
for i = 1:m
    x = double(Trial_2(i).eCodes(find(Trial_2(i).eCodes>=4200 & Trial_2(i).eCodes<=4800, 1, 'first'))) - 4500;
    if x == a
        TrialTargInRF = TrialTargInRF + 1;
    elseif x == -a
        TrialTargOutRF = TrialTargOutRF + 1;
    end
end

PropTrialInRF = TrialTargInRF/(TrialTargInRF+TrialTargOutRF);
PropTrialOutRF = TrialTargOutRF/(TrialTargInRF+TrialTargOutRF);

%%
%Cleaning out Trial_4 to save
takeoutfields = {'eCodes', 'timeindexes', 'times', 'first_saccade_time', 'first_saccade_index', 'all_vertical_positions_sized'...
    'all_horizontal_positions_sized', 'alltimes_index_sized', 'velocity_vertical', 'velocity_horizontal', 'velocity_diagonal', 'accelerations',...
    'FPOff', 'TGOn', 'DelayPeriod'};
Trial_2 = rmfield(Trial_2, takeoutfields);
end