function [Trial_4,coh,filename,y_avg_med_RT, RTStruct,YAndParams] = decisiontask_rawdataprocessor(coh_choice,files_inputted, RSC, LSC, RT, DT);
%-------------------- Set Coherences in Descending Order --------------------%

% Set the appropriate coherences based on the spot file used for the experiment conducted, in descending order
% Use the examples below for reference
switch coh_choice  
    case 1
        coh = [36 24 17 10 5 3 0]; % coherences for SC muscimol project for 3 muscimol injections 
    case 2
        coh = [50 36 24 17 10 5 0]; % coherences for SC muscimol project for the rest of the injections
end

continueAction_concat = 1;
files_concat = 1;
DataStruct = struct();
Trial_4 = [];

%%

NEV = 1; % 1 if NEV file was recorded; 0 if NEV file not recorded and REX file was recorded instead
REX = 0; % 1 if REX file was recorded; 0 if REX file was not recorded and NEV was recorded instead
while(continueAction_concat)
    
    if NEV == 1
        % Pass the raw data through the openNEV function
        [fileName, pathname] = getFile('*.nev');
        filename = fileName(1,1:end-4);
        fileNEV_in = strcat(pathname,fileName);
        [Trial_1] = openNEV(fileNEV_in);
        
        file_NSx = [filename '.ns5'];
        file_in = strcat(pathname,file_NSx);

        Trial = openNSx(file_in, 'report', 'c:2:3'); %c:3:4 for BT 12/10 injection for pre/post; c:1:2 for some of sparkys data, c:17:18 for some of sparky's data TT recorded, c:33:34 for 11/19 monkey B data, 
        
        % Create Data Structure (Trial_3)
        Trial_3 = struct('eCodes', [], 'correct',[], 'difficulty', [], 'coherence', [], 'target', [], 'VRNum', [],'timeindexes', [], 'times', [], 'RT', [], 'PeakVelocity', [], 'RFx', [], 'RFy', [], 'first_saccade_time', [], 'first_saccade_index', []);
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
        
    elseif REX == 1
        [fileName pathname] = getFile('*.*');
        filename = fileName(1,1:end-4);
        file_in = strcat(pathname,fileName);
        Trial = rex_unpack(file_in);
        % Create Data Structure (Trial_3)
        Trial_3 = struct('eCodes', [], 'correct',[], 'difficulty', [], 'coherence', [], 'target', [], 'VRNum', [],'timeindexes', [], 'times', [], 'RT', [], 'PeakVelocity', [], 'RFx', [], 'RFy', [], 'first_saccade_time', [], 'first_saccade_index', []);
        for i = 1:numel(Trial)
            Trial_3(i).eCodes = Trial(i).codes;
            Trial_3(i).timeindexes = Trial(i).t;
            Trial_3(i).times = Trial(i).t./2000;
        end 
    end
    

%%
% For-Loop for the identification of correct Trial_3s and removal of incorrect Trial_3s
remove = [];
trialtypes = struct('correct', 0, 'incorrect', 0, 'failtoaqfp', 0, 'anticipsac', 0, 'failedsactime', 0,'failtohldfp',0,'failhldtg', 0,...
    'violatedminLT', 0);

for i=1:numel(Trial_3)
    
    codes = Trial_3(i).eCodes; % Event Codes
    
    if ~isempty(find(codes == 5510)) % Saccade to correct target
        trialtypes.correct = trialtypes.correct + 1;
        Trial_3(i).correct = 1;
    elseif ~isempty(find(codes == 5006)) % Saccade to distractor
        Trial_3(i).correct = 0;
    elseif ~isempty(find(codes == 5001)) %Failed to acquire fix point
        trialtypes.failtoaqfp = trialtypes.failtoaqfp +1;
        remove = [remove i];
    elseif ~isempty(find(codes == 5004)) %Failed to hold target
        trialtypes.failtohldfp = trialtypes.failtohldfp +1;
        remove = [remove i];
    elseif ~isempty(find(codes == 5005)) %Anticipatory accade
        trialtypes.anticipsac = trialtypes.anticipsac + 1;
        remove = [remove i];
    elseif ~isempty(find(codes == 5007)) %Failed to make a saccade in time
        trialtypes.failedsactime = trialtypes.failedsactime + 1;
        remove = [remove i];
    elseif ~isempty(find(codes == 5009))  %fail to hold fp
        trialtypes.failtohldfp = trialtypes.failtohldfp + 1;
        remove = [remove i];
    elseif ~isempty(find(codes == 5010))   %violated 50ms min LT
        trialtypes.violatedminLT = trialtypes.violatedminLT + 1;
        remove = [remove i];
    else
        remove = [remove i];
    end
    
    %filter out trials where crucial ecodes didnt go through
    if isempty(find(codes == 3000)) || isempty(find(codes == 5000)) || numel(find(codes >= 4100 & codes <= 4100+numel(coh)))>1
        remove = [remove i];
    end
    
end

%remove invalid trials
Trial_3(unique(remove)) = [];  

% Determine time stamps for several stimuli
N = numel(Trial_3);
for i = 1:N
    
    % Event Codes
    codes = Trial_3(i).eCodes; 
    
    % Difficulty levels w.r.t coherence
    idx = find(codes >= 4100 & codes <= 4100+numel(coh));
    diff_code = codes(idx);
    diff_code = diff_code - 4100 + 1;
    Trial_3(i).difficulty = diff_code;
    
    % Coherence values
    coh_val = coh(diff_code);
    Trial_3(i).coherence = coh_val;
    
    % target code
    idx_1 = codes >= 4000 & codes <= 4001;
    tar_code = codes(idx_1);
    tar_code = tar_code - 4000 + 1;
    Trial_3(i).target = tar_code;     %1 is Left, 2 is Right 
    
    %VRnum -- variable ratio number for the trial
    idx_2 = find(codes >= 6000 & codes < 6100);
    if isnan(idx_2)     %if there are no codes for VRnum, then vrnum code is nan
        vrnum_code = nan;
    else
        vrnum_code = codes(idx_2);
    end 
    Trial_3(i).VRNum = vrnum_code-6000;
end


%%

NumTrials = numel(Trial_3);
[Trial_3,allh,allv,saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold, accelerationThreshold] = Saccade_Detector_GP_RTtask(Trial, Trial_3, RT, DT);
%%
%Getting all the RT's in the structure
outofrangelowRT = [];
outofrangehighPV = [];
outofrangelowPV = [];
remove = [];
N = numel(Trial_3);
for i=1:N   %for all the trials in valid trial struc
    
    codes = Trial_3(i).eCodes;  %get current trial codes
    Times = Trial_3(i).times;   %get current trial times for the codes
    
    idx = find(codes == 5000);    %get index for GPON ecode
    idx_RFx = find(codes >= 7200 & codes <= 7800); %get indices of RF ecodes
    idx_RFy = find(codes >= 8200 & codes <= 8800);
    idx_FPOFF = find(codes == 3000); %get index for FPOFF ecode    
    
    % GPOn and FPOff and saccade times
    Trial_3(i).GPOn = Times(1,idx);  %list of all the GPOn abosolute times in sec
    Trial_3(i).FPOff = Times(1,idx_FPOFF);
   
    
    % Storing times in Trial structure 
    Trial_3(i).RT =  Trial_3(i).first_saccade_time -  Trial_3(i).GPOn;   %RT is RT derived from saccade detector saccade times

    
    %if there was somehow a glitch in the ecodes and 2 or none GPOn ecodes
    %were recorded with 2 different timestamps, delete the trial
    if numel(Trial_3(i).RT) > 1 || Trial_3(i).GPOn == 0 
        remove = [remove, i];
    end 
    
    Trial_3(i).RFx = Trial_3(i).eCodes((idx_RFx))-7500;
    Trial_3(i).RFy = Trial_3(i).eCodes(idx_RFy)-8500;
    
end

Trial_3(remove) = [];
N = numel(Trial_3);

%%
% Cleaning up data -- invalid data points

LowRT = 0.15;

[Trial_3, remove_3] = SacDetectManualCheck_GP(Trial, Trial_3,DT, RT, saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold,accelerationThreshold,allh, allv);
stop = 1;
Trial_3(remove_3) = []; 
%print out the min, max of RT's and peak velocities just to make sure that
%they are in range and the saccade detector worked well
min([Trial_3.RT])
max([Trial_3.RT])
min([Trial_3.PeakVelocity])
max([Trial_3.PeakVelocity])
numel(remove_3) %print out of how many trials were removed because a saccade wasn't able to be detected well 
stop =1 ;
display('Last concatenated file was')
filename
prompt_0 = ('Concatenate with more GP files? (y/n)');   %concatenate more data into the same session?
    add_file = input(prompt_0, 's');
    if add_file == 'Y' || add_file == 'y'
        continueAction_concat = 1;
        files_concat = files_concat + 1;
    elseif add_file == 'N' || add_file == 'n'
        continueAction_concat = 0;
    end
    
Datastruct.(['Trial_3_' num2str(files_concat)]) = Trial_3;
Trial_4 = [Trial_4, Trial_3];

end 

        
%% Getting the parameters for psychometric function
N = numel(Trial_4);
% 2 by 7 correct and incorrect Trial_4s matrices. Raws indicate the target and columns the level of coherence
C=zeros(2,length(coh)); %% Correct Trial_4s
I=zeros(2,length(coh)); %% Incorrect Trial_4s

for i = 1:N
    
    if Trial_4(i).correct == 1
        C(Trial_4(i).target,Trial_4(i).difficulty)=C(Trial_4(i).target,Trial_4(i).difficulty)+1;
    elseif Trial_4(i).correct == 0
        I(Trial_4(i).target,Trial_4(i).difficulty)=I(Trial_4(i).target,Trial_4(i).difficulty)+1;
    end
    
end

% 2 by 7 matrix sum of all the trials, both incorrect and correct
A = C+I;

% Ynum trials chosen to the right
y_num_neg = A(1,:) - C(1,:);
y_num_pos = fliplr(C(2,:));
y_num_R = [y_num_neg y_num_pos];

%Ynum trials chosen correct (Acc) from coh0 to coh36
for i = 1:size(C, 2)
    y_num_acc(i) = sum(C(:,i)); 
    y_num_ntrial_acc(i) = sum(A(:,i));
    y_prop_acc(i) = y_num_acc(i)./y_num_ntrial_acc(i);
end 
% flip them so that it's hardes coh --> easiest coh to plot the accuracy
% PMF
y_num_acc = fliplr(y_num_acc);
y_num_ntrial_acc = fliplr(y_num_ntrial_acc);
y_prop_acc = fliplr(y_prop_acc);

% percentage of accuracy by coherence
acc = C./A;

x_pos = fliplr(coh);            % hardest to easiest for positive coh
y_pos = fliplr(acc(2,:));       % hardest to easiest for positive coh
trialn_pos = fliplr(A(2,:));    % total trials hardest to easiest for positive coh
x_neg = sort(-coh);
y_neg = 1-acc(1,:);             %to convert to choosing for the right, rather than accuracy
trialn_neg = A(1,:);            %total trials
x_both_R = [x_neg x_pos];
y_both = [y_neg y_pos];
trialn_both = [trialn_neg, trialn_pos];
x_both_Acc = fliplr(coh);

N=length(coh);

if sum(double(ismember(coh,0))) > 0         %if there's a 0 coherence in the coherence values
    
    x_both_R=[x_both_R(1:N-1),(x_both_R(N)+x_both_R(N+1)),x_both_R(N+2:end)];
    y_both=[y_both(1:N-1),(y_both(N)+y_both(N+1))/2,y_both(N+2:end)];       %proportion of trials chosen to right
    y_num_R = [y_num_R(1:N-1),(y_num_R(N)+y_num_R(N+1)),y_num_R(N+2:end)];            %number trials chosen to right
    trialn_both = [trialn_both(1:N-1),(trialn_both(N)+trialn_both(N+1)),trialn_both(N+2:end)];
    
end

y_both = y_both.*100;       % for the R vs L PMF (psychometric function)
y_prop_acc = y_prop_acc.*100;       %for the Acc PMF

%Convert data to proportion toIF choices PMF, where positive coherences are
%toIF coherences and negative coherences are awayIF coherences
if RSC == 1
    y_num_inact = fliplr(trialn_both-y_num_R);
    y_ntrials_inact = fliplr(trialn_both);
elseif LSC == 1
    y_num_inact = y_num_R;
    y_ntrials_inact = trialn_both;
end 
YPropInact = 100*(y_num_inact./y_ntrials_inact);

%%
% -------------------- Fit The Psychometric Data with least squares method --not in use ----------------------%

% % initial parameter estimates
% params_init  = [50 50 0.01 0];
% options = optimset('display', 'off');
% 
% % unconstrained fit
% params_fit = fminsearch(@fit_err,params_init,options,x_both,y_both);
% 
% % fitted curve
% x_fit = linspace(x_both(1),x_both(end),100);
% y_fit = sigm(params_fit,x_fit);

% -------------------- Fit The Psychometric Data with maximum likelihood method ----------------------%
%Now using palamedes toolbox to get psychometric function parameters, alpha
%and beta, using MLE
%Gamma is LambdaL and lambda is LambdaR
if ~exist('PAL_info','file')
    disp('This snippet requires the Palamedes toolbox.');
    disp('Visit www.palamedestoolbox.org');
    return;
end
options = PAL_minimize('options');

%Case of field names must match exactly. E.g., options.tolx = 1 e-12 will 
%be ignored (no warning will be issued). Options structure must be passed
%to PAL_PFML_Fit (see call to PAL_PFML_Fit below) or will go ignored.
options.TolX = 1e-50;   %increase desired precision (by lowering tolerance)
options.TolFun = 1e-50;
options.MaxIter = 10000;
options.MaxFunEvals = 10000;
options.Display = 'iter';   %Display progress of Nelder-mead iteration-by-
                            %iterationn
%Use the Logistic function
PF = @PAL_Logistic;

%Threshold (alpha) and Slope (beta) are free parameters, guess (gamma) and
%lapse (lambda) rate are fixed, in that order of input for 
paramsFreeR = [1 1 0 0];  %1: free parameter, 0: fixed parameter
paramsFreeAcc = [1 1 0 1];  %guess rate is always fixed to 0.5 in Acc PMFs

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGridR.alpha = -20:1:20;
searchGridR.beta = logspace(0,3,101);
searchGridR.gamma =  0:0.01:0.05; %LambdaL; scaar here (since fixed) but may be vectorl
searchGridR.lambda = 0:0.01:0.05;  %LambdaR
%YAcc 
searchGridAcc.alpha = 0:1:36;
searchGridAcc.beta = logspace(0,3,101); 
searchGridAcc.gamma = 0.5;  
searchGridAcc.lambda = 0:0.01:0.2;  


%Perform fit for Y prop R (PMF for proportion of choices to the right)
disp('Fitting function.....');
[paramsValuesYRight LL exitflag] = PAL_PFML_Fit(x_both_R,y_num_R, ...
    trialn_both,searchGridR,paramsFreeR,PF)

%Perform fit for Y Inact (PMF for the proportion of choices toIF)
disp('Fitting function.....');
[paramsValuesInact LL exitflag] = PAL_PFML_Fit(x_both_R,y_num_inact, ...
    y_ntrials_inact,searchGridR,paramsFreeR,PF)

%Perform fit for Y Acc
disp('Fitting function.....');
[paramsValuesYAcc LL exitflag] = PAL_PFML_Fit(x_both_Acc,y_num_acc, ...
    y_num_ntrial_acc,searchGridAcc,paramsFreeAcc,PF)

x_fit = linspace(x_both_R(1),x_both_R(end),100);                            %for both RvsL and Inact PMF
x_fitAcc = linspace(x_both_Acc(1), x_both_Acc(end), 100);                   %for Acc PMF

y_HatML_YRight = PAL_Logistic(paramsValuesYRight, x_fit);
y_HatML_YRight = y_HatML_YRight*100;

Y_HatML_Inact = PAL_Logistic(paramsValuesInact, x_fit);
Y_HatML_Inact = Y_HatML_Inact*100;

Y_HatML_YAcc = PAL_Logistic(paramsValuesYAcc, x_fitAcc);
Y_HatML_YAcc = Y_HatML_YAcc*100;


YAndParams = struct('DataRParams', [paramsValuesYRight], 'DataInactParams', paramsValuesInact, 'DataAccParams', [paramsValuesYAcc], 'DataYPropR', [y_both],...
    'DataYNumR', [y_num_R], 'DataYnTrialR', [trialn_both], 'DataYPropAcc', [y_prop_acc],...
    'DataYNumAcc', [y_num_acc], 'DataYnTrialAcc', [y_num_ntrial_acc], 'DataYPropInact', [YPropInact], 'DataYNumInact', [y_num_inact],...
    'DataYnTrialInact', [y_ntrials_inact], 'x_coh', [x_both_R], 'x_coh_Acc', [x_both_Acc], 'xfit', x_fit, 'xfitacc', x_fitAcc, 'Y_HatML_YRight',y_HatML_YRight, ...
    'Y_HatML_Inact', [Y_HatML_Inact], 'Y_HatML_YAcc', [Y_HatML_YAcc]);


%%
%1 is left, 2 is right, difficulty 1-7 is from easiest to hardest


%A structure that holds all the RT with coh = columns [-36, -24, -17, -10,
%-5, -3, 0, 3, 5, 10, 17, 24, 36]. 
%1st row of RTStruct is all the trials (correct + incorrect) together
RTStruct(1,1)= {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 1)).RT]};
RTStruct(1,2) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 1)).RT]};
RTStruct(1,3) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 1)).RT]};
RTStruct(1,4) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 1)).RT]};
RTStruct(1,5) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 1)).RT]};
RTStruct(1,6) = {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 1)).RT]};
%if there's a 0 coherence, lump all the 0 coherences from left and right
%together
if coh(end) == 0
    RTStruct(1,7) = {[[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2)).RT], [Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1)).RT]]};
    RTStruct(1,8)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,9) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,10) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,11) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,12) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,13) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2)).RT]};
    %Stats on the RTs: 1st row is the stats on all the trials together; 1st
    %column is mean, 2nd column is median
    y_avg_med_RT{1,1} = [mean(RTStruct{1,1}), mean(RTStruct{1,2}), mean(RTStruct{1,3}), mean(RTStruct{1,4}), mean(RTStruct{1,5}),...
    mean(RTStruct{1,6}), mean(RTStruct{1,7}), mean(RTStruct{1,8}), ... 
    mean(RTStruct{1,9}), mean(RTStruct{1,10}), mean(RTStruct{1,11}), mean(RTStruct{1,12}), mean(RTStruct{1,13})];   
    y_avg_med_RT{1,2} = [median(RTStruct{1,1}), median(RTStruct{1,2}), median(RTStruct{1,3}), median(RTStruct{1,4}), median(RTStruct{1,5}),...
    median(RTStruct{1,6}), median(RTStruct{1,7}), median(RTStruct{1,8}), ... 
    median(RTStruct{1,9}), median(RTStruct{1,10}), median(RTStruct{1,11}), median(RTStruct{1,12}), median(RTStruct{1,13})]; 
else
    RTStruct(1,7) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1)).RT]};
    RTStruct(1,8) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,9)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,10) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,11) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,12) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,13) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2)).RT]};
    RTStruct(1,14) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2)).RT]};
    y_avg_med_RT{1,1} = [mean(RTStruct{1,1}), mean(RTStruct{1,2}), mean(RTStruct{1,3}), mean(RTStruct{1,4}), mean(RTStruct{1,5}),...
        mean(RTStruct{1,6}), mean(RTStruct{1,7}), mean(RTStruct{1,8}), ...
        mean(RTStruct{1,9}), mean(RTStruct{1,10}), mean(RTStruct{1,11}), mean(RTStruct{1,12}), mean(RTStruct{1,13}), mean(RTStruct{1,14})];
    y_avg_med_RT{1,2} = [median(RTStruct{1,1}), median(RTStruct{1,2}), median(RTStruct{1,3}), median(RTStruct{1,4}), median(RTStruct{1,5}),...
        median(RTStruct{1,6}), median(RTStruct{1,7}), median(RTStruct{1,8}), ...
        median(RTStruct{1,9}), median(RTStruct{1,10}), median(RTStruct{1,11}), median(RTStruct{1,12}), median(RTStruct{1,13}), median(RTStruct{1,14})];
end

% For RTs of correct trials
% 2nd row of RT is the RTs of the correct trials 
RTStruct(2,1)= {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
RTStruct(2,2) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
RTStruct(2,3) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
RTStruct(2,4) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
RTStruct(2,5) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
RTStruct(2,6) = {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
if coh(end) == 0
    RTStruct(2,7) = {[[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT], [Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]]};
    RTStruct(2,8)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,9) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,10) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,11) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,12) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,13) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    y_avg_med_RT{2,1} = [mean(RTStruct{2,1}), mean(RTStruct{2,2}), mean(RTStruct{2,3}), mean(RTStruct{2,4}), mean(RTStruct{2,5}),...
    mean(RTStruct{2,6}), mean(RTStruct{2,7}), mean(RTStruct{2,8}), ... 
    mean(RTStruct{2,9}), mean(RTStruct{2,10}), mean(RTStruct{2,11}), mean(RTStruct{2,12}), mean(RTStruct{2,13})];   
    y_avg_med_RT{2,2} = [median(RTStruct{2,1}), median(RTStruct{2,2}), median(RTStruct{2,3}), median(RTStruct{2,4}), median(RTStruct{2,5}),...
    median(RTStruct{2,6}), median(RTStruct{2,7}), median(RTStruct{2,8}), ... 
    median(RTStruct{2,9}), median(RTStruct{2,10}), median(RTStruct{2,11}), median(RTStruct{2,12}), median(RTStruct{2,13})];   
else
    RTStruct(2,7) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,8) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,9)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,10) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,11) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,12) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,13) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    RTStruct(2,14) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 1)).RT]};
    y_avg_med_RT{2,1} = [mean(RTStruct{2,1}), mean(RTStruct{2,2}), mean(RTStruct{2,3}), mean(RTStruct{2,4}), mean(RTStruct{2,5}),...
        mean(RTStruct{2,6}), mean(RTStruct{2,7}), mean(RTStruct{2,8}), ...
        mean(RTStruct{2,9}), mean(RTStruct{2,10}), mean(RTStruct{2,11}), mean(RTStruct{2,12}), mean(RTStruct{2,13}), mean(RTStruct{2,14})];
    y_avg_med_RT{2,2} = [median(RTStruct{2,1}), median(RTStruct{2,2}), median(RTStruct{2,3}), median(RTStruct{2,4}), median(RTStruct{2,5}),...
        median(RTStruct{2,6}), median(RTStruct{2,7}), median(RTStruct{2,8}), ...
        median(RTStruct{2,9}), median(RTStruct{2,10}), median(RTStruct{2,11}), median(RTStruct{2,12}), median(RTStruct{2,13}), median(RTStruct{2,14})];
end

% For incorrect
RTStruct(3,1)= {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
RTStruct(3,2) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
RTStruct(3,3) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
RTStruct(3,4) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
RTStruct(3,5) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
RTStruct(3,6) = {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
if coh(end) == 0
    RTStruct(3,7) = {[[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT], [Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]]};
    RTStruct(3,8)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,9) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,10) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,11) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,12) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,13) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    y_avg_med_RT{3,1} = [mean(RTStruct{3,1}), mean(RTStruct{3,2}), mean(RTStruct{3,3}), mean(RTStruct{3,4}), mean(RTStruct{3,5}),...
    mean(RTStruct{3,6}), mean(RTStruct{3,7}), mean(RTStruct{3,8}), ... 
    mean(RTStruct{3,9}), mean(RTStruct{3,10}), mean(RTStruct{3,11}), mean(RTStruct{3,12}), mean(RTStruct{3,13})];   
    y_avg_med_RT{3,2} = [median(RTStruct{3,1}), median(RTStruct{3,2}), median(RTStruct{3,3}), median(RTStruct{3,4}), median(RTStruct{3,5}),...
    median(RTStruct{3,6}), median(RTStruct{3,7}), median(RTStruct{3,8}), ... 
    median(RTStruct{3,9}), median(RTStruct{3,10}), median(RTStruct{3,11}), median(RTStruct{3,12}), median(RTStruct{3,13})];   
else
    RTStruct(3,7) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 1 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,8) = {[Trial_4(find([Trial_4(:).difficulty] == 7 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,9)= {[Trial_4(find([Trial_4(:).difficulty] == 6 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,10) = {[Trial_4(find([Trial_4(:).difficulty] == 5 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,11) = {[Trial_4(find([Trial_4(:).difficulty] == 4 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,12) = {[Trial_4(find([Trial_4(:).difficulty] == 3 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,13) = {[Trial_4(find([Trial_4(:).difficulty] == 2 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    RTStruct(3,14) = {[Trial_4(find([Trial_4(:).difficulty] == 1 & [Trial_4(:).target] == 2 & [Trial_4(:).correct] == 0)).RT]};
    y_avg_med_RT{3,1} = [mean(RTStruct{3,1}), mean(RTStruct{3,2}), mean(RTStruct{3,3}), mean(RTStruct{3,4}), mean(RTStruct{3,5}),...
        mean(RTStruct{3,6}), mean(RTStruct{3,7}), mean(RTStruct{3,8}), ...
        mean(RTStruct{3,9}), mean(RTStruct{3,10}), mean(RTStruct{3,11}), mean(RTStruct{3,12}), mean(RTStruct{3,13}), mean(RTStruct{3,14})];
    y_avg_med_RT{3,2} = [median(RTStruct{3,1}), median(RTStruct{3,2}), median(RTStruct{3,3}), median(RTStruct{3,4}), median(RTStruct{3,5}),...
        median(RTStruct{3,6}), median(RTStruct{3,7}), median(RTStruct{3,8}), ...
        median(RTStruct{3,9}), median(RTStruct{3,10}), median(RTStruct{3,11}), median(RTStruct{3,12}), median(RTStruct{3,13}), median(RTStruct{3,14})];
end

%%
%Cleaning out Trial_4 to save 
takeoutfields = {'eCodes', 'timeindexes', 'times', 'RFx', 'RFy', 'first_saccade_time', 'first_saccade_index', 'all_vertical_positions_sized'...
    'all_horizontal_positions_sized', 'alltimes_index_sized', 'velocity_vertical', 'velocity_horizontal', 'velocity_diagonal', 'accelerations',...
    'GPOn', 'FPOff'};
Trial_4 = rmfield(Trial_4, takeoutfields);


end 