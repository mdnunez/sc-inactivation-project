%This script is essentially the same as the decision_task_all_session_mus.m
%except the parts of the code that deal with 2 separate sets of coherences.
%For the saline injections, only 1 set of coherences were performed
clc;
clear all;
close all;

%When you load the GPData (processed from
%decisiontask_masterscript_plotter.m), load the files in chronological
%order 
%---------MODIFY MAIN SETTINGS HERE --------------------------------------%
SP = 0; %1 if all the data are from monkey S, else 0
BT = 0; %1 if all the data are from monkey B, else 0
Both = 1; %1 if all the data are from monkey B and S, else 0
%-------------------------------------------------------------------------%
%% Loading Data
% Global variables
continueAction = 1;
files_inputted = 1;
i = 1;
C = {'r','b','g','d','y'}; % Cell array of colors for the line plotted.
M = {'o','s','d','h','p'}; % Cell array of markers shape.

%Open previously loaded data?
if input('Open Previously loaded data?', 's') == 'y'
    clc
    [fileName, pathname] = getFile('*.*');
    file = [pathname fileName];
    load(file)
else
    %Load Data
    while(continueAction)
        clc
        disp('Open GP Data')
        [fileName, pathname] = getFile('*.*');
        file = [pathname fileName];
        AllGPdata{i,1} = load(file);
        
        if input('Open one more file? y/n', 's') == 'y'
            files_inputted = files_inputted +1;
            continueAction = 1;
            i = i+1;
        else
            continueAction = 0;
        end
    end
    
end

%% Extracting out data from all the injections

coh1_inact = [-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36]; %only 1 set of coherences were performed with saline injections
coh1_acc = [0, 3, 5, 10, 17, 24, 36];
conditions = 3; % pre, post, rec

for i = 1:conditions
    AllPV_inRF_mean_list{i,1} = [];
    AllPV_inRF_std_list{i,1} = [];
    AllPV_outRF_mean_list{i,1} = [];
    AllPV_outRF_std_list{i,1} = [];
    AllPV_inRF_list{i,1} = [];
    AllPV_outRF_list{i,1} = [];
    
    AllPV_inRF_mean_list_RT{i,1} = [];
    AllPV_inRF_std_list_RT{i,1} = [];
    AllPV_outRF_mean_list_RT{i,1} = [];
    AllPV_outRF_std_list_RT{i,1} = [];
    AllPV_inRF_list_RT{i,1} = [];
    AllPV_outRF_list_RT{i,1} = [];
    
    AllRT_inRF_mean_list{i,1} = [];
    AllRT_inRF_std_list{i,1} = [];
    AllRT_outRF_mean_list{i,1} = [];
    AllRT_outRF_std_list{i,1} = [];
    AllRT_inRF_list{i,1} = [];
    AllRT_outRF_list{i,1} = [];
    
    AllAcc_inRF_list{i,1} = [];
    AllAcc_outRF_list{i,1} = [];
    
    AllAlpha{i,1} = [];
    AllBeta{i,1} = [];
    
    %for all the mean values of YInact of each experiment for PMF
    for c_inact = 1:numel(coh1_inact)    %for each coherence in inact PMF (signed)
        AllYInactProp_coh1{i,c_inact} = [];
        Ultim_YInactProp_coh1{i,c_inact} = [];
        AllYInactNum_coh1{i,c_inact} = [];
        Ultim_YInactNum_coh1{i,c_inact} = [];
        AllYInactTottrials_coh1{i,c_inact} = [];
        Ultim_YInactTottrials_coh1{i, c_inact} = [];
        AllRT_coh1{i,c_inact} = [];
        Ultim_RT_coh1{1,c_inact} = [];
    end
    %for all the mean values of YAcc of each experiment for AccPMF
    for c_acc = 1:numel(coh1_acc)     %for each coherence in acc PMF (unsigned)
        AllAccProp_coh1{i,c_acc} = [];
        Ultim_AccProp_coh1{i,c_acc} = [];
        AllYAccNum_coh1{i, c_acc} = [];
        Ultim_YAccNum_coh1{i,c_acc} = [];
        AllYAccTottrials_coh1{i,c_acc} = [];
        Ultim_YAccTottrials_coh1{i, c_acc} = [];
    end
end

for i = 1:numel(AllGPdata)
    for c = 1:numel(coh1_inact)
        RT_Pre{i,c} = [];
        RT_Post{i,c} = [];
        RT_Rec{i,c} = [];
    end
end

filenames_coh1 = [];

coh_dc = unique(abs(coh1_inact));
coh_dc = coh_dc(2:end);
for s = 1:conditions
    givenawayIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    awayIFgivenawayIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    giventoIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    awayIFgiventoIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    
end

for i = 1:numel(AllGPdata)  %for each experiment
    
    %Find the side of the RF
    if AllGPdata{i,1}.GPData(1).RFx < 0    %negative RF is left RF
        RF_side = 1;     % 1 for left RF
    elseif AllGPdata{i,1}.GPData(1).RFx > 0     %positive RF is right RF
        RF_side = 2;     % 2 for right RF
    end
    
    for s = 1:conditions     %for pre, post and rec
        PV_inRF_list = [];
        PV_inRF_list_RT = [];
        PV_outRF_list = [];
        PV_outRF_list_RT = [];
        RT_inRF_list = [];
        RT_outRF_list = [];
        Outcome_inRF_list = [];
        Outcome_outRF_list = [];
        
        for t = 1:numel(AllGPdata{i,1}.GPData(s).Trial_data_array)     %for each trial
            
            %separate peak velocities by inRF and outRF
            %inRF
            if AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side          %if the tg is in the RF
                %Getting PV for inRF
                if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                        PV_inRF_list = [PV_inRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];
                    end
                end
                %Getting RT for inRF
                if AllGPdata{i,1}.GPData(s).Task_version == 'RT'
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT)
                        RT_inRF_list = [RT_inRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                        if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                            PV_inRF_list_RT = [PV_inRF_list_RT, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];
                        end
                        end 
                    end
                end
                %outRF
            elseif AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side
                %Getting PV for outRF
                if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                        PV_outRF_list = [PV_outRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];
                    end
                end
                %Getting RT and PV for outRF
                if AllGPdata{i,1}.GPData(s).Task_version == 'RT'
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT)
                        RT_outRF_list = [RT_outRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                        if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                            PV_outRF_list_RT = [PV_outRF_list_RT, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];
                        end
                        end 
                    end
                end
                
            end
            
            %Getting accuracy for toIF and awayIF trials
            if AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side   %if toIF trial
                %Getting outcome of all non-0 coherence trials; 0 incorrect, 1 correct
                if AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence ~= 0
                    Outcome_inRF_list = [Outcome_inRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct];
                end
            elseif AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side  %if awayIF trial
                %Getting outcome; 0 incorrect, 1 correct
                if AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence ~= 0
                    Outcome_outRF_list = [Outcome_outRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct];
                end
            end
            
            
            
            %getting correct RT's for chronometric functions
            if (AllGPdata{i,1}.GPData(s).Task_version == 'RT')
                if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1
                    %find signed coherence value
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side  %if toIF
                        trial_coh = AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence; %coherence is positive
                    else AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side %if awayIF
                        trial_coh = -AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence; %coherence is negative
                    end
                    
                    if s == 1
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Pre{i,c} = [RT_Pre{i,c},AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                            end
                        end
                    elseif s == 2
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Post{i,c} = [RT_Post{i,c}, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                            end
                        end
                    elseif s == 3
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Rec{i,c} = [RT_Rec{i,c}, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                            end
                        end
                    end
                    
                end
            end
            
            
            
            %Getting awayIF choices when given awayIF stim and toIF stim
            for c = 1:numel(coh_dc)  % for each unsigned coherence,
                if AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence == coh_dc(c)  %check if the current trial has the specified unsigned coherence
                    
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side %if tg is NOT in rf side (if tg is NOT set to toIF, if tg is set to awayIF)
                        givenawayIF{s}(i,c) = givenawayIF{s}(i,c) + 1;  %keep a running count of all give awayIF trials per coherence
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1  % and got it correct awayIF, answered awayIF given awayIF stim
                            awayIFgivenawayIF{s}(i,c) = awayIFgivenawayIF{s}(i,c) + 1;  %keep a running count of all awayIF choices given awayIF stim
                        end
                        
                    elseif AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side  %if tg is set toIF
                        giventoIF{s}(i,c) = giventoIF{s}(i,c) + 1;  %total num of toIF set trials
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 0      %got incorrect toIF, answered awayIF given toIF
                            awayIFgiventoIF{s}(i,c) = awayIFgiventoIF{s}(i,c) + 1; %keep a running count of all awayIF choices given toIF
                        end
                    end
                    
                end
            end
            
        end
        
        %getting a list of all the data across experiments

        %toIF (inRF) 
        %PVs
        AllPV_inRF_list{s,1} = [AllPV_inRF_list{s,1}, PV_inRF_list];        % _list is all the peak velocities, collapsed together across experiments, separated by pre, post, and rec 
        AllPV_inRF_array{s,i} = PV_inRF_list;                               % _array is all the peak velocities, separated into pre, post, rec, and each experiment
        AllPV_inRF_mean_list{s,1} = [AllPV_inRF_mean_list{s,1}, mean(PV_inRF_list)];    
        AllPV_inRF_std_list{s,1} = [AllPV_inRF_std_list{s,1}, std(PV_inRF_list)];
        AllPV_inRF_list_RT{s,1} = [AllPV_inRF_list_RT{s,1}, PV_inRF_list_RT];        
        AllPV_inRF_array_RT{s,i} = PV_inRF_list_RT;                               
        AllPV_inRF_mean_list_RT{s,1} = [AllPV_inRF_mean_list_RT{s,1}, mean(PV_inRF_list_RT)];   
        AllPV_inRF_std_list_RT{s,1} = [AllPV_inRF_std_list_RT{s,1}, std(PV_inRF_list_RT)];
        %RTs
        AllRT_inRF_list{s,1} = [AllRT_inRF_list{s,1}, RT_inRF_list];        
        AllRT_inRF_array{s,i} = RT_inRF_list;                               
        AllRT_inRF_mean_list{s,1} = [AllRT_inRF_mean_list{s,1}, mean(RT_inRF_list)];    
        AllRT_inRF_std_list{s,1} = [AllRT_inRF_std_list{s,1}, std(RT_inRF_list)];
        %Acc
        AllAcc_inRF_list{s,1} = [AllAcc_inRF_list{s,1}, sum(Outcome_inRF_list)/numel(Outcome_inRF_list)];
        AllOutcome_inRF_arrays{s,i} = Outcome_inRF_list;
        
        %awayIF (outRF)
        %PVs
        AllPV_outRF_list{s,1} = [AllPV_outRF_list{s,1}, PV_outRF_list];
        AllPV_outRF_array{s,i} = PV_outRF_list;
        AllPV_outRF_mean_list{s,1} = [AllPV_outRF_mean_list{s,1}, mean(PV_outRF_list)];
        AllPV_outRF_std_list{s,1} = [AllPV_outRF_std_list{s,1}, std(PV_outRF_list)];
        AllPV_outRF_list_RT{s,1} = [AllPV_outRF_list_RT{s,1}, PV_outRF_list_RT];
        AllPV_outRF_array_RT{s,i} = PV_outRF_list_RT;
        AllPV_outRF_mean_list_RT{s,1} = [AllPV_outRF_mean_list_RT{s,1}, mean(PV_outRF_list_RT)];
        AllPV_outRF_std_list_RT{s,1} = [AllPV_outRF_std_list_RT{s,1}, std(PV_outRF_list_RT)];
        %RTs
        AllRT_outRF_list{s,1} = [AllRT_outRF_list{s,1}, RT_outRF_list];
        AllRT_outRF_array{s,i} = RT_outRF_list;
        AllRT_outRF_mean_list{s,1} = [AllRT_outRF_mean_list{s,1}, mean(RT_outRF_list)];
        AllRT_outRF_std_list{s,1} = [AllRT_outRF_std_list{s,1}, std(RT_outRF_list)];
        %Accs
        AllAcc_outRF_list{s,1} = [AllAcc_outRF_list{s,1}, sum(Outcome_outRF_list)/numel(Outcome_outRF_list)];
        AllOutcome_outRF_arrays{s,i} = Outcome_outRF_list;
        
        %Getting together all the YInactnum (number of choices toIF) and
        %YInactTottrials (total number of trials) for PMF fitting and the
        %YAcc (number of trials correct) for accuracy 
        for c_inact = 1:numel(coh1_inact)
            AllYInactNum_coh1{s,c_inact} = [AllYInactNum_coh1{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact)];
            AllYInactTottrials_coh1{s,c_inact} = [AllYInactTottrials_coh1{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact)];
            AllYInactProp_coh1{s,c_inact} = [AllYInactProp_coh1{s,c_inact}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact))];
        end
        %Get the AccMeans for pre/post/rec for individual experiments
        for c_acc = 1:numel(coh1_acc)
            AllYAccNum_coh1{s,c_acc} = [AllYAccNum_coh1{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc)];
            AllYAccTottrials_coh1{s,c_acc} = [AllYAccTottrials_coh1{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc)];
            AllAccProp_coh1{s,c_acc} = [AllAccProp_coh1{s,c_acc}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc))];
        end
        filenames_coh1{i,1} = AllGPdata{i,1}.GPData(s).Filenames;
        
        %Get Alpha and Beta and Lambda Values for each
        AllAlpha{s} = [AllAlpha{s}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams(1)];
        AllBeta{s} = [AllBeta{s}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams(2)];
    end
end

%Getting mean RT's per coherence for CMF (chronometric function)
mean_RT_Pre = zeros(4,13);  % 4 = number of RT task experiments for saline injections, 13 = number of coherences
mean_RT_Post = zeros(4,13);
mean_RT_Rec = zeros(4,13);
x = 1;
%Getting mean RT's
for i = 1:numel(AllGPdata)
    if ~isempty(RT_Pre{i,1})
        for c = 1:numel(coh1_inact)
            mean_RT_Pre(x, c) = mean(RT_Pre{i,c});
            mean_RT_Post(x, c) = mean(RT_Post{i,c});
            mean_RT_Rec(x, c) = mean(RT_Rec{i,c});
        end
        
        x = x + 1;
    end
end
%Getting the ultimate mean RTs for each coherence collapsed across all
%sessions
for c = 1:numel(coh1_inact)
    Ultim_mean_RT_Pre(c) = mean(mean_RT_Pre(:,c));
    Ultim_mean_RT_Post(c) = mean(mean_RT_Post(:,c));
    Ultim_mean_RT_Rec(c) = mean(mean_RT_Rec(:,c));
end
Ultim_mean_RT = [Ultim_mean_RT_Pre; Ultim_mean_RT_Post; Ultim_mean_RT_Rec];


%for computing d' and criterion collapsed across all sessions in SDT
%framework
mu = 0;
sigma = 1;
pd = makedist('Normal',mu,sigma);
%for computing d' and c for each session
for s = 1:conditions
    for c = 1:size(givenawayIF{s},2)   %add up all numbers across sessions but not across coherences for plotting d' and c vs coherence collapsed across sessions
        total_givenawayIF{s}(1,c) = sum(givenawayIF{s}(:,c));
        total_giventoIF{s}(1,c) = sum(giventoIF{s}(:,c));
        total_awayIFgivenawayIF{s}(1,c) = sum(awayIFgivenawayIF{s}(:,c));
        total_awayIFgiventoIF{s}(1,c) = sum(awayIFgiventoIF{s}(:,c));
        
    end
    
    for i = 1:size(givenawayIF{s},1) %add up all numbers across coherences but not across sessions, for plotting individual session's d' and c collapsed across coherences
        givenawayIF_sessions{s}(1,i) = sum(givenawayIF{s}(i,:));
        giventoIF_sessions{s}(1,i) = sum(giventoIF{s}(i,:));
        awayIFgivenawayIF_sessions{s}(1,i) = sum(awayIFgivenawayIF{s}(i,:));
        awayIFgiventoIF_sessions{s}(1,i) = sum(awayIFgiventoIF{s}(i,:));
    end
    
    %collapsing across all coherences, exlcuding out the 0 coherence --
    %for the criterion measurement
    total_givenawayIF_col{s} = sum(total_givenawayIF{s});
    total_giventoIF_col{s} = sum(total_giventoIF{s});
    total_awayIFgivenawayIF_col{s} = sum(total_awayIFgivenawayIF{s});
    total_awayIFgiventoIF_col{s} = sum(total_awayIFgiventoIF{s});
    
    % Gettting the proportion of awayIF/givenawayIF ("Hits") and
    % proportion awayIF/toIF ("FA") for each coherence but collapsed
    % across sessions
    total_proportion_awayIFgivenawayIF{s} = total_awayIFgivenawayIF{s}./total_givenawayIF{s};
    total_proportion_awayIFgiventoIF{s} = total_awayIFgiventoIF{s}./total_giventoIF{s};
    %Getting the proportion for each session but collapsed across
    %coherence
    proportion_awayIFgivenawayIF_sessions{s} = awayIFgivenawayIF_sessions{s}./givenawayIF_sessions{s};
    proportion_awayIFgiventoIF_sessions{s} = awayIFgiventoIF_sessions{s}./giventoIF_sessions{s};
    % Getting proportion for collapsed data across all sessions and
    % coherence
    total_proportion_awayIFgivenawayIF_col{s} = total_awayIFgivenawayIF_col{s}/total_givenawayIF_col{s};
    total_proportion_awayIFgiventoIF_col{s} = total_awayIFgiventoIF_col{s}./total_giventoIF_col{s};
    
    d_prime{s} = (1/sqrt(2)).*(icdf(pd,total_proportion_awayIFgivenawayIF{s}) - icdf(pd,total_proportion_awayIFgiventoIF{s}));
    d_prime_sessions{s} = (1/sqrt(2)).*(icdf(pd,proportion_awayIFgivenawayIF_sessions{s}) - icdf(pd,proportion_awayIFgiventoIF_sessions{s}));
    d_prime_col{s} = (1/sqrt(2)).*(icdf(pd,total_proportion_awayIFgivenawayIF_col{s}) - icdf(pd,total_proportion_awayIFgiventoIF_col{s}));
    criterion{s} = -0.5.*(icdf(pd,total_proportion_awayIFgivenawayIF{s}) + icdf(pd,total_proportion_awayIFgiventoIF{s}));
    criterion_sessions{s} = -0.5.*(icdf(pd,proportion_awayIFgivenawayIF_sessions{s}) + icdf(pd,proportion_awayIFgiventoIF_sessions{s}));
    criterion_col{s} = -0.5.*(icdf(pd,total_proportion_awayIFgivenawayIF_col{s}) + icdf(pd,total_proportion_awayIFgiventoIF_col{s}));
    
end
dprimec = struct('dprime_coh', [d_prime], 'c_coh', [criterion], 'dprime_collapsed', [d_prime_col],'c_collapsed', [criterion_col], 'dprime_sessions', [d_prime_sessions],'criterion_sesssions', [criterion_sessions]);

%%
%Ultimate Means across all experiments
Ultim_MeanPV_inRF_pre = mean(AllPV_inRF_list{1,1});
Ultim_MeanPV_inRF_post = mean(AllPV_inRF_list{2,1});
Ultim_MeanPV_inRF_rec = mean(AllPV_inRF_list{3,1});
Ultim_MeanPV_outRF_pre = mean(AllPV_outRF_list{1,1});
Ultim_MeanPV_outRF_post = mean(AllPV_outRF_list{2,1});
Ultim_MeanPV_outRF_rec = mean(AllPV_outRF_list{3,1});

Ultim_MeanPV_inRF_pre_RT = mean(AllPV_inRF_list_RT{1,1});
Ultim_MeanPV_inRF_post_RT = mean(AllPV_inRF_list_RT{2,1});
Ultim_MeanPV_inRF_rec_RT = mean(AllPV_inRF_list_RT{3,1});
Ultim_MeanPV_outRF_pre_RT = mean(AllPV_outRF_list_RT{1,1});
Ultim_MeanPV_outRF_post_RT = mean(AllPV_outRF_list_RT{2,1});
Ultim_MeanPV_outRF_rec_RT = mean(AllPV_outRF_list_RT{3,1});

Ultim_MeanRT_inRF_pre = mean(AllRT_inRF_list{1,1});
Ultim_MeanRT_inRF_post = mean(AllRT_inRF_list{2,1});
Ultim_MeanRT_inRF_rec = mean(AllRT_inRF_list{3,1});
Ultim_MeanRT_outRF_pre = mean(AllRT_outRF_list{1,1});
Ultim_MeanRT_outRF_post = mean(AllRT_outRF_list{2,1});
Ultim_MeanRT_outRF_rec = mean(AllRT_outRF_list{3,1});

Ultim_MeanAcc_inRF_pre = mean(AllAcc_inRF_list{1,1});
Ultim_MeanAcc_inRF_post = mean(AllAcc_inRF_list{2,1});
Ultim_MeanAcc_outRF_pre = mean(AllAcc_outRF_list{1,1});
Ultim_MeanAcc_outRF_post = mean(AllAcc_outRF_list{2,1});

for s = 1:conditions
    for c_inact = 1:numel(coh1_inact)
        Ultim_YInactNum{s,c_inact} = sum(AllYInactNum_coh1{s,c_inact});
        Ultim_YInactTottrials{s,c_inact} = sum(AllYInactTottrials_coh1{s,c_inact});
        Ultim_YInactProp{s,c_inact} = mean(AllYInactProp_coh1{s,c_inact});
        Ultim_RT{s,c_inact} = mean(AllRT_coh1{s,c_inact});
    end
    for c_acc = 1:numel(coh1_acc)
        Ultim_YAccNum{s,c_acc} = sum(AllYAccNum_coh1{s,c_acc});
        Ultim_YAccTottrials{s,c_acc} = sum(AllYAccTottrials_coh1{s,c_acc});
        Ultim_AccProp{s,c_acc} = mean(AllAccProp_coh1{s,c_acc});
    end
    
end
%% Get the logistic fits to the Ultimate mean PMF
%Use the Logistic function
PF = @PAL_Logistic;
%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = -36:1:36;
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma =  0;  %LambdaL; scalar here (since fixed) but may be vector
searchGrid.lambda = 0;  %LambdaR; ditto

x_fit_inact_36 = linspace(coh1_inact(1),coh1_inact(end),100);


%INACT
%Fit for every pre post rec in every exp - for all in RIGHT
for s = 1:3
    disp('Fitting Ultimate Inact pre/post/rec function.....');
    [paramsValuesUltimInact(s,:) LL exitflag] = PAL_PFML_Fit(coh1_inact,[Ultim_YInactNum{s,:}], ...
        [Ultim_YInactTottrials{s,:}],searchGrid,paramsFree,PF)
    y_HatML_inact(s,:) = (PAL_Logistic(paramsValuesUltimInact(s,:), x_fit_inact_36)).*100;
    
end
stop = 1;
%For all PMF fits with coherences (36 etc..)
for i = 1:numel(AllGPdata)
    for s = 1:3
        AllGPdata{i,1}.GPData(s).YAndParams_array.Y_HatML_Inact = (PAL_Logistic(AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams, x_fit_inact_36)).*100;
    end
end
%% Fitting linefit for RT line

x_fit_inactcoh1_1 = linspace(coh1_inact(1), coh1_inact(7),100);
x_fit_inactcoh1_2 = linspace(coh1_inact(7), coh1_inact(end), 100);

%get the RT's, organized by coherence (c) and condition (s)
for s = 1:conditions
    for c = 1:numel(coh1_inact)
        if s == 1
            AllRT_coh1{s,c} = mean_RT_Pre(:,c)';
        elseif s == 2
            AllRT_coh1{s,c} = mean_RT_Post(:,c)';
        elseif s == 3
            AllRT_coh1{s,c} = mean_RT_Rec(:,c)';
        end
    end
end
%reorganizing RT structure for linear fitting
AllRT_1{1,1} = [];
AllRT_1{2,1} = [];
AllRT_1{3,1} = [];
for c = 1:size(AllRT_coh1,2)
    for s = 1:size(AllRT_coh1, 1)
        AllRT_1{s,1} = [AllRT_1{s,1};AllRT_coh1{s,c}];
    end
end

linfit_yint_pre = [];
linfit_slope_pre= [];
linfit_yint_post = [];
linfit_slope_post = [];
linfit_yint_rec = [];
linfit_slope_rec = [];

%fit each individual session
for i = 1:4 %for each saline injection RT experiment
    for s = 1:3 %for each condition (pre, post, rec)
        %get RT's
        RT{s,1} = [];
        RT{s,1} = AllRT_1{s,1}(:,i)';
        %get linear fits
        x_1 = [ones(length(coh1_inact(1:7)),1),coh1_inact(1:7)'];
        x_2 = [ones(length(coh1_inact(7:end)),1),coh1_inact(7:end)'];
        b_1(s,:) = (x_1\RT{s,1}(1:7)')';
        b_2(s,:) = (x_2\RT{s,1}(7:end)')';
        y_linfit_1(s,:) = b_1(s,1) + b_1(s,2)*x_fit_inactcoh1_1;
        y_linfit_2(s,:) = b_2(s,1) + b_2(s,2)*x_fit_inactcoh1_2;
    end
    y_linefit_all{i,1} = y_linfit_1;
    y_linefit_all{i,2} = y_linfit_2;
    linfit_yint_pre = [linfit_yint_pre; b_1(1,1), b_2(1,1)];            %1st col is outRF, 2nd col is inRF
    linfit_slope_pre  = [linfit_slope_pre; b_1(1,2), b_2(1,2)];
    linfit_yint_post = [linfit_yint_post; b_1(2,1), b_2(2,1)];
    linfit_slope_post = [linfit_slope_post; b_1(2,2), b_2(2,2)];
    linfit_yint_rec = [ linfit_yint_rec; b_1(3,1), b_2(3,1)];
    linfit_slope_rec = [linfit_slope_rec; b_1(3,2), b_2(3,2)];
end
RT_linfit_params = struct('Slope', [], 'Int', []);
RT_linfit_params(1).Slope = linfit_slope_pre;
RT_linfit_params(2).Slope = linfit_slope_post;
RT_linfit_params(3).Slope = linfit_slope_rec;
RT_linfit_params(1).Int = linfit_yint_pre;
RT_linfit_params(2).Int = linfit_yint_post;
RT_linfit_params(3).Int = linfit_yint_rec;

%fit the mean RT line
for s = 1:3
    %get linear fits
    x_1 = [ones(length(coh1_inact(1:7)),1),coh1_inact(1:7)'];
    x_2 = [ones(length(coh1_inact(7:end)),1),coh1_inact(7:end)'];
    b_1 = x_1\[Ultim_mean_RT(s,1:7)]';      %first element is y-int, second element is slope
    b_2 = x_2\[Ultim_mean_RT(s,7:end)]';
    y_linfit_1_ultimate(s,:) = b_1(1) + b_1(2)*x_fit_inactcoh1_1;
    y_linfit_2_ultimate(s,:) = b_2(1) + b_2(2)*x_fit_inactcoh1_2;
end



%% Plotting parameters 

Linewidth = 8;
Linewidth_1 = 4;
dotsz = 100;
dotsz_1 = 70;
pre_color_mean = [0 0 0];
pre_color_indiv = [86/255, 86/255, 86/255];
post_color_mean = [0.4 0.4, 0.7];
post_color_indiv = [179/255, 179/255, 230/255];
rec_color_mean = [232/255, 104/255, 0];
rec_color_indiv = [255/255, 177/255, 96/255];

Linewidth_1 = 3;

toIF_color_mean = [0, 135/255, 255/255];
toIF_color_indiv = [153/255, 255/255, 255/255];
awayIF_color_mean = [205/255, 57/255, 136/255];
awayIF_color_indiv = [255/255, 180/255, 204/255];
if SP == 1
    marker = '^';
elseif BT == 1
    marker = 'o';
else
    marker = 'o'
end

%% Plot Peak Velocity inRF (toIF) vs OutRF (awayIF) for all injections

x1 = zeros(1,numel(AllPV_inRF_mean_list{1,1}))+1;
x2 = x1 +1;
%Plotter for PV inRF vs OutRF
figure(1);
hold on
scatter(x1, AllPV_inRF_mean_list{1,1}, 50, toIF_color_indiv, marker,'filled');
scatter(x1, AllPV_outRF_mean_list{1,1},50, awayIF_color_indiv, marker,'filled');
scatter(x2, AllPV_inRF_mean_list{2,1},50, toIF_color_indiv,marker, 'filled');
scatter(x2, AllPV_outRF_mean_list{2,1},50, awayIF_color_indiv,marker,'filled');
scatter(1,Ultim_MeanPV_inRF_pre, 100, toIF_color_mean,marker, 'filled')
scatter(1,Ultim_MeanPV_outRF_pre, 100, awayIF_color_mean,marker, 'filled')
scatter(2, Ultim_MeanPV_inRF_post, 100, toIF_color_mean,marker, 'filled');
scatter(2, Ultim_MeanPV_outRF_post, 100, awayIF_color_mean,marker, 'filled');
plot([1,2], [Ultim_MeanPV_inRF_pre, Ultim_MeanPV_inRF_post], marker, 'color', toIF_color_mean, 'linestyle', '-','linewidth', Linewidth_1);
plot([1,2], [Ultim_MeanPV_outRF_pre, Ultim_MeanPV_outRF_post],marker, 'color', awayIF_color_mean,'linestyle', '-', 'linewidth', Linewidth_1);
set(gca, 'xtick', [1,2], 'xticklabel', {'Pre', 'Post'}, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('PV', 'Fontsize', 30, 'FontWeight', 'bold');
xlim([0.7, 2.3]);
ylim([200, 1000])
legend('toIF', 'awayIF', 'location', 'south');
hold off;


%% Plot Accuracy inRF (toIF) vs OutRF (awayIF) for all injections

x1 = zeros(1,numel(AllAcc_inRF_list{1,1}))+1;
x2 = x1 +1;
%Plotter for Acc inRF vs OutRF
figure(36);
% plot([1,1], [Ultim_MeanAcc_inRF_pre,Ultim_MeanAcc_outRF_pre], 'linewidth', 2, 'color', pre_color_mean);
hold on
% plot([2,2], [Ultim_MeanAcc_inRF_post,Ultim_MeanAcc_outRF_post], 'linewidth', 2, 'color', pre_color_mean);
scatter(x1, AllAcc_inRF_list{1,1}, 50, toIF_color_indiv, marker,'filled');
scatter(x1, AllAcc_outRF_list{1,1},50, awayIF_color_indiv, marker,'filled');
scatter(x2, AllAcc_inRF_list{2,1},50, toIF_color_indiv,marker, 'filled');
scatter(x2, AllAcc_outRF_list{2,1},50, awayIF_color_indiv,marker,'filled');
scatter(1,Ultim_MeanAcc_inRF_pre, 100, toIF_color_mean,marker, 'filled')
scatter(1,Ultim_MeanAcc_outRF_pre, 100, awayIF_color_mean,marker, 'filled')
scatter(2, Ultim_MeanAcc_inRF_post, 100, toIF_color_mean,marker, 'filled');
scatter(2, Ultim_MeanAcc_outRF_post, 100, awayIF_color_mean,marker, 'filled');
plot([1,2], [Ultim_MeanAcc_inRF_pre, Ultim_MeanAcc_inRF_post], marker, 'color', toIF_color_mean, 'linestyle', '-','linewidth', Linewidth_1);
plot([1,2], [Ultim_MeanAcc_outRF_pre, Ultim_MeanAcc_outRF_post],marker, 'color', awayIF_color_mean,'linestyle', '-', 'linewidth', Linewidth_1);
set(gca, 'xtick', [1,2], 'xticklabel', {'Pre', 'Post'}, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Acc', 'Fontsize', 30, 'FontWeight', 'bold');
xlim([0.7, 2.3]);
ylim([0.3, 1])
legend('toIF', 'awayIF', 'location', 'south');
hold off;



%% PMF Plot -- PRE vs POST
stop = 1;
figure(4);
set(gcf,'color','w');
hold on
plot(x_fit_inact_36, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact_36, y_HatML_inact(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, [Ultim_YInactProp{1,:}]*100, dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, [Ultim_YInactProp{2,:}]*100, dotsz, post_color_mean,  'filled');
for i = 1:numel(AllGPdata)
    plot(AllGPdata{i,1}.GPData(1).YAndParams_array.xfit, AllGPdata{i,1}.GPData(1).YAndParams_array.Y_HatML_Inact , 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(AllGPdata{i,1}.GPData(2).YAndParams_array.xfit, AllGPdata{i,1}.GPData(2).YAndParams_array.Y_HatML_Inact , 'color', post_color_indiv, 'Linewidth', Linewidth_1);
    scatter(AllGPdata{i,1}.GPData(1).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(1).YAndParams_array.DataYPropInact, dotsz_1, pre_color_indiv,  'filled');
    scatter(AllGPdata{i,1}.GPData(2).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(2).YAndParams_array.DataYPropInact, dotsz_1, post_color_indiv,  'filled');
    pause = 1;
end
plot(x_fit_inact_36, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact_36, y_HatML_inact(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, [Ultim_YInactProp{1,:}]*100, dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, [Ultim_YInactProp{2,:}]*100, dotsz, post_color_mean,  'filled');
s = plot([0 0],[0 100], ':k');
m = plot([-40, 40], [50 50], ':k');
ticks = [0:10:100];
ylabel('Choices to IF (%)', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold');
ylim([0 100]);
title('Pre VS Post','FontWeight', 'bold');
ylim([0 100]);
set(gca, 'ytick', [0, 20,40, 60, 80, 100], 'yticklabels', [0, 20, 40, 60, 80, 100]);
stop = 1;
%% PMF Plot -- PRE vs REC
figure(6);
set(gcf,'color','w');
hold on
plot(x_fit_inact_36, y_HatML_inact(1,:), 'color',pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact_36, y_HatML_inact(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, [Ultim_YInactProp{1,:}]*100, dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, [Ultim_YInactProp{3,:}]*100, dotsz, rec_color_mean,  'filled');
for i = 1:numel(AllGPdata)
    plot(AllGPdata{i,1}.GPData(1).YAndParams_array.xfit, AllGPdata{i,1}.GPData(1).YAndParams_array.Y_HatML_Inact , 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(AllGPdata{i,1}.GPData(3).YAndParams_array.xfit, AllGPdata{i,1}.GPData(3).YAndParams_array.Y_HatML_Inact , 'color', rec_color_indiv, 'Linewidth', Linewidth_1);
    scatter(AllGPdata{i,1}.GPData(1).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(1).YAndParams_array.DataYPropInact, dotsz_1, pre_color_indiv,  'filled');
    scatter(AllGPdata{i,1}.GPData(3).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(3).YAndParams_array.DataYPropInact, dotsz_1, rec_color_indiv,  'filled');
    pause = 1;
end
plot(x_fit_inact_36, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact_36, y_HatML_inact(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, [Ultim_YInactProp{1,:}]*100, dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, [Ultim_YInactProp{3,:}]*100, dotsz, rec_color_mean,  'filled');
s = plot([0 0],[0 100], ':k');
m = plot([-40, 40], [50 50], ':k');
ticks = [0:10:100];
ylabel('Choices to IF (%)', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold', 'Fontsize', 30);
ylim([0 100]);
title('Pre VS Rec (Saline)','FontWeight', 'bold', 'Fontsize', 30);
fsize=30;
ylim([0 100]);
set(gca, 'ytick', [0, 20,40, 60, 80, 100], 'yticklabels', [0, 20, 40, 60, 80, 100]);
stop = 1;
%% RT chronometric function  Plots -- PRE vs POST
figure(7);
set(gcf,'color','w');
hold on
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, Ultim_mean_RT(1,:), dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, Ultim_mean_RT(2,:), dotsz, post_color_mean,  'filled');

for i = 1:4
    plot(x_fit_inactcoh1_1, y_linefit_all{i,1}(1,:), 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_2, y_linefit_all{i,2}(1,:), 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_1, y_linefit_all{i,1}(2,:), 'color', post_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_2, y_linefit_all{i,2}(2,:), 'color', post_color_indiv, 'Linewidth', Linewidth_1);
end
for i = 1:numel(coh1_inact)
    scatter(zeros(1, numel(AllRT_coh1{1,i}))+coh1_inact(i),AllRT_coh1{1,i}, dotsz_1, pre_color_indiv,  'filled');
    scatter(zeros(1, numel(AllRT_coh1{1,i}))+coh1_inact(i),AllRT_coh1{2,i}, dotsz_1, post_color_indiv,  'filled');
end
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, Ultim_mean_RT(1,:), dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, Ultim_mean_RT(2,:), dotsz, post_color_mean,  'filled');

% fiducials
s = plot([0 0],[0 100], ':k');
set(gca, 'xtick', coh1_inact, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ticks = [0:10:100];
ylabel('RT', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold', 'Fontsize', 30);
ylim([0.4 1.3]);
title('Pre VS Post','FontWeight', 'bold', 'Fontsize', 30);

%% RT chronometric function  Plots -- PRE vs REC
figure(8);
set(gcf,'color','w');
hold on
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, Ultim_mean_RT(1,:), dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, Ultim_mean_RT(3,:), dotsz, rec_color_mean,  'filled');
for i = 1:numel(coh1_inact)
    scatter(zeros(1, numel(AllRT_coh1{1,i}))+coh1_inact(i),AllRT_coh1{1,i}, dotsz_1, pre_color_indiv,  'filled');
    scatter(zeros(1, numel(AllRT_coh1{1,i}))+coh1_inact(i),AllRT_coh1{3,i}, dotsz_1, rec_color_indiv,  'filled');
end
for i = 1:4
    plot(x_fit_inactcoh1_1, y_linefit_all{i,1}(1,:), 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_2, y_linefit_all{i,2}(1,:), 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_1, y_linefit_all{i,1}(3,:), 'color', rec_color_indiv, 'Linewidth', Linewidth_1);
    plot(x_fit_inactcoh1_2, y_linefit_all{i,2}(3,:), 'color', rec_color_indiv, 'Linewidth', Linewidth_1);
end
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_1, y_linfit_1_ultimate(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inactcoh1_2, y_linfit_2_ultimate(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(coh1_inact, Ultim_mean_RT(1,:), dotsz, pre_color_mean, 'filled');
scatter(coh1_inact, Ultim_mean_RT(3,:), dotsz, rec_color_mean,  'filled');
s = plot([0 0],[0 100], ':k');
set(gca, 'xtick', coh1_inact, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ticks = [0:10:100];
ylabel('RT', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold', 'Fontsize', 30);
ylim([0.4 1.3]);
title('Pre VS Rec','FontWeight', 'bold', 'Fontsize', 30);
fsize=20;
title('Pre VS Rec','FontWeight', 'bold', 'Fontsize', 30);
%% Alpha pre vs post vs rec plot 
figure(9);
scatter(AllAlpha{1}, AllAlpha{2}, dotsz, post_color_indiv, 'filled')
hold on
scatter(AllAlpha{1}, AllAlpha{3}, dotsz, rec_color_indiv, 'filled')
scatter(median(AllAlpha{1}), median(AllAlpha{2}), dotsz, post_color_mean, 'filled')
scatter(median(AllAlpha{1}), median(AllAlpha{3}), dotsz, rec_color_mean,'filled')
ylim([-10 30])
xlim([-10 30])
plot([-10:30],[-10:30], ':k');
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold'); 
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Alpha (Saline)','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post', 'Rec'}, 'Fontsize', 20);
%% Beta Plot pre vs post vs rec plot
figure(10);
scatter(AllBeta{1}, AllBeta{2}, dotsz, post_color_indiv, 'filled')
hold on
scatter(AllBeta{1}, AllBeta{3}, dotsz, rec_color_indiv, 'filled')
scatter(median(AllBeta{1}), median(AllBeta{2}), dotsz, post_color_mean, 'filled')
scatter(median(AllBeta{1}), median(AllBeta{3}), dotsz, rec_color_mean,'filled')
ylim([0.05 0.3])
xlim([0.05 0.3])
plot([0.05:0.01:0.3],[0.05:0.01:0.3], ':k');
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Beta (Saline)','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post', 'Rec'}, 'Fontsize', 20);


%% RT slope plot - toIF and awayIF - pre vs post 
figure(15);
scatter(linfit_slope_pre(:,1), linfit_slope_post(:,1), dotsz, awayIF_color_indiv, 's', 'filled');
hold on
scatter(linfit_slope_pre(:,2), linfit_slope_post(:,2),dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_slope_pre(:,1)), median(linfit_slope_post(:,1)), dotsz,awayIF_color_mean,'s', 'filled');
scatter(median(linfit_slope_pre(:,2)), median(linfit_slope_post(:,2)),dotsz, toIF_color_mean, 'filled');

plot([-0.008 0.008 ], [-0.008 0.008 ], ':k')
plot([0, 0 ], [-0.008, 0.008] ,':k')
plot([-0.008, 0.008], [0, 0 ], ':k')
ylim([-0.008 , 0.008 ])
xlim([-0.008 , 0.008 ])

ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Slope of RT linear fit - Pre vs Post (Saline)','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);

%% RT slope plot - toIF and awayIF - pre vs rec 
figure(16);
hold on
scatter(linfit_slope_pre(:,1), linfit_slope_rec(:,1),dotsz,awayIF_color_indiv, 's', 'filled');
scatter(linfit_slope_pre(:,2), linfit_slope_rec(:,2),dotsz,toIF_color_indiv, 'filled');

scatter(median(linfit_slope_pre(:,1)), median(linfit_slope_rec(:,1)),dotsz,awayIF_color_mean, 's', 'filled');
scatter(median(linfit_slope_pre(:,2)), median(linfit_slope_rec(:,2)),dotsz,toIF_color_mean, 'filled');

plot([-0.008 0.008 ], [-0.008 0.008 ], ':k')
plot([0, 0 ], [-0.008, 0.008] ,':k')
plot([-0.008, 0.008], [0, 0 ], ':k')
ylim([-0.008 , 0.008 ])
xlim([-0.008 , 0.008 ])

ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Slope of RT linear fit - Pre vs Rec (Saline)','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);

%%  RT intercept plot - toIF and awayIF - pre vs post 
figure(17);
scatter(linfit_yint_pre(:,1), linfit_yint_post(:,1), dotsz, awayIF_color_indiv, 's', 'filled');
hold on
scatter(linfit_yint_pre(:,2), linfit_yint_post(:,2), dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_yint_pre(:,1)), median(linfit_yint_post(:,1)), dotsz, awayIF_color_mean, 's', 'filled');
scatter(median(linfit_yint_pre(:,2)), median(linfit_yint_post(:,2)), dotsz, toIF_color_mean, 'filled');

plot([0.5 , 1.2 ], [0.5 , 1.2 ], ':k')
ylim([0.5 , 1.2 ])
xlim([0.5 , 1.2 ])
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Intercept of RT linear fit - Pre vs Post(Saline)','FontWeight', 'bold', 'Fontsize', 15);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);

%%  RT intercept plot - toIF and awayIF - pre vs rec 
figure(18);
hold on
scatter(linfit_yint_pre(:,1), linfit_yint_rec(:,1), dotsz,awayIF_color_indiv, 's', 'filled');
scatter(linfit_yint_pre(:,2), linfit_yint_rec(:,2),dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_yint_pre(:,1)), median(linfit_yint_rec(:,1)), dotsz, awayIF_color_mean, 's', 'filled');
scatter(median(linfit_yint_pre(:,2)), median(linfit_yint_rec(:,2)), dotsz,toIF_color_mean, 'filled');

plot([0.5 , 1.2 ], [0.5 , 1.2 ], ':k')
ylim([0.5 , 1.2 ])
xlim([0.5 , 1.2 ])
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Intercept of RT linear fit - Pre vs Rec (Saline)','FontWeight', 'bold', 'Fontsize', 15);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);


%% Plot D' and C in SDT -- separated by coherences and collapsed over coherence 

C = {pre_color_mean,post_color_mean, rec_color_mean};
figure(33);
hold on
for s = 1:conditions
    plot(coh_dc, d_prime{s}, 'Linewidth', Linewidth, 'color', C{s})
    scatter(coh_dc, d_prime{s}, dotsz, C{s},'filled')
    ylim([0 3.5])
    xlim([ 0, 39]);
end
xticks(coh_dc);
ylabel('D-Prime', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Stimulus Difficulty','FontWeight', 'bold', 'Fontsize', 30);
title('D-Prime Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
stop = 1;


figure(34);
hold on
for s = 1:conditions
    plot(coh_dc, criterion{s}, 'Linewidth', Linewidth, 'color', C{s})
    scatter(coh_dc, criterion{s}, dotsz, C{s},'filled')
    ylim([-0.7 0.25])
    xlim([ 0, 39]);
end
ylabel('Criterion', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Stimulus Difficulty','FontWeight', 'bold', 'Fontsize', 30);
title('Criterion Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
stop = 1;



figure(35)
hold on
for s = 1:conditions
    bar(s, mean(d_prime_sessions{s}),'Facecolor',C{s});
end
for i = 1:size(d_prime_sessions{1},2)
    plot([1,2,3], [d_prime_sessions{1}(i), d_prime_sessions{2}(i), d_prime_sessions{3}(i)], '-o');
end
ylim([0 1.8])
ylabel('D-Prime', 'Fontsize', 30, 'FontWeight', 'bold');
xticks([1,2,3])
xticklabels({'Pre', 'Post', 'Rec'});
title('D-Prime Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
legend('Pre', 'Post', 'Rec', 'Fontsize', 30);
stop = 1;

figure(37)
hold on
for s = 1:conditions
    bar(s, mean(criterion_sessions{s}),'Facecolor',C{s});
end
for i = 1:size(criterion_sessions{1},2)
    plot([1,2,3], [criterion_sessions{1}(i), criterion_sessions{2}(i), criterion_sessions{3}(i)], '-o');
end
ylim([-1 0.5])
ylabel('Criterion', 'Fontsize', 30, 'FontWeight', 'bold');
% xlabel('Stimulus Difficulty','FontWeight', 'bold', 'Fontsize', 30);
xticks([1,2,3])
xticklabels({'Pre', 'Post', 'Rec'});
title('Criterion Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
legend('Pre', 'Post', 'Rec', 'Fontsize', 30);
stop = 1;