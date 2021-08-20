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

coh1_inact = [-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36];   %1st set of coherences performed from all sessions from all injections (muscimol,saline) except 3 muscimol injections from monkey S
coh2_inact = [-50, -36, -24, -17, -10, -5, 0, 5, 10, 17, 24, 36, 50]; %2nd set of coherences performed on all session from 3 muscimol injections from monkey S
cohall_inact = [-50, -36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36, 50];  %coherences combined into 1 set for plotting
coh1_acc = [0, 3, 5, 10, 17, 24, 36];
coh2_acc = [0, 5, 10, 17, 24, 36, 50];
cohall_acc = [0, 3, 5, 10, 17, 24, 26, 50];

conditions = 3; % pre, post, rec

%make empty data arrays to fill in (PV = peak velocity, RT = reaction time,
%Acc = accuracy, PMF = psychometric function)
for i = 1:conditions
    AllPV_inRF_mean_list{i,1} = [];
    AllPV_inRF_mean1_list{i,1} = [];
    AllPV_inRF_median_list{i,1} = [];
    AllPV_inRF_median1_list{i,1} = [];
    AllPV_inRF_std_list{i,1} = [];
    AllPV_outRF_mean_list{i,1} = [];
    AllPV_outRF_mean1_list{i,1} = [];
    AllPV_outRF_median_list{i,1} = [];
    AllPV_outRF_median1_list{i,1} = [];
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
    
    %for all the mean values of YInact (proportion of choices toIF for inactivation experiments) of each experiment for PMF
    for c_inact = 1:numel(coh1_inact)    %for each coherence in inact PMF (signed)
        AllYInactProp_coh1{i,c_inact} = []; %list of all the proportion of choices toIF for coh1 for 1st set of coherences for all injections
        AllYInactProp_coh2{i,c_inact} = []; %list of all the proportion of choices toIF coh2 for 2nd set of coherences for all injections 
        Ultim_YInactProp_coh1{i,c_inact} = [];
        Ultim_YInactProp_coh2{i,c_inact} = [];
        AllYInactNum_coh1{i,c_inact} = [];
        AllYInactNum_coh2{i,c_inact} = [];
        Ultim_YInactNum_coh1{i,c_inact} = [];
        Ultim_YInactNum_coh2{i,c_inact} = [];
        AllYInactTottrials_coh1{i,c_inact} = [];
        AllYInactTottrials_coh2{i,c_inact} = [];
        Ultim_YInactTottrials_coh1{i, c_inact} = [];
        Ultim_YInactTottrials_coh2{i, c_inact} = [];
        AllRT_coh1{i,c_inact} = []; %there's only coh1 because only the 1st set of coherences were performed for RT tasks
        Ultim_RT_coh1{1,c_inact} = [];
    end
    %for all the mean values of YAcc (proportion correct) of each
    %experiment for accuracy PMF
    for c_acc = 1:numel(coh1_acc)     %for each coherence in acc PMF (unsigned coherences)
        AllAccProp_coh1{i,c_acc} = [];
        AllAccProp_coh2{i,c_acc} = [];
        Ultim_AccProp_coh1{i,c_acc} = [];
        Ultim_AccProp_coh2{i,c_acc} = [];
        AllYAccNum_coh1{i, c_acc} = [];
        AllYAccNum_coh2{i, c_acc} = [];
        Ultim_YAccNum_coh1{i,c_acc} = [];
        Ultim_YAccNum_coh2{i,c_acc} = [];
        AllYAccTottrials_coh1{i,c_acc} = [];
        AllYAccTottrials_coh2{i,c_acc} = [];
        Ultim_YAccTottrials_coh1{i, c_acc} = [];
        Ultim_YAccTottrials_coh2{i, c_acc} = [];
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
filenames_coh2 = [];
%for d' and criterion calculation
coh_dc = unique(abs(cohall_inact));
coh_dc = coh_dc(2:end);
for s = 1:conditions
    givenawayIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    awayIFgivenawayIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    giventoIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
    awayIFgiventoIF{s} = zeros(numel(AllGPdata),numel(unique(abs(coh_dc))));
end

%% Extracting out data from all injections
for i = 1:numel(AllGPdata)  %for each experiment
    
    %Find the side of the RF
    if AllGPdata{i,1}.GPData(1).RFx < 0    %negative RF is left RF
        RF_side = 1;     % 1 for left RF
    elseif AllGPdata{i,1}.GPData(1).RFx > 0     %positive RF is right RF
        RF_side = 2;     % 2 for right RF
    end
    
    %Find which monkey the data belongs to
    if strfind(AllGPdata{i,1}.GPData(1).Filenames, 'SP')
        SP_file = 1;
        BT_file = 0;
    elseif strfind(AllGPdata{i,1}.GPData(1).Filenames, 'BT')
        SP_file = 0;
        BT_file = 1;
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
            
            %separate peak velocities and RT by toIF and awayIF saccades
            %toIF
            if (AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side &&  AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1) || (AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side &&  AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 0)      %if sac is to IF
                %Getting PV for toIF
                if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity) %if peak velocity was able to be recorded in this trial
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100  %make sure peak velocity is in correct range 
                        PV_inRF_list = [PV_inRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity]; %list of PV values made toIF (inRF)
                    end
                end
                %Getting RT for toIF
                if AllGPdata{i,1}.GPData(s).Task_version == 'RT'
                    
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT) %if RT was able to be recorded in this trial
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT < 50 %no RTs less than 50ms
                            RT_inRF_list = [RT_inRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT]; %list of RT values made toIF (inRF)
                        end
                    end
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                            PV_inRF_list_RT = [PV_inRF_list_RT, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];  %list of peak velocities only for RT trials made toIF (inRF)
                        end
                    end
                    
                end
                
            %awayIF
            elseif (AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side &&  AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1) || (AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side &&  AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 0)
                %Getting PV for outRF
                if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                        PV_outRF_list = [PV_outRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity]; %list of PV values for awayIF (outRF)
                    end
                end
                %Getting RT and PV for outRF
                if AllGPdata{i,1}.GPData(s).Task_version == 'RT'
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT)
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT < 50
                            RT_outRF_list = [RT_outRF_list, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT];
                        end
                    end
                    if ~isempty(AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity)
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity < 1600 && AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity > 100
                            PV_outRF_list_RT = [PV_outRF_list_RT, AllGPdata{i,1}.GPData(s).Trial_data_array(t).PeakVelocity];
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
                %Getting outcome of all non-0 coherence trials; 0 incorrect, 1 correct
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
                    else AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side; %if awayIF
                        trial_coh = -AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence; %coherence is negative
                    end 
                    
                    if s == 1 %if pre session
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Pre{i,c} = [RT_Pre{i,c},AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT]; %list of all correct RT's for pre
                            end 
                        end 
                    elseif s == 2 %if post session
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Post{i,c} = [RT_Post{i,c}, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT]; %list of all correct RT's for post
                            end 
                        end 
                    elseif s == 3 %if rec session
                        for c = 1:numel(coh1_inact)
                            if trial_coh == coh1_inact(c)
                                RT_Rec{i,c} = [RT_Rec{i,c}, AllGPdata{i,1}.GPData(s).Trial_data_array(t).RT]; %list of all correct RT's for rec
                            end 
                        end 
                    end 
                    
                end 
            end 

            
            %Getting awayIF choices when given awayIF stimulus and toIF
            %stimulus for d' and criterion calculation in SDT 
            for c = 1:numel(coh_dc)  % for each unique unsigned coherence,
                if AllGPdata{i,1}.GPData(s).Trial_data_array(t).coherence == coh_dc(c)  %check if the current trial has the specified unsigned coherence
                    
                    if AllGPdata{i,1}.GPData(s).Trial_data_array(t).target ~= RF_side %if tg is NOT in rf side (if tg is NOT set to toIF, if tg is set to awayIF)
                        givenawayIF{s}(i,c) = givenawayIF{s}(i,c) + 1;  %keep a running count of all give awayIF trials per coherence
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1  % and got it correct awayIF, answered awayIF given awayIF stim
                            awayIFgivenawayIF{s}(i,c) = awayIFgivenawayIF{s}(i,c) + 1;  %keep a running count of all awayIF choices given awayIF stim
                        end
                        givenawayIF{s}(i,c) = givenawayIF{s}(i,c) + 1;  %keep a running count of all give awayIF trials per coherence
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 1  % and got it correct awayIF, answered awayIF given awayIF stim
                            awayIFgivenawayIF{s}(i,c) = awayIFgivenawayIF{s}(i,c) + 1;  %keep a running count of all awayIF choices given awayIF stim
                        end
                        
                    elseif AllGPdata{i,1}.GPData(s).Trial_data_array(t).target == RF_side  %if tg is set toIF
                        giventoIF{s}(i,c) = giventoIF{s}(i,c) + 1;  %total num of toIF set trials
                        if AllGPdata{i,1}.GPData(s).Trial_data_array(t).correct == 0      %got incorrect toIF, answered awayIF given toIF
                            awayIFgiventoIF{s}(i,c) = awayIFgiventoIF{s}(i,c) + 1; %keep a running count of all awayIF choices given toIF
                        end
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
        AllPV_inRF_median_list{s,1} = [AllPV_inRF_median_list{s,1}, median(PV_inRF_list)];
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
        AllPV_outRF_median_list{s,1} = [AllPV_outRF_median_list{s,1}, median(PV_outRF_list)];
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
        if i <= numel(AllGPdata)-3  %for all experiments but the last 3 experiments (had coh1 set of coherences)
            for c_inact = 1:numel(coh1_inact) %for each signed coherence
                AllYInactNum_coh1{s,c_inact} = [AllYInactNum_coh1{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact)]; %number of choices toIF per coherence
                AllYInactTottrials_coh1{s,c_inact} = [AllYInactTottrials_coh1{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact)]; %total number trials per coherence
                AllYInactProp_coh1{s,c_inact} = [AllYInactProp_coh1{s,c_inact}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact))]; %proportion of choices toIF per coherence 
            end
            for c_acc = 1:numel(coh1_acc) %for each unsigned coherence
                AllYAccNum_coh1{s,c_acc} = [AllYAccNum_coh1{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc)]; %number of trials correct per coherence (unsigned)
                AllYAccTottrials_coh1{s,c_acc} = [AllYAccTottrials_coh1{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc)]; %total number of trials per coherence (unsigned)
                AllAccProp_coh1{s,c_acc} = [AllAccProp_coh1{s,c_acc}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc))]; %proportion correct per coherence (unsigned)
            end
            filenames_coh1{i,1} = AllGPdata{i,1}.GPData(s).Filenames;
        elseif i>numel(AllGPdata)-3  %do same thing for the last 3 experiments (were done with coh2 set of coherences)
            for c_inact = 1:numel(coh1_inact)
                AllYInactNum_coh2{s,c_inact} = [AllYInactNum_coh2{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact)];
                AllYInactTottrials_coh2{s,c_inact} = [AllYInactTottrials_coh2{s,c_inact}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact)];
                AllYInactProp_coh2{s,c_inact} = [AllYInactProp_coh2{s,c_inact}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumInact(c_inact))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialInact(c_inact))];
            end
            for c_acc = 1:numel(coh1_acc)
                AllYAccNum_coh2{s,c_acc} = [AllYAccNum_coh2{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc)];
                AllYAccTottrials_coh2{s,c_acc} = [AllYAccTottrials_coh2{s,c_acc}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc)];
                AllAccProp_coh2{s,c_acc} = [AllAccProp_coh2{s,c_acc}, (AllGPdata{i,1}.GPData(s).YAndParams_array.DataYNumAcc(c_acc))/(AllGPdata{i,1}.GPData(s).YAndParams_array.DataYnTrialAcc(c_acc))];
            end
            filenames_coh2{i,1} = AllGPdata{i,1}.GPData(s).Filenames;
        end
        
        %Get Alpha and Beta and Lambda Values for each
        AllAlpha{s} = [AllAlpha{s}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams(1)];
        AllBeta{s} = [AllBeta{s}, AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams(2)];
    end
   
end

%Getting mean RT's per coherence for CMF (chronometric function)
mean_RT_Pre = zeros(9,13); % 9 = number of RT task experiments, 13 = number of coherences
mean_RT_Post = zeros(9,13);
mean_RT_Rec = zeros(9,13);
x = 1;
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
    Ultim_mean_RT_Pre(c) = mean(mean_RT_Pre(find(mean_RT_Pre(:,c) > 0), c));
    Ultim_mean_RT_Post(c) = mean(mean_RT_Post(find(mean_RT_Post(:,c) > 0), c));
    Ultim_mean_RT_Rec(c) = mean(mean_RT_Rec(find(mean_RT_Rec(:,c) > 0), c));
end 
Ultim_mean_RT = [Ultim_mean_RT_Pre; Ultim_mean_RT_Post; Ultim_mean_RT_Rec];

%for computing d' and criterion collapsed across all sessions in SDT
%framework
mu = 0;
sigma = 1;
pd = makedist('Normal',mu,sigma);
for s = 1:conditions %for pre, post, and rex
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
%get rid of all nan values, where the peak velocity was not recorded or
%could not be calculated 
AllPV_inRF_mean_list{1,1}(find(isnan(AllPV_inRF_mean_list{1,1}))) = [];
AllPV_inRF_mean_list{2,1}(find(isnan(AllPV_inRF_mean_list{2,1}))) = [];
AllPV_inRF_mean_list{3,1}(find(isnan(AllPV_inRF_mean_list{3,1}))) = [];
AllPV_outRF_mean_list{1,1}(find(isnan(AllPV_outRF_mean_list{1,1}))) = [];
AllPV_outRF_mean_list{2,1}(find(isnan(AllPV_outRF_mean_list{2,1}))) = [];
AllPV_outRF_mean_list{3,1}(find(isnan(AllPV_outRF_mean_list{3,1}))) = [];

%Get ultimate Means across all experiments
%PV
Ultim_MeanPV_inRF_pre = mean(AllPV_inRF_mean_list{1,1});
Ultim_MeanPV_inRF_post = mean(AllPV_inRF_mean_list{2,1});
Ultim_MeanPV_inRF_rec = mean(AllPV_inRF_mean_list{3,1});
Ultim_MeanPV_outRF_pre = mean(AllPV_outRF_mean_list{1,1});
Ultim_MeanPV_outRF_post = mean(AllPV_outRF_mean_list{2,1});
Ultim_MeanPV_outRF_rec = mean(AllPV_outRF_mean_list{3,1});
%PV for RT experiments
Ultim_MeanPV_inRF_pre_RT = mean(AllPV_inRF_list_RT{1,1});
Ultim_MeanPV_inRF_post_RT = mean(AllPV_inRF_list_RT{2,1});
Ultim_MeanPV_inRF_rec_RT = mean(AllPV_inRF_list_RT{3,1});
Ultim_MeanPV_outRF_pre_RT = mean(AllPV_outRF_list_RT{1,1});
Ultim_MeanPV_outRF_post_RT = mean(AllPV_outRF_list_RT{2,1});
Ultim_MeanPV_outRF_rec_RT = mean(AllPV_outRF_list_RT{3,1});
%RT
Ultim_MeanRT_inRF_pre = mean(AllRT_inRF_list{1,1});
Ultim_MeanRT_inRF_post = mean(AllRT_inRF_list{2,1});
Ultim_MeanRT_inRF_rec = mean(AllRT_inRF_list{3,1});
Ultim_MeanRT_outRF_pre = mean(AllRT_outRF_list{1,1});
Ultim_MeanRT_outRF_post = mean(AllRT_outRF_list{2,1});
Ultim_MeanRT_outRF_rec = mean(AllRT_outRF_list{3,1});
%Accuracy
Ultim_MeanAcc_inRF_pre = mean(AllAcc_inRF_list{1,1});
Ultim_MeanAcc_inRF_post = mean(AllAcc_inRF_list{2,1});
Ultim_MeanAcc_outRF_pre = mean(AllAcc_outRF_list{1,1});
Ultim_MeanAcc_outRF_post = mean(AllAcc_outRF_list{2,1});
%Number of choices toIF (YInactNum), total number of trials
%(YInactTottrials), and proportion of trials chosen toIF (YInactProp),
%separately for 1st coherence set (coh1) and 2nd coherence set (coh2), and
%those for the accuracy stats
for s = 1:conditions
    for c_inact = 1:numel(coh1_inact)
        Ultim_YInactNum_coh1{s,c_inact} = sum(AllYInactNum_coh1{s,c_inact});
        Ultim_YInactNum_coh2{s,c_inact} = sum(AllYInactNum_coh2{s,c_inact});
        Ultim_YInactTottrials_coh1{s,c_inact} = sum(AllYInactTottrials_coh1{s,c_inact});
        Ultim_YInactTottrials_coh2{s,c_inact} = sum(AllYInactTottrials_coh2{s,c_inact});
        Ultim_YInactProp_coh1{s,c_inact} = mean(AllYInactProp_coh1{s,c_inact});
        Ultim_YInactProp_coh2{s,c_inact} = mean(AllYInactProp_coh2{s,c_inact});
    end
    for c_acc = 1:numel(coh1_acc)
        Ultim_YAccNum_coh1{s,c_acc} = sum(AllYAccNum_coh1{s,c_acc});
        Ultim_YAccNum_coh2{s,c_acc} = sum(AllYAccNum_coh2{s,c_acc});
        Ultim_YAccTottrials_coh1{s,c_acc} = sum(AllYAccTottrials_coh1{s,c_acc});
        Ultim_YAccTottrials_coh2{s,c_acc} = sum(AllYAccTottrials_coh2{s,c_acc});
        Ultim_AccProp_coh1{s,c_acc} = mean(AllAccProp_coh1{s,c_acc});
        Ultim_AccProp_coh2{s,c_acc} = mean(AllAccProp_coh2{s,c_acc});
    end
    
end
%% group the data into one combined coherence set, combining the 1st and 2nd set of coherences (total of 15 coherences)
Ultim_YInactNum = zeros(3,15);
Ultim_YInactTottrials = zeros(3,15);

%coh1 = -50, coh2 = -36, coh3 = -24, coh4 = -17, coh5 = -10, coh6 = -5, coh
%7 = -3, coh8 = 0, coh9 = 3, coh10 = 5, coh11 = 10, coh 12 = 17, coh 13 =
%24, coh14 = 36, coh15 = 50
Ultim_YInactNum(:,1) = [Ultim_YInactNum_coh2{:,1}];     %-50 coh
Ultim_YInactTottrials(:,1) = [Ultim_YInactTottrials_coh2{:,1}];     %-50 coh
for c = 2:6
    b = c-1;
    Ultim_YInactNum(:,c) = [Ultim_YInactNum_coh2{:,c}] + [Ultim_YInactNum_coh1{:,b}];
    Ultim_YInactTottrials(:,c) = [Ultim_YInactTottrials_coh2{:,c}] + [Ultim_YInactTottrials_coh1{:,b}];
end
Ultim_YInactNum(:,7) = [Ultim_YInactNum_coh1{:,6}];     %-3 coh
Ultim_YInactTottrials(:,7) = [Ultim_YInactTottrials_coh1{:,6}];
Ultim_YInactNum(:,8) = [Ultim_YInactNum_coh2{:,7}] + [Ultim_YInactNum_coh1{:,7}]; %0 coherence
Ultim_YInactTottrials(:,8) = [Ultim_YInactTottrials_coh2{:,7}] + [Ultim_YInactTottrials_coh1{:,7}];
Ultim_YInactNum(:,9) = [Ultim_YInactNum_coh1{:,8}];     %3 coh
Ultim_YInactTottrials(:,9) = [Ultim_YInactTottrials_coh1{:,8}];
for c = 10:14
    b = c-1;
    a = c-2;
    Ultim_YInactNum(:,c) = [Ultim_YInactNum_coh2{:,a}] + [Ultim_YInactNum_coh1{:,b}];
    Ultim_YInactTottrials(:,c) = [Ultim_YInactTottrials_coh2{:,a}] + [Ultim_YInactTottrials_coh1{:,b}];
end
Ultim_YInactNum(:,15) = [Ultim_YInactNum_coh2{:,13}]; %50 coh
Ultim_YInactTottrials(:,15) = [Ultim_YInactTottrials_coh2{:,13}];

Ultim_YInactProp = Ultim_YInactNum./Ultim_YInactTottrials;

%% Get the logistic fits to the Ultimate mean PMF -- using palamedes toolbox (www.palamedestoolbox.org)
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

x_fit_inact = linspace(cohall_inact(1),cohall_inact(end),100); %coherence values for logistic fit for the mean PMF
x_fit_inact_36 = linspace(coh1_inact(1),coh1_inact(end),100); %coherence values for the logistic fit for the PMF's of experiments with coherence set 1 (starting with 36, etc...)
x_fit_inact_50 = linspace(coh2_inact(1),coh2_inact(end),100); %coherence values for the logistic fit for the PMF's of experiments with coherence set 1 (starting with 50, etc...)

%Fit for every pre post rec for every experiment  
for s = 1:3
    disp('Fitting Ultimate Inact pre/post/rec function.....');
    [paramsValuesUltimInact(s,:) LL exitflag] = PAL_PFML_Fit(cohall_inact,[Ultim_YInactNum(s,:)], ...
        [Ultim_YInactTottrials(s,:)],searchGrid,paramsFree,PF)
    y_HatML_inact(s,:) = (PAL_Logistic(paramsValuesUltimInact(s,:), x_fit_inact)).*100;
end
for i = 1:numel(AllGPdata) %
    for s = 1:3
        %For all PMF fits with coherences (start with 36 etc..)
        if i <= numel(AllGPdata)-3
            AllGPdata{i,1}.GPData(s).YAndParams_array.Y_HatML_Inact = (PAL_Logistic(AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams, x_fit_inact_36)).*100;
        %For all PMF fits with coherences (start with 50 etc...)
        elseif i > numel(AllGPdata)-3
            AllGPdata{i,1}.GPData(s).YAndParams_array.Y_HatML_Inact = (PAL_Logistic(AllGPdata{i,1}.GPData(s).YAndParams_array.DataInactParams, x_fit_inact_50)).*100;
        end
    end
end

%% Fitting a linear fit for RT chronometric function
%
x_fit_inactcoh1_1 = linspace(coh1_inact(1), coh1_inact(7),100); %coherence values for linear fit; _1 is awayIF RT chronometric linear fit
x_fit_inactcoh1_2 = linspace(coh1_inact(7), coh1_inact(end), 100); %coherence values for linear fit; _2 is toIF RT chronometric linear fit 

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
AllRT{1,1} = [];
AllRT{2,1} = [];
AllRT{3,1} = [];
for c = 1:size(AllRT_coh1,2)
    for s = 1:size(AllRT_coh1, 1)
        AllRT{s,1} = [AllRT{s,1};AllRT_coh1{s,c}];
    end
end

linfit_yint_pre = [];
linfit_slope_pre= [];
linfit_yint_post = [];
linfit_slope_post = [];
linfit_yint_rec = [];
linfit_slope_rec = [];

%fit each individual RT task experiment
for i = 1:9
    for s = 1:3
        %get RT's
        RT{s,1} = [];
        RT{s,1} = AllRT{s,1}(:,i)';
        %get linear fits
        x_1 = [ones(length(coh1_inact(1:7)),1),coh1_inact(1:7)'];
        x_2 = [ones(length(coh1_inact(7:end)),1),coh1_inact(7:end)'];
        b_1(s,:) = (x_1\RT{s,1}(1:7)')'; %intercept and slope; _1 is awayIF
        b_2(s,:) = (x_2\RT{s,1}(7:end)')'; %intercept and slope; _2 is toIF
        y_linfit_1(s,:) = b_1(s,1) + b_1(s,2)*x_fit_inactcoh1_1; %y values for the linear fit; _1 is awayIF
        y_linfit_2(s,:) = b_2(s,1) + b_2(s,2)*x_fit_inactcoh1_2; %y values for the linear fit; _2 is toIF
    end
    y_linefit_all{i,1} = y_linfit_1;
    y_linefit_all{i,2} = y_linfit_2;
    linfit_yint_pre = [linfit_yint_pre; b_1(1,1), b_2(1,1)];            %RT intercept paramater - 1st col is awayIF, 2nd col is toIF
    linfit_slope_pre  = [linfit_slope_pre; b_1(1,2), b_2(1,2)];         %RT slope paramater - 1st col is awayIF, 2nd col is toIF
    linfit_yint_post = [linfit_yint_post; b_1(2,1), b_2(2,1)];
    linfit_slope_post = [linfit_slope_post; b_1(2,2), b_2(2,2)];
    linfit_yint_rec = [ linfit_yint_rec; b_1(3,1), b_2(3,1)];
    linfit_slope_rec = [linfit_slope_rec; b_1(3,2), b_2(3,2)];
end

%Put the slope and intercept RT parameters for pre, post, and rec in a
%data structure
RT_linfit_params = struct('Slope', [], 'Int', []);
RT_linfit_params(1).Slope = linfit_slope_pre;
RT_linfit_params(2).Slope = linfit_slope_post;
RT_linfit_params(3).Slope = linfit_slope_rec;
RT_linfit_params(1).Int = linfit_yint_pre;
RT_linfit_params(2).Int = linfit_yint_post;
RT_linfit_params(3).Int = linfit_yint_rec;

%fit the mean RT liner fit
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
    marker = 'o';
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
hold on
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
plot(x_fit_inact, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact, y_HatML_inact(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(cohall_inact, Ultim_YInactProp(1,:)*100, dotsz, pre_color_mean, 'filled');
scatter(cohall_inact, Ultim_YInactProp(2,:)*100, dotsz, post_color_mean,  'filled');
for i = 1:numel(AllGPdata)
    plot(AllGPdata{i,1}.GPData(1).YAndParams_array.xfit, AllGPdata{i,1}.GPData(1).YAndParams_array.Y_HatML_Inact , 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(AllGPdata{i,1}.GPData(2).YAndParams_array.xfit, AllGPdata{i,1}.GPData(2).YAndParams_array.Y_HatML_Inact , 'color', post_color_indiv, 'Linewidth', Linewidth_1);
    scatter(AllGPdata{i,1}.GPData(1).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(1).YAndParams_array.DataYPropInact, dotsz, pre_color_indiv,  'filled');
    scatter(AllGPdata{i,1}.GPData(2).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(2).YAndParams_array.DataYPropInact, dotsz, post_color_indiv,  'filled');
    pause = 1;
end
plot(x_fit_inact, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact, y_HatML_inact(2,:), 'color', post_color_mean, 'Linewidth', Linewidth);
scatter(cohall_inact, Ultim_YInactProp(1,:)*100, dotsz, pre_color_mean, 'filled');
scatter(cohall_inact, Ultim_YInactProp(2,:)*100, dotsz, post_color_mean,  'filled');
% fiducials
s = plot([0 0],[0 100], ':k');
m = plot([-50, 50], [50 50], ':k');
ticks = [0:10:100];
ylabel('Choices to IF (%)', 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold');
ylim([0 100]);
title('Pre VS Post','FontWeight', 'bold');
fsize=30;
ylim([0 100]);
set(gca, 'ytick', [0, 20,40, 60, 80, 100], 'yticklabels', [0, 20, 40, 60, 80, 100]);

%% PMF Plot -- PRE vs REC
figure(6);
set(gcf,'color','w');
hold on
plot(x_fit_inact, y_HatML_inact(1,:), 'color',pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact, y_HatML_inact(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(cohall_inact, Ultim_YInactProp(1,:)*100, dotsz, pre_color_mean, 'filled');
scatter(cohall_inact, Ultim_YInactProp(3,:)*100, dotsz, rec_color_mean,  'filled');
for i = 1:numel(AllGPdata)
    plot(AllGPdata{i,1}.GPData(1).YAndParams_array.xfit, AllGPdata{i,1}.GPData(1).YAndParams_array.Y_HatML_Inact , 'color', pre_color_indiv, 'Linewidth', Linewidth_1);
    plot(AllGPdata{i,1}.GPData(3).YAndParams_array.xfit, AllGPdata{i,1}.GPData(3).YAndParams_array.Y_HatML_Inact , 'color', rec_color_indiv, 'Linewidth', Linewidth_1);
    scatter(AllGPdata{i,1}.GPData(1).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(1).YAndParams_array.DataYPropInact, dotsz, pre_color_indiv,  'filled');
    scatter(AllGPdata{i,1}.GPData(3).YAndParams_array.x_coh, AllGPdata{i,1}.GPData(3).YAndParams_array.DataYPropInact, dotsz, rec_color_indiv,  'filled');
    pause = 1;
end
plot(x_fit_inact, y_HatML_inact(1,:), 'color', pre_color_mean, 'Linewidth', Linewidth);
plot(x_fit_inact, y_HatML_inact(3,:), 'color', rec_color_mean, 'Linewidth', Linewidth);
scatter(cohall_inact, Ultim_YInactProp(1,:)*100, dotsz, pre_color_mean, 'filled');
scatter(cohall_inact, Ultim_YInactProp(3,:)*100, dotsz, rec_color_mean,  'filled');
s = plot([0 0],[0 100], ':k');
m = plot([-40, 40], [50 50], ':k');
ticks = [0:10:100];
ylabel('Choices to IF (%)', 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold');
ylim([0 100]);
title('Pre VS Rec','FontWeight', 'bold');
fsize=30;
ylim([0 100]);
set(gca, 'ytick', [0, 20,40, 60, 80, 100], 'yticklabels', [0, 20, 40, 60, 80, 100]);

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

for i = 1:9
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

s = plot([0 0],[0 100], ':k');
set(gca, 'xtick', coh1_inact, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ticks = [0:10:100];
ylabel('RT', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)','FontWeight', 'bold', 'Fontsize', 30);
ylim([0.4 1.3]);
title('Pre VS Post','FontWeight', 'bold', 'Fontsize', 30);
fsize=20;
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
for i = 1:9
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
scatter(median(AllAlpha{1}), median(AllAlpha{3}), dotsz, rec_color_mean, 'filled')
ylim([-10 30])
xlim([-10 30])
plot([-10:30],[-10:30], ':k');
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Alpha','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post', 'Rec'}, 'Fontsize', 30);
%% Beta Plot pre vs post vs rec plot
figure(10);
scatter(AllBeta{1}, AllBeta{2}, dotsz,post_color_indiv, 'filled')
hold on
scatter(AllBeta{1}, AllBeta{3}, dotsz,rec_color_indiv, 'filled')
scatter(median(AllBeta{1}), median(AllBeta{2}), dotsz,post_color_mean, 'filled')
scatter(median(AllBeta{1}), median(AllBeta{3}), dotsz,rec_color_mean, 'filled')
ylim([0.05 0.3])
xlim([0.05 0.3])
% fiducials
plot([0.05:0.01:0.3],[0.05:0.01:0.3], ':k');
% data
% set axis and figure properties
ylabel('Post and Rec', 'Fontsize', 30, 'FontWeight', 'bold'); %  correct/total per coherence
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
% fig=gcf;
% set(findall(fig,'-property','fontsize'),'fontsize',fsize);
title('Beta','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post', 'Rec'}, 'Fontsize', 30);

%% RT slope plot - toIF and awayIF - pre vs post
figure(14);
scatter(linfit_slope_pre(:,1), linfit_slope_post(:,1), dotsz, awayIF_color_indiv, 'filled');
hold on
scatter(linfit_slope_pre(:,2), linfit_slope_post(:,2),dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_slope_pre(:,1)), median(linfit_slope_post(:,1)), dotsz, awayIF_color_mean, 'filled');
scatter(median(linfit_slope_pre(:,2)), median(linfit_slope_post(:,2)),dotsz, toIF_color_mean, 'filled');

plot([-0.008 0.008 ], [-0.008 0.008 ], ':k')
plot([0, 0 ], [-0.008, 0.008] ,':k')
plot([-0.008, 0.008], [0, 0 ], ':k')
ylim([-0.008 , 0.008 ])
xlim([-0.008 , 0.008 ])

ylabel('Post', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Slope of RT linear fit - Pre vs Post','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);


%% RT slope plot - toIF and awayIF - pre vs rec
figure(36);
hold on
scatter(linfit_slope_pre(:,1), linfit_slope_rec(:,1),dotsz,awayIF_color_indiv, 'filled');
scatter(linfit_slope_pre(:,2), linfit_slope_rec(:,2),dotsz,toIF_color_indiv, 'filled');

scatter(median(linfit_slope_pre(:,1)), median(linfit_slope_rec(:,1)),dotsz,awayIF_color_mean, 'filled');
scatter(median(linfit_slope_pre(:,2)), median(linfit_slope_rec(:,2)),dotsz,toIF_color_mean, 'filled');

plot([-0.008 0.008 ], [-0.008 0.008 ], ':k')
plot([0, 0 ], [-0.008, 0.008] ,':k')
plot([-0.008, 0.008], [0, 0 ], ':k')
ylim([-0.008 , 0.008 ])
xlim([-0.008 , 0.008 ])

ylabel('Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Slope of RT linear fit - Pre vs Rec','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);
%% RT intercept plot - toIF and awayIF - pre vs post
figure(15);
scatter(linfit_yint_pre(:,1), linfit_yint_post(:,1), dotsz, awayIF_color_indiv, 'filled');
hold on
scatter(linfit_yint_pre(:,2), linfit_yint_post(:,2), dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_yint_pre(:,1)), median(linfit_yint_post(:,1)), dotsz, awayIF_color_mean, 'filled');
scatter(median(linfit_yint_pre(:,2)), median(linfit_yint_post(:,2)), dotsz, toIF_color_mean, 'filled');

plot([0.5 , 1.2 ], [0.5, 1.2 ], ':k')
ylim([0.5 , 1.2 ])
xlim([0.5 , 1.2 ])
ylabel('Post', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Intercept of RT linear fit - Pre vs Post','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);
%% RT intercept plot - toIF and awayIF - pre vs  rec 
figure(37);
hold on
scatter(linfit_yint_pre(:,1), linfit_yint_rec(:,1), dotsz,awayIF_color_indiv, 'filled');
scatter(linfit_yint_pre(:,2), linfit_yint_rec(:,2),dotsz, toIF_color_indiv, 'filled');

scatter(median(linfit_yint_pre(:,1)), median(linfit_yint_rec(:,1)), dotsz, awayIF_color_mean, 'filled');
scatter(median(linfit_yint_pre(:,2)), median(linfit_yint_rec(:,2)), dotsz,toIF_color_mean, 'filled');

plot([0.5 , 1.2 ], [0.5, 1.2 ], ':k')
ylim([0.5 , 1.2 ])
xlim([0.5 , 1.2 ])
ylabel('Rec', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Pre','FontWeight', 'bold', 'Fontsize', 30);
title('Intercept of RT linear fit - Pre vs Rec','FontWeight', 'bold', 'Fontsize', 30);
legend({'Post outIF', 'Rec outIF', 'Post inIF', 'Rec inIF' }, 'Fontsize', 30);
%% Plot D' and C in SDT -- separated by coherences and collapsed over coherence 
C = {pre_color_mean,post_color_mean, rec_color_mean};
figure(33);
hold on
for s = 1:conditions
    plot(coh_dc, d_prime{s}, 'Linewidth', Linewidth, 'color', C{s})
    scatter(coh_dc, d_prime{s}, dotsz, C{s},'filled')
    xlim([ 0, 53]);
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
    xlim([ 0, 53]);
end
xticks(coh_dc);
ylabel('Criterion', 'Fontsize', 30, 'FontWeight', 'bold');
xlabel('Stimulus Difficulty','FontWeight', 'bold', 'Fontsize', 30);
title('Criterion Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
% stop = 1;

clear x
x{1} = zeros(1,size(d_prime_sessions{1},2)) +1;
x{2} = zeros(1,size(d_prime_sessions{1},2)) +2;
x{3} = zeros(1,size(d_prime_sessions{1},2)) +3;

figure(35)
hold on
for s = 1:conditions
    bar(s, mean(d_prime_sessions{s}),'Facecolor',C{s});
end
for i = 1:size(d_prime_sessions{1},2)
    plot([1,2,3], [d_prime_sessions{1}(i), d_prime_sessions{2}(i), d_prime_sessions{3}(i)], '-o');
end
ylabel('D-Prime', 'Fontsize', 30, 'FontWeight', 'bold');
xticks([1,2,3])
xticklabels({'Pre', 'Post', 'Rec'});
title('D-Prime Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
legend('Pre', 'Post', 'Rec', 'Fontsize', 30);
stop = 1;

figure(38)
hold on
for s = 1:conditions
    bar(s, mean(criterion_sessions{s}),'Facecolor',C{s});
    %     ylim([-0.7 0.25])
end
for i = 1:size(criterion_sessions{1},2)
    plot([1,2,3], [criterion_sessions{1}(i), criterion_sessions{2}(i), criterion_sessions{3}(i)], '-o');
end
ylabel('Criterion', 'Fontsize', 30, 'FontWeight', 'bold');
% xlabel('Stimulus Difficulty','FontWeight', 'bold', 'Fontsize', 30);
xticks([1,2,3])
xticklabels({'Pre', 'Post', 'Rec'});
title('Criterion Pre VS Post VS Rec','FontWeight', 'bold', 'Fontsize', 30);
% legend('Pre', 'Post', 'Rec', 'Fontsize', 30);
stop = 1;


