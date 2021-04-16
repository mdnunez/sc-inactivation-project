% Script to use once all the individual session data have been processed.
% This script analyzes/plots all the injection session data together 
clc;
clear all;

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

prev_load = input('Open Previously loaded data?', 's');
%Open previously loaded data?
if prev_load == 'y'
    clc
    [fileName, pathname] = getFile('*.*');
    file = [pathname fileName];
    load(file)
elseif prev_load == 'n'
    %Load Data
    while(continueAction)
        clc
        disp('Open LM Data')
        [fileName, pathname] = getFile('*.*');
        file = [pathname fileName];
        AllLMdata{i,1} = load(file);
        AllLMdata{i,2} = fileName;
        
        if input('Open one more file? y/n', 's') == 'y'
            files_inputted = files_inputted +1;
            continueAction = 1;
            i = i+1;
        else
            continueAction = 0;
        end
    end
    
end
%%

conditions = 3; %pre, post, rec
%make empty data arrays to fill in (PV = peak velocity, RT = reaction time)
for i = 1:conditions
    AllPV_inRF_mean_list{i,1} = [];
    AllPV_inRF_mean1_list{i,1} = [];
    AllPV_inRF_median_list{i,1} = [];
    AllPV_outRF_mean_list{i,1} = [];
    AllPV_outRF_mean1_list{i,1} = [];
    AllPV_outRF_median_list{i,1} = [];
    
    AllRT_inRF_mean_list{i,1} = [];
    AllRT_outRF_mean_list{i,1} = [];
    
end

for i = 1:conditions
    AllPV_inRF_list{i,1} = [];
    AllPV_outRF_list{i,1} = [];
end

%% Extracting out data from all injections
for i = 1:size(AllLMdata,1)        %for each experiment
    
    currentRFx = AllLMdata{i,1}.LMData.RFx; %get current horizontal degree position of RF
    currentRFy = AllLMdata{i,1}.LMData.RFy; %get current vertical degree position of RF
    currentcontraRFx = -(AllLMdata{i,1}.LMData.RFx); %same for contralateral to the RF
    currentcontraRFy = AllLMdata{i,1}.LMData.RFy;
    
    for s = 1:conditions                   %for each session: pre, post, rec
        
        PV_inRF_list = [];
        PV_outRF_list = [];
        RT_inRF_list = [];
        RT_outRF_list = [];
        
        for t = 1:numel(AllLMdata{i,1}.LMData.Trial_3{s,1})   %for each trial
            
            %separate peak velocities and RT by toIF and awayIF saccades
            %toIF
            if AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RFx == currentRFx && AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RFy == currentRFy      
                %if target is inRF (toIF)
                %Get the peak velocity 
                if ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                    if ~isempty(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                        if AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity < 1600 && ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                            PV_inRF_list = [PV_inRF_list, AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity];
                        end
                    end
                end
                %Get the reaction time
                if ~isempty(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT)
                    if AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT > 0.1 && ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT)
                        RT_inRF_list = [RT_inRF_list, AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT];
                    end
                end
            elseif AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RFx == currentcontraRFx && AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RFy == currentcontraRFy  
                %if target is outRF (awayIF)
                %Get Peak velocity
                if ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                    if ~isempty(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                        if AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity < 1600 && ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity)
                            PV_outRF_list = [PV_outRF_list, AllLMdata{i,1}.LMData.Trial_3{s,1}(t).PeakVelocity];
                        end
                    end
                end
                %Get RT
                if AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT > 0.1 && ~isempty(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT)
                    if ~isnan(AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT)
                        RT_outRF_list = [RT_outRF_list, AllLMdata{i,1}.LMData.Trial_3{s,1}(t).RT];
                    end
                end
            end
            
            
        end
        
        
        %getting a list of all the data across experiments, all PV values
        %separated by pre, post, and rec, but collapsed over experiments
        AllPV_inRF_list{s,1} = [AllPV_inRF_list{s,1}, PV_inRF_list]; %all PV inRF values
        AllPV_outRF_list{s,1} = [AllPV_outRF_list{s,1}, PV_outRF_list]; %allPV outRF values 
        AllPV_inRF_mean_list{s,1} = [AllPV_inRF_mean_list{s,1}, AllLMdata{i,1}.LMData.MeanPVInRF{s,1}];    %mean of the PV inRF for each session
        AllPV_outRF_mean_list{s,1} = [AllPV_outRF_mean_list{s,1}, AllLMdata{i,1}.LMData.MeanPVOutRF{s,1}]; %mean of the PV outRF for each session
        AllRT_inRF_mean_list{s,1} = [AllRT_inRF_mean_list{s,1},  AllLMdata{i,1}.LMData.MeanRTInRF{s,1}];    %same for RT (mean RT)
        AllRT_outRF_mean_list{s,1} = [AllRT_outRF_mean_list{s,1}, AllLMdata{i,1}.LMData.MeanRTOutRF{s,1}];
        
        AllPV_inRF_median_list{s,1} = [AllPV_inRF_median_list{s,1}, AllLMdata{i,1}.LMData.MedianPVInRF{s,1}];    %median of the PV inRF for each session
        AllPV_outRF_median_list{s,1} = [AllPV_outRF_median_list{s,1}, AllLMdata{i,1}.LMData.MedianPVOutRF{s,1}]; %median of the PV outRF for each session
        
    end
    
end
%%
%Ultimate Means across all experiments
Ultim_MeanPV_inRF_pre = mean(AllPV_inRF_mean_list{1,1});
Ultim_MeanPV_inRF_post = mean(AllPV_inRF_mean_list{2,1});
Ultim_MeanPV_inRF_rec = mean(AllPV_inRF_mean_list{3,1});
Ultim_MeanPV_outRF_pre = mean(AllPV_outRF_mean_list{1,1});
Ultim_MeanPV_outRF_post = mean(AllPV_outRF_mean_list{2,1});
Ultim_MeanPV_outRF_rec = mean(AllPV_outRF_mean_list{3,1});

Ultim_Mean1PV_inRF_pre = mean(AllPV_inRF_mean1_list{1,1});
Ultim_Mean1PV_inRF_post = mean(AllPV_inRF_mean1_list{2,1});
Ultim_Mean1PV_inRF_rec = mean(AllPV_inRF_mean1_list{3,1});
Ultim_Mean1PV_outRF_pre = mean(AllPV_outRF_mean1_list{1,1});
Ultim_Mean1PV_outRF_post = mean(AllPV_outRF_mean1_list{2,1});
Ultim_Mean1PV_outRF_rec = mean(AllPV_outRF_mean1_list{3,1});

Ultim_MedianPV_inRF_pre = median(AllPV_inRF_median_list{1,1});
Ultim_MedianPV_inRF_post = median(AllPV_inRF_median_list{2,1});
Ultim_MedianPV_inRF_rec = median(AllPV_inRF_median_list{3,1});
Ultim_MedianPV_outRF_pre = median(AllPV_outRF_median_list{1,1});
Ultim_MedianPV_outRF_post = median(AllPV_outRF_median_list{2,1});
Ultim_MedianPV_outRF_rec = median(AllPV_outRF_median_list{3,1});

Ultim_MeanRT_inRF_pre = mean(AllRT_inRF_mean_list{1,1});
Ultim_MeanRT_inRF_post = mean(AllRT_inRF_mean_list{2,1});
Ultim_MeanRT_inRF_rec = mean(AllRT_inRF_mean_list{3,1});
Ultim_MeanRT_outRF_pre = mean(AllRT_outRF_mean_list{1,1});
Ultim_MeanRT_outRF_post = mean(AllRT_outRF_mean_list{2,1});
Ultim_MeanRT_outRF_rec = mean(AllRT_outRF_mean_list{3,1});

%%
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
else
    marker = 'o';
end



%% Plot Peak Velocity inRF (toIF) vs OutRF (awayIF) for all injections
x1 = zeros(1,numel(AllPV_inRF_mean_list{1,1}))+1;
x2 = x1 +1;
%Plotter for PV inRF vs OutRF
figure(37);
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