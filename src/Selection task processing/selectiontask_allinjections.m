%Script to plot and analyze all the injection data together for the saccade
%selection task
%% Loading Data
close all 
clear all
continueAction = 1;
files_inputted = 1;
i = 1;

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
        disp('Open Selection Task Data')
        [fileName, pathname] = getFile('*.*');
        file = [pathname fileName];
        AllSTOdata{i,1} = load(file);
        AllSTOdata{i,2} = fileName;
        
        if input('Open one more file? y/n', 's') == 'y'
            files_inputted = files_inputted +1;
            continueAction = 1;
            i = i+1;
        else
            continueAction = 0;
        end
    end
    
end

if input('Is this monkey S data?', 's') == 'y'
    SP = 1;
    BT = 0;
else
    SP = 0;
    BT = 1;
end 
conditions = 3; %for pre, post, and rec

%%
%PV = peak velocity; RT = reaction time; outRF = awayIF; inRF = toIF
%create empty arrays to fill in
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
    
    AllRT_inRF_mean_list{i,1} = [];
    AllRT_inRF_std_list{i,1} = [];
    AllRT_outRF_mean_list{i,1} = [];
    AllRT_outRF_std_list{i,1} = [];
    AllRT_inRF_list{i,1} = [];
    AllRT_outRF_list{i,1} = [];
    
    AllAcc_inRF_list{i,1} = [];
    AllAcc_outRF_list{i,1} = [];
end


for i = 1:size(AllSTOdata,1)        %for each experiment
    
    for s = 1:conditions                   %for each session
        PV_inRF_list = [];
        PV_outRF_list = [];
        RT_inRF_list = [];
        RT_outRF_list = [];
        
        for t = 1:numel(AllSTOdata{i,1}.STOData(s).Trial_1{1,1})   %for each trial
            
            if (AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).TGinRF == 1 && AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).outcome == 1) || (AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).TGinRF ~= 1 && AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).outcome == 0)   %if choice is inRF (toIF)
                %Get Peak velocity
                if ~isempty(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity)
                    if AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity < 1600 && ~isnan(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity)
                        PV_inRF_list = [PV_inRF_list, AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity];
                    end
                end
                %Get RT
                if ~isempty(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT)
                    if ~isnan(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT)
                        RT_inRF_list = [RT_inRF_list, AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT];
                    end
                end
            elseif(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).TGinRF ~= 1 && AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).outcome == 1) || (AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).TGinRF == 1 && AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).outcome == 0)   %if choice is outRF (awayIF)
                %Get Peak velocity
                if ~isempty(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity)
                    if AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity < 1600 && ~isnan(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity)
                        PV_outRF_list = [PV_outRF_list, AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).PeakVelocity];
                    end
                end
                %Get RT
                if ~isempty(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT)
                    if ~isnan(AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT)
                        RT_outRF_list = [RT_outRF_list, AllSTOdata{i,1}.STOData(s).Trial_1{1,1}(t).RT];
                    end
                end
            end
        end
        
        AllPV_inRF_list{s,1} = [AllPV_inRF_list{s,1}, PV_inRF_list];        %separated into pre, post, rec
        AllPV_inRF_array{s,i} = PV_inRF_list;                               %separated into pre, post, rec, for each experiment

        AllPV_inRF_mean_list{s,1} = [AllPV_inRF_mean_list{s,1}, mean(PV_inRF_list)];    %separated into pre,post,rec
        AllPV_inRF_median_list{s,1} = [AllPV_inRF_median_list{s,1}, median(PV_inRF_list)];
        AllPV_inRF_std_list{s,1} = [AllPV_inRF_std_list{s,1}, std(PV_inRF_list)];
        AllPV_outRF_list{s,1} = [AllPV_outRF_list{s,1}, PV_outRF_list];
        AllPV_outRF_array{s,i} = PV_outRF_list;
        AllPV_outRF_mean_list{s,1} = [AllPV_outRF_mean_list{s,1}, mean(PV_outRF_list)];
        AllPV_outRF_median_list{s,1} = [AllPV_outRF_median_list{s,1}, median(PV_outRF_list)];
        AllPV_outRF_std_list{s,1} = [AllPV_outRF_std_list{s,1}, std(PV_outRF_list)];

        
        AllRT_inRF_list{s,1} = [AllRT_inRF_list{s,1}, RT_inRF_list];        %separated into pre, post, rec
        AllRT_inRF_array{s,i} = RT_inRF_list;                               %separated into pre, post, rec, for each experiment
        AllRT_inRF_mean_list{s,1} = [AllRT_inRF_mean_list{s,1}, mean(RT_inRF_list)];    %separated into pre,post,rec
        AllRT_inRF_std_list{s,1} = [AllRT_inRF_std_list{s,1}, std(RT_inRF_list)];
        AllRT_outRF_list{s,1} = [AllRT_outRF_list{s,1}, RT_outRF_list];
        AllRT_outRF_array{s,i} = RT_outRF_list;
        AllRT_outRF_mean_list{s,1} = [AllRT_outRF_mean_list{s,1}, mean(RT_outRF_list)];
        AllRT_outRF_std_list{s,1} = [AllRT_outRF_std_list{s,1}, std(RT_outRF_list)];

        AllAcc_outRF_list{s,1} =  [AllAcc_outRF_list{s,1}, AllSTOdata{i,1}.STOData(s).Acc_array{1,2}];
        AllAcc_inRF_list{s,1} =  [AllAcc_inRF_list{s,1}, AllSTOdata{i,1}.STOData(s).Acc_array{2,2}];
        
    end 

end 

original_AllPV_inRF_mean_list = AllPV_inRF_mean_list;
original_AllPV_outRF_mean_list = AllPV_outRF_mean_list;

%%

Ultim_Mean2PV_inRF_pre = mean(AllPV_inRF_mean_list{1,1});
Ultim_Mean2PV_inRF_post = mean(AllPV_inRF_mean_list{2,1});
Ultim_Mean2PV_inRF_rec = mean(AllPV_inRF_mean_list{3,1});
Ultim_Mean2PV_outRF_pre = mean(AllPV_outRF_mean_list{1,1});
Ultim_Mean2PV_outRF_post = mean(AllPV_outRF_mean_list{2,1});
Ultim_Mean2PV_outRF_rec = mean(AllPV_outRF_mean_list{3,1});

Ultim_MeanAcc_inRF_pre = mean(AllAcc_inRF_list{1,1});
Ultim_MeanAcc_inRF_post = mean(AllAcc_inRF_list{2,1});
Ultim_MeanAcc_inRF_rec = mean(AllAcc_inRF_list{3,1});
Ultim_MeanAcc_outRF_pre = mean(AllAcc_outRF_list{1,1});
Ultim_MeanAcc_outRF_post = mean(AllAcc_outRF_list{2,1});
Ultim_MeanAcc_outRF_rec = mean(AllAcc_outRF_list{3,1});

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
end
%%
x1 = zeros(1,numel(AllPV_inRF_mean_list{1,1}))+1;
x2 = x1 +1;
x3 = zeros(1,numel(AllPV_inRF_mean1_list{1,1})) +1;
x4 = x3 +1;

%% Peak velocity plot
x1 = zeros(1,numel(AllPV_inRF_mean_list{1,1}))+1;
x2 = x1 +1;
%Plotter for PV inRF vs OutRF
figure(1);
hold on
scatter(x1, AllPV_inRF_mean_list{1,1}, 50, toIF_color_indiv, marker,'filled');
scatter(x1, AllPV_outRF_mean_list{1,1},50, awayIF_color_indiv, marker,'filled');
scatter(x2, AllPV_inRF_mean_list{2,1},50, toIF_color_indiv,marker, 'filled');
scatter(x2, AllPV_outRF_mean_list{2,1},50, awayIF_color_indiv,marker,'filled');
scatter(1,Ultim_Mean2PV_inRF_pre, 100, toIF_color_mean,marker, 'filled')
scatter(1,Ultim_Mean2PV_outRF_pre, 100, awayIF_color_mean,marker, 'filled')
scatter(2, Ultim_Mean2PV_inRF_post, 100, toIF_color_mean,marker, 'filled');
scatter(2, Ultim_Mean2PV_outRF_post, 100, awayIF_color_mean,marker, 'filled');
plot([1,2], [Ultim_Mean2PV_inRF_pre, Ultim_Mean2PV_inRF_post], marker, 'color', toIF_color_mean, 'linestyle', '-','linewidth', Linewidth_1);
plot([1,2], [Ultim_Mean2PV_outRF_pre, Ultim_Mean2PV_outRF_post],marker, 'color', awayIF_color_mean,'linestyle', '-', 'linewidth', Linewidth_1);
set(gca, 'xtick', [1,2], 'xticklabel', {'Pre', 'Post'}, 'tickdir','out', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('PV', 'Fontsize', 30, 'FontWeight', 'bold');
xlim([0.7, 2.3]);
ylim([200, 1000])
legend('toIF', 'awayIF', 'location', 'south');
hold off;

%% Accuracy plot
x1 = zeros(1,numel(AllAcc_inRF_list{1,1}))+1;
x2 = x1 +1;
%Plotter for Acc inRF vs OutRF
figure(2);
plot([1,1], [Ultim_MeanAcc_inRF_pre,Ultim_MeanAcc_outRF_pre], 'linewidth', 2, 'color', pre_color_mean);
hold on
plot([2,2], [Ultim_MeanAcc_inRF_post,Ultim_MeanAcc_outRF_post], 'linewidth', 2, 'color', pre_color_mean);
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