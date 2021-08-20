clc;
clear all;

%File loading variables
continueAction = 1;
files_inputted = 1;
i = 1;

%Color/marker variables for plotting
C = {'r','b','g','d','y'}; % Cell array of colors for the line plotted.
M = {'o','s','d','h','p'}; % Cell array of markers shape.


% Coherence choice prompt
prompt = ['select coherences:'...
    '\n 1) coh [36 24 17 10 5 3 0] coherences'...
    '\n 2) coh [50 36 24 17 10 5 0] coherences'];

coh_choice = input(prompt);

if input('Delay or RT task? (type "d" or "r")', 's') == 'd'
    DT = 1;
    RT = 0;
else
    RT = 1;
    DT = 0;
end 

rfx = input('Input RFx: ');
rfy = input('Input RFy: ');
if rfx < 0          %If RF is negative, on the left side, then it's a right SC inactivation
    RSC = 1;
    LSC = 0;
else 
    LSC = 1;
    RSC = 0;
end 

a = input('Is this pre or post? 1 for pre, 2 for post, 3 for rec: ');
if a == 1
    pre = 1;
    post = 0;
    rec = 0;
elseif a == 2
    pre = 0;
    post = 1;
    rec = 0;
elseif a == 3
    pre = 0;
    post = 0;
    rec = 1;
end 



% data arrays for the individual files analyzed.
while(continueAction)       
    [Trial_4,coh,filename,y_avg_med_RT, RTStruct,YAndParams] = decisiontask_rawdataprocessor(coh_choice,files_inputted, RSC, LSC, RT, DT);
    i;
    Trial_data_array(1,i) = {Trial_4};
    y_avg_med_RT_array(1,i) = {y_avg_med_RT};
    RTStruct_array(1,i) = {RTStruct};
    filename_array(1,i) = {filename};
    YAndParams_array(1,i) = {YAndParams};
    i = i + 1;
    
    display('Last opened file was');
    filename
    prompt_0 = ('Open another more session (either pre, post, or rec)? (y/n)');
    add_file = input(prompt_0, 's');
    if add_file == 'Y' || add_file == 'y'
        continueAction = 1;
        files_inputted = files_inputted +1;
    elseif add_file == 'N' || add_file == 'n'
        continueAction = 0;
    end
end

%%
%------------ Plot PMF Fitted Data ------------%
close all
figure(1);
set(gcf,'color','w'); 
clf;
hold on;

% fiducials
s = plot([0 0],[0 100], ':k');
m = plot([YAndParams_array{1,1}.x_coh(1), YAndParams_array{1,1}.x_coh(end)], [50 50], ':k');

% fit
for z = 1:files_inputted
plot(YAndParams_array{1,z}.xfit, YAndParams_array{1,z}.Y_HatML_YRight,'color', C{z}, 'linewidth',1);
end 
for z = 1:files_inputted
plot(YAndParams_array{1,z}.x_coh,YAndParams_array{1,z}.DataYPropR,M{z},...
    'markersize',8, ...
    'markeredgecolor','k',...
    'markerfacecolor',C{z});   %color of dots 
end 


% set legend
    if files_inputted == 3
        set(get(get(s(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(m(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('Location','best','Pre Injection', 'Post Injection', 'Recovery')
    elseif files_inputted == 2
        set(get(get(s(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(m(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('Location','best','Pre Injection', 'Post Injection')
    elseif files_inputted == 1
    end

% set axis and figure properties
xlim(1.05*[YAndParams_array{1,z}.x_coh(1), YAndParams_array{1,z}.x_coh(end)]);
set(gca, 'xtick', sort(unique(YAndParams_array{1,1}.x_coh)));
ticks = 0:10:100;
set(gca,'tickdir','out');
ylim([0 105]);
ylabel('Fraction of choices to the right(%)'); %  correct/total per coherence
hold off;
xlabel('Coherence (%)');
fig=gcf;
fsize=12;
set(findall(fig,'-property','fontsize'),'fontsize',fsize);
% text(7,25,{'Alpha/Beta Pre (RED)' num2str(YAndParams_array{1,1}.DataRParams(1:2))}, 'FontSize', 8);
% text(7,15,{'Alpha/Beta Post (BLUE)' num2str(YAndParams_array{1,2}.DataRParams(1:2))}, 'FontSize', 8);
% text(7,5,{'Alpha/Beta Rec (GREEN)' num2str(YAndParams_array{1,3}.DataRParams(1:2))}, 'FontSize', 8);

%% Plot RTs
if RT == 1
        lb = 0.4;
        ub = 2;
     
        
        C = {'g', 'b'};
        E = {'r', 'c'};
        %Plotting RT mean curve
        figure(3);
        hold on
        % data
        for y = 1:files_inputted
            scatter(YAndParams_array{1,1}.x_coh, [y_avg_med_RT_array{1, 1}{2,1}], 10, C{y}, 'filled', 'o');
            scatter(YAndParams_array{1,1}.x_coh, [y_avg_med_RT_array{1, 1}{3,1}], 10, E{y}, 'filled', 'o');
        end
        ylim([lb, ub])
        xlim(1.05*[YAndParams_array{1,1}.x_coh(1), YAndParams_array{1,1}.x_coh(end)]);
        set(gca, 'xtick', sort(unique(YAndParams_array{1,1}.x_coh)));
        ylabel('RT'); %  correct/total per coherence
        xlabel('Coherence (%)');
        title('Mean RT - ALL TRIALS');
        
        hold on
        plot([0 0],[0 1.5], ':k');
        legend('Correct RT_1', 'Incorrect RT_1');
        stop = 1;
end
%% Reformat data for modeling script

    datastruct(:,1) = [Trial_4(1:end).coherence]';
    datastruct(:,2) = [Trial_4(1:end).target]';
    for z = 1:numel(Trial_4)
        if Trial_4(z).correct == 1 
            datastruct(z,3) = [Trial_4(z).target]';          %choice is in the same direction as tg direction when correct
        elseif Trial_4(z).correct == 0 && Trial_4(z).target == 1
            datastruct(z,3) = 2;
        elseif Trial_4(z).correct == 0 && Trial_4(z).target == 2
            datastruct(z,3) = 1;
        end 
    end 
    datastruct(:,4) = [Trial_4(1:end).correct]';
    datastruct(:,5) = [Trial_4(1:end).RT]';
    datastruct(:,6) = zeros(z,1) + 1;
    if LSC == 1         %if left side inactivation
        datastruct(:,7) = zeros(z,1) + 1;
    elseif RSC == 1     %if right side inactivation
        datastruct(:,7) = zeros(z,1) + 2;
    end 
    if pre == 1           %if pre data
        datastruct(:,8) = zeros(z,1) + 1;
    elseif post == 1       %if post data
        datastruct(:,8) = zeros(z,1) + 2;
    elseif rec == 1
        datastruct(:,8) = zeros(z,1) + 3;
    end
    datastruct(:,9) = zeros(z,1) + 1;


%Save variables for modeling
if files_inputted == 1
    if input('Save variables for modeling? y/n', 's') == 'y'
        savingFilename = input('Name the saved variable with initial of monkey, task, date, and "_VarModel": ', 's'); % Name of file
        save([savingFilename], 'datastruct'); % Save the file
    else
    end
end
stop = 1;
%% Save all the data that will be needed 

if DT == 1
    task_v = 'DT';
else 
    task_v = 'RT';
end 
GPData = struct('Filenames', filename_array, 'Task_version', [task_v], 'Side_of_Inactivation', [LSC,RSC], 'RFx', rfx, 'RFy', rfy, 'Trial_data_array', Trial_data_array, 'RT_struct_array', RTStruct_array,...
    'YAndParams_array', YAndParams_array, 'y_avg_med_RT_array', y_avg_med_RT_array);

% Save the GPData
if files_inputted == 3          %If all the pre, post, and rec files have been analyzed, then save them all into one GPdata structures to later process all experiments
    if input('Save all the info in GPData? y/n', 's') == 'y'
        savingFilename = input('Name the saved variable with initial of monkey, task, date, and "_GPData": ', 's'); % Name of file
        save([savingFilename], 'GPData'); % Save the file
    else
    end
else
end