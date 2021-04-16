%To process the raw data of the simple selection task 
clc;
clear all;
close all;

%----------MODIFY MAIN SETTINGS HERE -------------------------------------
NEV = 1; %if the data was recorded with a NEV file, =1
REX = 0; %if the data was recorded with a REX file, =1
%-------------------------------------------------------------------------

continueAction = 1;                                       
i = 1;                                                   
files_inputted = 1;

disp('%%---------------------------------------------%%');
prompt_1 = 'Input rf parameters x and y for the file. x: ';
a = input(prompt_1);
prompt_2 = 'y: ';
b = input(prompt_2);

while(continueAction)
[saccAwayRF, saccToRF, Trial_1, TrialTargInRF,TrialTargOutRF,PropTrialInRF,PropTrialOutRF, percCin, percCopp,filename] = selectiontask_rawdataprocessor(a,b,files_inputted, NEV, REX);    
    display('Last added file was');
    filename
    prompt_0 = ('Would you like to add a file? (y/n)');     %variable for prompt
    add_file = input(prompt_0, 's');                        %varibale for user input string into command line after prompt is shown
    i;                                                      %i exists
    if a > 0                                                % if a is positive, RFx is positive
        left_right_array(1,1) = {'AwayRF (Left)'};          %away RF is left 
        left_right_array(2,1) = {'ToRF (Right)'};           %toRF is right or else
        Trial_distribution_array(1,1) = {'AwayRF (Left)'};
        Trial_distribution_array(2,1) = {'ToRF (Right)'};
        Acc_array(1,1) = {'AwayRF (Left)'};
        Acc_array(2,1) = {'ToRF (Right)'};
    else                                                    %a is negative, RFx is negative
        left_right_array(1,1) = {'AwayRF (Right)'};         %away RF is right
        left_right_array(2,1) = {'ToRF (Left)'};            %toRF is left
        Trial_distribution_array(1,1) = {'AwayRF (Right)'};
        Trial_distribution_array(2,1) = {'ToRF (Left)'};
        Acc_array(1,1) = {'AwayRF (Right)'};
        Acc_array(2,1) = {'ToRF (Left)'};
    end
    left_right_array(i,2) = {saccAwayRF};
    left_right_array(i+1,2) = {saccToRF};
    Trial_distribution_array(i,2) = {TrialTargOutRF};
    Trial_distribution_array(i+1,2) = {TrialTargInRF};
    Trial_distribution_array(i,3) = {PropTrialOutRF};
    Trial_distribution_array(i+1,3) = {PropTrialInRF};
    Acc_array(i,2) = {percCopp};
    Acc_array(i+1,2) = {percCin};
    Trial_1_array{files_inputted,1} = Trial_1;
    
    i = i + 2;          %index for the next added data plot
    if add_file == 'Y' || add_file == 'y'
        continueAction = 1;
        files_inputted = files_inputted + 1;
    elseif add_file == 'N' || add_file == 'n'
        continueAction = 0;
        continue;
    end
end

%%
%-------------------- Plot Histogram --------------------%
%for all the histograms, pre-injection data gets plotted first (on the
%left) and then post-injection data gets plotted after (on the right) 

%Plotting for proportion made to the left/right 
m = length(left_right_array)/2;     %number of panels
j = 1;  
k = 2;

figure ('units', 'normalized', 'outerposition', [0 0 1 1]); 
ylim([0 0.6]);% Y-Axis Limits

close all
for s = 1 : m                                                                       %vector of panel numbers
    s_1 = subplot(1,m,s);                                                           %creates subplots within a figure and assigning a variable to it; subplots(rows, columns, position)
    bar ([left_right_array{j:k,2}]);
    set(gca, 'Xtick', [1 2],'xticklabels',left_right_array(1:2,1),'Fontsize',10);   %set axis labels 
    dim4 = [0.2, 0.7, 0.9, 0.3];
    str4 = 'Proportion of trials of saccades made away or towards the RF';
%     annotation('textbox', dim4, 'String',str4,'FitBoxToText','on', 'Linestyle', 'none');
    ylim([0 0.6]);% Y-Axis Limits
    text(1,0.9, {'awayRF =' num2str(left_right_array{j,2}) 'toRF=' num2str(left_right_array{k,2})}, 'color', 'r');
    j = j + 2;
    k = k + 2;
    
end

%Plotting for Trial Distribution
figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
h = 1;
p = 2;
ylim([0 1]);% Y-Axis Limits

for s = 1 : m                                                                       %vector of panel numbers
    s_1 = subplot(1,m,s);                                                           %creates subplots within a figure and assigning a variable to it; subplots(rows, columns, position)
   bar ([Trial_distribution_array{h:p,3}], 'g');
    set(gca, 'Xtick', [1 2],'xticklabels',{'TrialTargOutRF', 'TrialTargInRF'},'Fontsize',10);   %set axis labels 
    dim4 = [0.2, 0.7, 0.9, 0.3];
    str4 = {'Proportion of trials with targets set away or towards the RF'};
%     annotation('textbox', dim4, 'String',str4,'FitBoxToText','on', 'Linestyle', 'none');
    ylim([0 1]);% Y-Axis Limits
    text(1,0.9, {'awayRF =' num2str(Trial_distribution_array{h,2}), num2str(Trial_distribution_array{h,3}), 'toRF=' num2str(Trial_distribution_array{p,2}), num2str(Trial_distribution_array{p,3})}, 'color', 'r');     
    h = h + 2;
    p = p + 2;
   
end

%Plotting for Accuracy 
q = 1;  
w = 2;
figure ('units', 'normalized', 'outerposition', [0 0 1 1]); 
ylim([0 1.1]);% Y-Axis Limits

for s = 1 : m                                                                       %vector of panel numbers
    s_1(m) = subplot(1,m,s);                                                           %creates subplots within a figure and assigning a variable to it; subplots(rows, columns, position)
    bar ([Acc_array{q:w,2}], 'y');
    set(gca, 'Xtick', [1 2],'xticklabels',{'OutRF' 'inRF'},'Fontsize',30);   %set axis labels 
    dim4 = [0.2, 0.7, 0.9, 0.3];
    str4 = 'Proportion correct when target is away from RF or towards RF';
%     annotation('textbox', dim4, 'String',str4,'FitBoxToText','on', 'Linestyle', 'none');
    ylim([0 1.1]);% Y-Axis Limits
    text(1,0.9, {'OutRF =' num2str(Acc_array{q,2}) 'InRF=' num2str(Acc_array{w,2})}, 'color', 'k','Fontsize',30);
    q = q + 2;
    w = w + 2;
end


%% Saving data structure

STOData = struct();
STOData.Left_Right_Array = left_right_array;
STOData.Trial_Distribution_Array = Trial_distribution_array;
STOData.Acc_array = Acc_array;
STOData.Trial_1 = Trial_1_array;

if input('Save all the info in STOData? y/n', 's') == 'y'
        savingFilename = input('Name the saved variable with initial of monkey, task, date, and "_STOData": ', 's'); % Name of file
        save(['/' savingFilename], 'STOData'); % Save the file
end 