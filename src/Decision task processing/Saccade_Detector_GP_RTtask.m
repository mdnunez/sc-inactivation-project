function [Trial_3,allh,allv,saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold, accelerationThreshold] = Saccade_Detector_GP_RTtask(Trial, Trial_3, RT, DT)


%Parameters to smooth out eye signals
sigma_saccadeSmoothing = 0.001;
sigmaValues = -sigma_saccadeSmoothing*3:.001:sigma_saccadeSmoothing*3;
kernel = normpdf(sigmaValues, 0, sigma_saccadeSmoothing) * .001;

%Saccade Detector window (Time from GPOn code in RT tasks)
saccadeTimeFrame_0 = 0;
saccadeTimeFrame_1 = 2500;

matrix_columns_size = ((saccadeTimeFrame_1 - saccadeTimeFrame_0)*30)+1;

velocityThreshold = 20;
accelerationThreshold = 5000;

% Shorthands:
% vel = velocity
% acc = acceleration
% smo = smoothened
% uns = unsmoothened
% sac = saccade
% ard = around
% ver = vertical
% hor = horizontal

% ========================================================= %
% === PART 1/2: Unpack the NEV and/or NSx files === %
% ========================================================= %

% -- Checking for pauses in the recording session -- %
% Some NSx files can contain more then one recording session. This is due
% to the recording being paused in between sessions on central, rather than stopped and
% saved as an individual file. The two different cases return different
% data structures. Paused sessions are stored in cell arrays, Stopped
% sessions are recorded in an matrix where the rows represent the
% channels recorded and the colums the number of data points. The code
% below checks for the presence of cells and chooses the appropriate method
% for accessing the data.


ranges = [];
clear x;
index = 1;

if iscell(Trial.Data)
    
    for x = 1:size(Trial.Data{1,1},1)
        ranges = [ranges; range(Trial.Data{1,1}(x,:))]; 
    end 
    hv = find(ranges > 25000);
    
    for i = 1: length(Trial.Data)
        
        if i == 1                                                           %if the file has been paused and there's more than 1 cell for Trial.data
            index_1 = length(Trial.Data{1,i});
        else
            index_1 = length(Trial.Data{1,i}) + Trial.MetaTags.Timestamp(1,i) - 1;
            index = Trial.MetaTags.Timestamp(1,i);
        end
        allh(1,index:index_1) = Trial.Data{1,i}(hv(1),:);
        allv(1,index:index_1)= Trial.Data{1,i}(hv(2),:);
    end
else
    for x = 1:size(Trial.Data,1)
        ranges = [ranges; range(Trial.Data(x,:))]; 
    end 
    hv = find(ranges > 25000);   %makes sure this is analog signal channel we're getting, it ranges to high bits, greater than 25000 
    allh = Trial.Data(hv(1),:);
    allv = Trial.Data(hv(2),:);
end

clear x i;
NumTrials = numel(Trial_3);
%defining matrix size
all_vertical_positions = zeros(NumTrials,matrix_columns_size);
all_horizontal_positions = zeros(NumTrials, matrix_columns_size);
alltimes_index = zeros(NumTrials, matrix_columns_size);

for n = 1: NumTrials         %This is excluding the last trial; For all the trials, get the horiztonal/vertical eye signals within the saccade window defined above after GPon; last trials have invalid data when the recording stops mid-trial
    
    codes = Trial_3(n).eCodes;
    col1 = find(codes == 5000);     %structure index for the 5000 codes; GPOn for RT task

    all_vertical_positions(n,1:matrix_columns_size) = allv(1,Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_3(n).timeindexes(col1)+saccadeTimeFrame_1*30);
    all_horizontal_positions(n,1:matrix_columns_size) = allh(1,Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_3(n).timeindexes(col1)+saccadeTimeFrame_1*30);
    alltimes_index(n,1:matrix_columns_size) = Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_3(n).timeindexes(col1)+ saccadeTimeFrame_1*30;
    
    
end
%print out the ranges of the bits of the eye signal to double check that
%they are eye signals; ranges are usually > 25000
range_allh = range(allh)
range_allv = range(allv)

%Conversion factor of volts to degrees for the eye signal
all_vertical_positions = double((all_vertical_positions+2800).*(5/3400));
all_horizontal_positions = double((all_horizontal_positions+2800).*(5/3400));


%% Putting data into ms time resolution

[row0,col0] = size(all_vertical_positions);

col_sized_factor = floor(col0/30);

if col0 == 30*col_sized_factor
    col_sized = 30*col_sized_factor;
else
    col_sized = 30*col_sized_factor + 1;
end

for b = 1:row0
    q = 1;
    for i = 1: 30: col_sized
        all_vertical_positions_sized(b,q) = double(all_vertical_positions(b,i));
        all_horizontal_positions_sized(b,q) = double(all_horizontal_positions(b,i));
        alltimes_index_sized(b,q) = alltimes_index(b,i);
        
        q = q + 1;
    end
end

%Storing these values into the struc
for j = 1:NumTrials
    Trial_3(j).all_vertical_positions_sized = all_vertical_positions_sized(j,:);
    Trial_3(j).all_horizontal_positions_sized = all_horizontal_positions_sized(j,:);
    Trial_3(j).alltimes_index_sized = alltimes_index_sized(j,:);
end

%%

% % % ================================================================= %
% % % === PART 2/2: Calculate the saccade from the NEV and NSx file === %
% % % ================================================================= %

% The number of ms per trial in the eye position matrix
NumTrial = size(all_vertical_positions_sized,1);
ms_positions = size(all_vertical_positions_sized,2);
% For loop that goes through each trial to calculate the acceleration
stop = 1;

for i = 1:NumTrial
    
    i_readout = i   %just a readout of where the loop is at for detecting a saccade within the trial
    
    vertical_smoothened = conv(Trial_3(i).all_vertical_positions_sized(:),kernel, 'same');
    horizontal_smoothened = conv(Trial_3(i).all_horizontal_positions_sized(:),kernel, 'same');
    
    % ------------------------------ %
    % --- Calculate the velocity --- %
    % ------------------------------ %
    
    % For loop that goes through each ms of the position arrays except the last one
    for j = 1:(ms_positions-1)
       
        % Vertical velocity
        position1 = vertical_smoothened(j);   % Position of the eye at ms
        position2 = vertical_smoothened(j+1); % Position of the eye at ms+1
        current_velocity_vertical = abs(position2 - position1)*1000;
        Trial_3(i).velocity_vertical(j) = current_velocity_vertical; % Store the velocity
        
        % Horizontal velocity
        position1 = horizontal_smoothened(j);   % Position of the eye at ms
        position2 = horizontal_smoothened(j+1); % Position of the eye at ms+1
        current_velocity_horizontal = abs(position2 - position1)*1000;
        Trial_3(i).velocity_horizontal(j) = current_velocity_horizontal;  % Store the velocity
        
        % Diagonal velocity
        velocity_diagonals(j)= sqrt((current_velocity_vertical^2) + (current_velocity_horizontal^2));
        
    end % End of for loop that goes through each ms of the position arrays (j)
    
    %Smoothing
    Trial_3(i).velocity_diagonal = conv(velocity_diagonals(:),kernel, 'same');
    
    % ---------------------------------- %
    % --- Calculate the acceleration --- %
    % ---------------------------------- %
    
    % The number of ms per trial in the velocity matrix
    ms_velocity = numel(Trial_3(i).velocity_diagonal);
    
    % For loop that goes through each ms of the velocity arrays except the
    % last window size + 1 to calculate acceleration
    for j = 1:(ms_velocity - 1) % -1 because we have 1 less ms after calculating acceleration from velocity
        
        velocity1 = Trial_3(i).velocity_diagonal(j);        % Velocity at ms
        velocity2 = Trial_3(i).velocity_diagonal(j+1); % Velocity at ms + window
        current_acceleration = abs(velocity2 - velocity1)*1000;
        acceleration(j) = current_acceleration; % Store the acceleration
    end
    
    %Smoothing
    Trial_3(i).accelerations = conv(acceleration(:),kernel, 'same');
    
    %Making positions/velocity/accelerations even
    Trial_3(i).all_vertical_positions_sized(end-1:end) = [];
    Trial_3(i).all_horizontal_positions_sized(end-1:end) = [];
    Trial_3(i).alltimes_index_sized(end-1:end) = [];
    Trial_3(i).velocity_vertical(end) = [];
    Trial_3(i).velocity_horizontal(end) = [];
    Trial_3(i).velocity_diagonal(end) = [];

    %Getting the first pass saccade thresholds to detect any big changes in the eye
    %movements after GPOn
    FirstPassThresholds = [400, 300, 200, 100, 50];
    saccade_index = [];
    e = 1;
    stop = 1;
    while isempty(saccade_index)
        CellsAboveFirstPassVThreshold = Trial_3(i).velocity_diagonal > FirstPassThresholds(e);
        saccade_index = find(CellsAboveFirstPassVThreshold(50:end-100) == 1, 1, 'first')+49;    %to exclude out the first 50ms from analysis since theres no way he can make a sac 50ms after GPOn, and there's a lot of noise at the beginning
        if max(Trial_3(i).velocity_diagonal(saccade_index:saccade_index+100)) > 1700         %to exclude out any blinks that were falsely caught, they have big peak velocities so use that to identify blinks
            saccade_index = [];
        end 
        e = e + 1;
        if e > numel(FirstPassThresholds)           %if the last pass threshold still didn't catch the sac, make the saccade_index == nan to exit the loop
            saccade_index = nan;
        end
    end
    %If you can't detect any general saccade, make it nan
    if isnan(saccade_index)
        Trial_3(i).first_saccade_index= nan;
        Trial_3(i).first_saccade_time = nan;
        Trial_3(i).PeakVelocity = nan;
    else 
        %But if you do detect a general saccade, get a more precise
        %position of the saccade initiation and peak velocity using your official
        %velocity/acceleration values
        vel_accel_around_sac = [Trial_3(i).velocity_diagonal(saccade_index - 40: saccade_index)';Trial_3(i).accelerations(saccade_index - 40: saccade_index)'];
         precise_sac_index = (size(vel_accel_around_sac,2) - (find(fliplr(vel_accel_around_sac(1,:)<velocityThreshold) & fliplr(vel_accel_around_sac(2,:)<accelerationThreshold), 1, 'first')-1))+1;
         if isempty(precise_sac_index)
             Trial_3(i).first_saccade_index = nan;
             Trial_3(i).first_saccade_time = nan;
             Trial_3(i).PeakVelocity = nan;
         else
             Trial_3(i).first_saccade_index = saccade_index - 40 + precise_sac_index -1;
             Trial_3(i).first_saccade_time = (Trial_3(i).alltimes_index_sized(Trial_3(i).first_saccade_index))/30000;
             Trial_3(i).PeakVelocity = max(Trial_3(i).velocity_diagonal(Trial_3(i).first_saccade_index:Trial_3(i).first_saccade_index+75));
         end 
    end     
    
        
    
end % End of for loop that goes through each trial (i)


done =1;


end
