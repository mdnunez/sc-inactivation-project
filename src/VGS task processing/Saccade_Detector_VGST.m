function [Trial_4,allh,allv,saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold, accelerationThreshold] = Saccade_Detector_VGST(Trial, Trial_3)

Trial_4 = Trial_3;
%Parameters to smooth out eye signals
sigma_saccadeSmoothing = 0.001;
sigmaValues = -sigma_saccadeSmoothing*3:.001:sigma_saccadeSmoothing*3;
kernel = normpdf(sigmaValues, 0, sigma_saccadeSmoothing) * .001;

%Saccade Detector window (Time from FPOFF code)
saccadeTimeFrame_0 = 0;
saccadeTimeFrame_1 = 1100;

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
    hv = find(ranges > 25000);
    allh = Trial.Data(hv(1),:);
    allv = Trial.Data(hv(2),:);
end

clear x i;

%defining matrix size
NumTrials = numel(Trial_4);

all_vertical_positions = zeros(NumTrials,matrix_columns_size);
all_horizontal_positions = zeros(NumTrials, matrix_columns_size);
alltimes_index = zeros(NumTrials, matrix_columns_size);

for n = 1: NumTrials          %This is excluding the last trial
    
    codes = Trial_4(n).eCodes;
    col1 = find(codes == 3000);     %structure index for the 3000 codes
    
    all_vertical_positions(n,1:matrix_columns_size) = allv(1,Trial_4(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_4(n).timeindexes(col1)+saccadeTimeFrame_1*30);
    all_horizontal_positions(n,1:matrix_columns_size) = allh(1,Trial_4(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_4(n).timeindexes(col1)+saccadeTimeFrame_1*30);
    alltimes_index(n,1:matrix_columns_size) = Trial_4(n).timeindexes(col1) + saccadeTimeFrame_0*30: Trial_4(n).timeindexes(col1)+ saccadeTimeFrame_1*30;
    
    %Storing these values in the structure
    Trial_4(n).all_vertical_positions(1:matrix_columns_size) = all_vertical_positions(n,1:matrix_columns_size);
    Trial_4(n).all_horizontal_positions(1:matrix_columns_size) = all_horizontal_positions(n,1:matrix_columns_size);
    Trial_4(n).alltimes_index(1:matrix_columns_size) = alltimes_index(n,1:matrix_columns_size);
    
end
range_allh = range(allh)
range_allv = range(allv)

%Conversion factor of volts to degrees
all_vertical_positions = double((all_vertical_positions+2800).*(5/3400));
all_horizontal_positions = double((all_horizontal_positions+2800).*(5/3400));

for n = 1: NumTrials
    % Converting eye position from A/D units to degree of visual angle
    Trial_4(n).all_vertical_positions(:) = (Trial_4(n).all_vertical_positions(:)+2800).*(5/3400);
    Trial_4(n).all_horizontal_positions(:) = (Trial_4(n).all_horizontal_positions(:)+2800).*(5/3400);
end

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
    Trial_4(j).all_vertical_positions_sized = all_vertical_positions_sized(j,:);
    Trial_4(j).all_horizontal_positions_sized = all_horizontal_positions_sized(j,:);
    Trial_4(j).alltimes_index_sized = alltimes_index_sized(j,:);
end
%%
% % % ================================================================= %
% % % === PART 2/2: Calculate the saccade from the NEV and NSx file === %
% % % ================================================================= %

% The number of ms per trial in the eye position matrix
NumTrial = size(all_vertical_positions_sized,1);
ms_positions = size(all_vertical_positions_sized,2);
% For loop that goes through each trial to calculate the acceleration
for i = 1:NumTrial
    
    i_readout = i
    
    Trial_4(i).vertical_smoothened = conv(Trial_4(i).all_vertical_positions_sized(:),kernel, 'same');
    Trial_4(i).horizontal_smoothened = conv(Trial_4(i).all_horizontal_positions_sized(:),kernel, 'same');
    
    
    % ------------------------------ %
    % --- Calculate the velocity --- %
    % ------------------------------ %
    
    % For loop that goes through each ms of the position arrays except the last one
    for j = 1:(ms_positions-1)                
        
        % Vertical velocity
        position1 = Trial_4(i).vertical_smoothened(j);   % Position of the eye at ms
        position2 = Trial_4(i).vertical_smoothened(j+1); % Position of the eye at ms+1
        current_velocity_vertical = abs(position2 - position1)*1000;
        Trial_4(i).velocity_vertical(j) = current_velocity_vertical; % Store the velocity
        
        % Horizontal velocity
        position1 = Trial_4(i).horizontal_smoothened(j);   % Position of the eye at ms
        position2 = Trial_4(i).horizontal_smoothened(j+1); % Position of the eye at ms+1
        current_velocity_horizontal = abs(position2 - position1)*1000;
        Trial_4(i).velocity_horizontal(j) = current_velocity_horizontal;  % Store the velocity
        
        % Diagonal velocity
        Trial_4(i).velocity_diagonals(j)= sqrt((current_velocity_vertical^2) + (current_velocity_horizontal^2));
        
    end % End of for loop that goes through each ms of the position arrays (j)
    
    %Smoothing
    Trial_4(i).velocity_diagonal = conv(Trial_4(i).velocity_diagonals(:),kernel, 'same');
    
    % ---------------------------------- %
    % --- Calculate the acceleration --- %
    % ---------------------------------- %
    
    % The number of ms per trial in the velocity matrix
    ms_velocity = numel(Trial_4(i).velocity_diagonal);
    
    % For loop that goes through each ms of the velocity arrays except the last window size + 1
    for j = 1:(ms_velocity - 1) % -1 because we have 1 less ms after calculating acceleration from velocity
        velocity1 = Trial_4(i).velocity_diagonal(j);        % Velocity at ms
        velocity2 = Trial_4(i).velocity_diagonal(j+1); % Velocity at ms + window
        current_acceleration = abs(velocity2 - velocity1)*1000;
        Trial_4(i).acceleration(j) = current_acceleration; % Store the acceleration
    end
    
    %Smoothing
    Trial_4(i).accelerations = conv(Trial_4(i).acceleration(:),kernel, 'same');
    
    %Making positions/velocity/accelerations even
    Trial_4(i).all_vertical_positions_sized(end-1:end) = [];
    Trial_4(i).all_horizontal_positions_sized(end-1:end) = [];
    Trial_4(i).alltimes_index_sized(end-1:end) = [];
    Trial_4(i).vertical_smoothened(end-1:end) = [];
    Trial_4(i).horizontal_smoothened(end-1:end) = [];
    
    Trial_4(i).velocity_vertical(end) = [];
    Trial_4(i).velocity_horizontal(end) = [];
    Trial_4(i).velocity_diagonals(end) = [];
    Trial_4(i).velocity_diagonal(end) = [];

    cells_above_accel_threshold = Trial_4(i).accelerations > accelerationThreshold;
    cells_above_vel_threshold =  Trial_4(i).velocity_diagonal > velocityThreshold;
    vel_accel_cells_logical = cells_above_accel_threshold + cells_above_vel_threshold;
    vel_accel_cells_logical(1:5) = [];      %getting rid of the first 5 points bc they're distorted by smoothing
    logical_index = vel_accel_cells_logical > 1;
    first_saccade_index = find(logical_index == 1, 1, 'first') + 5;
    
    first_saccade_index_list = find(logical_index == 1);
    first_saccade_index_list = first_saccade_index_list +5;
        
    while isempty(Trial_4(i).first_saccade_time)
        for z = 1:numel(first_saccade_index_list)       %for all the first saccade candidates that passed the vel, accel threshold
                z_readout = z
                if numel(Trial_4(i).velocity_diagonal(first_saccade_index_list(z):end)) > 10 %if the next 20 velocity indexes, from the first saccade, are a number
                    cells_above_vel_threshold_time = [];
                    for x = 1:10                               %for the next 20 1ms data points, make sure the velocity is above the threshold
                        cells_above_vel_threshold_time = [cells_above_vel_threshold_time, Trial_4(i).velocity_diagonal(first_saccade_index_list(z)+x) > velocityThreshold];
                    end
                    
                    if sum(cells_above_vel_threshold_time) >= numel(cells_above_vel_threshold_time)
                        Trial_4(i).first_saccade_time = Trial_4(i).alltimes_index_sized(first_saccade_index_list(z))/30000;
                        Trial_4(i).first_saccade_index = first_saccade_index_list(z);
                        if numel(Trial_4(i).velocity_diagonal) > first_saccade_index_list(z)+100
                            Trial_4(i).PeakVelocity = max(Trial_4(i).velocity_diagonal(first_saccade_index_list(z):first_saccade_index_list(z)+100));
                        else
                            Trial_4(i).PeakVelocity = max(Trial_4(i).velocity_diagonal(first_saccade_index_list(z):first_saccade_index_list(z)+25));
                        end 
                        break
                    else
                    end
                    
                else
                end
                

            if z == numel(first_saccade_index_list)
                Trial_4(i).first_saccade_time = nan;
                Trial_4(i).first_saccade_index = nan;
                Trial_4(i).PeakVelocity = nan;
            end 
            
        end
    end
  
end % End of for loop that goes through each trial (i)

end
