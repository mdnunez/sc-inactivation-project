function [first_saccade_time_spec, peak_vel_spec, Trial_3] = Saccade_Detector_GP_RTtask_spec1A(Trial, Trial_3, saccadeTimeFrame_0_spec,saccadeTimeFrame_1_spec, velocityThreshold_spec, accelerationThreshold_spec, w, allh, allv, DT, RT, FirstPassThresholds_spec)


%Parameters to smooth out eye signals
sigma_saccadeSmoothing = 0.001;
sigmaValues = -sigma_saccadeSmoothing*3:.001:sigma_saccadeSmoothing*3;
kernel = normpdf(sigmaValues, 0, sigma_saccadeSmoothing) * .001;

matrix_columns_size = ((saccadeTimeFrame_1_spec - saccadeTimeFrame_0_spec)*30)+1;

first_saccade_time_spec = [];
peak_vel_spec = [];
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

%%

all_vertical_positions = zeros(1,matrix_columns_size);
all_horizontal_positions = zeros(1, matrix_columns_size);
alltimes_index = zeros(1, matrix_columns_size);

n = w         %This is excluding the last trial
    
    codes = Trial_3(n).eCodes;
    col1 = find(codes == 5000);     %structure index for the 5000 codes; GPOn
    
    
    all_vertical_positions(1,1:matrix_columns_size) = allv(1,Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0_spec*30: Trial_3(n).timeindexes(col1)+saccadeTimeFrame_1_spec*30);
    all_horizontal_positions(1,1:matrix_columns_size) = allh(1,Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0_spec*30: Trial_3(n).timeindexes(col1)+saccadeTimeFrame_1_spec*30);
    alltimes_index(1,1:matrix_columns_size) = Trial_3(n).timeindexes(col1) + saccadeTimeFrame_0_spec*30: Trial_3(n).timeindexes(col1)+ saccadeTimeFrame_1_spec*30;
    

%Conversion factor of volts to degrees
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

    all_vertical_positions_sized = all_vertical_positions_sized(:);
    all_horizontal_positions_sized = all_horizontal_positions_sized(:);
    alltimes_index_sized = alltimes_index_sized(:);


%     Trial_3(w).all_vertical_positions_sized = all_vertical_positions_sized(w,:);
%     Trial_3(w).all_horizontal_positions_sized = all_horizontal_positions_sized(w,:);
%     Trial_3(w).alltimes_index_sized = alltimes_index_sized(w,:);
%%

% % % ================================================================= %
% % % === PART 2/2: Calculate the saccade from the NEV and NSx file === %
% % % ================================================================= %

% The number of ms per trial in the eye position matrix

ms_positions = size(all_vertical_positions_sized,1);
% For loop that goes through each trial to calculate the acceleration

        
    vertical_smoothened = conv(all_vertical_positions_sized(:),kernel, 'same');
    horizontal_smoothened = conv(all_horizontal_positions_sized(:),kernel, 'same');
   
    
    % ------------------------------ %
    % --- Calculate the velocity --- %
    % ------------------------------ %
    
    % For loop that goes through each ms of the position arrays except the last one
    for e = 1:(ms_positions-1)
        
        %  --- Smoothened ---  %
        
        
        % Vertical velocity
        position1 = vertical_smoothened(e);   % Position of the eye at ms
        position2 = vertical_smoothened(e+1); % Position of the eye at ms+1
        current_velocity_vertical = abs(position2 - position1)*1000;
        velocity_vertical(e) = current_velocity_vertical; % Store the velocity
        
        % Horizontal velocity
        position1 = horizontal_smoothened(e);   % Position of the eye at ms
        position2 = horizontal_smoothened(e+1); % Position of the eye at ms+1
        current_velocity_horizontal = abs(position2 - position1)*1000;
        velocity_horizontal(e) = current_velocity_horizontal;  % Store the velocity
        
        % Diagonal velocity
        velocity_diagonals(e)= sqrt((current_velocity_vertical^2) + (current_velocity_horizontal^2));
        
    end % End of for loop that goes through each ms of the position arrays (e)
    
    %Smoothing
    velocity_diagonal = conv(velocity_diagonals(:),kernel, 'same');
    
    %Getting Peak Velocity for each trial for Velocity Map
%     Trial_3(i).PeakVelocity = max(Trial_3(i).velocity_diagonal(6:end-6));
    
    % ---------------------------------- %
    % --- Calculate the acceleration --- %
    % ---------------------------------- %
    
    % The number of ms per trial in the velocity matrix
    ms_velocity = numel(velocity_diagonal);
    
    % For loop that goes through each ms of the velocity arrays except the last window size + 1
    for e = 1:(ms_velocity - 1) % -1 because we have 1 less ms after calculating acceleration from velocity
        
        velocity1 = velocity_diagonal(e);        % Velocity at ms
        velocity2 = velocity_diagonal(e+1); % Velocity at ms + window
        current_acceleration = abs(velocity2 - velocity1)*1000;
        acceleration(e) = current_acceleration; % Store the acceleration
    end
    
    %Smoothing
    accelerations = conv(acceleration(:),kernel, 'same');
    
    
    %Making positions/velocity/accelerations even
   all_vertical_positions_sized(end-1:end) = [];
    all_horizontal_positions_sized(end-1:end) = [];
   alltimes_index_sized(end-1:end) = [];
    velocity_vertical(end) = [];
   velocity_horizontal(end) = [];
   velocity_diagonal(end) = [];
    
    
    %Getting the first pass thresholds
%     FirstPassThresholds_spec = [400, 300, 200, 100, 50];
    saccade_index = [];
    e = 1;

    while isempty(saccade_index)
        CellsAboveFirstPassVThreshold = velocity_diagonal > FirstPassThresholds_spec(e);
        saccade_index = find(CellsAboveFirstPassVThreshold(50:end-50) == 1, 1, 'first')+49;   
        if max(velocity_diagonal(saccade_index:saccade_index+50)) > 1700         %to exclude out any blinks that were falsely caught, they have big peak velocities so use that to identify blinks
            saccade_index = [];
        end 
        e = e + 1;
        if e > numel(FirstPassThresholds_spec)           %if the last pass threshold still didn't catch the sac, make the saccade_index == nan to exit the loop
            saccade_index = nan;
        end
    end
    %If you can't detect any general saccade, make it nan
    if isnan(saccade_index)
        first_saccade_index= nan;
        first_saccade_time_spec = nan;
        peak_vel_spec = nan;
    else 
        %But if you do detect a general saccade, get a more precise
        %position of the saccade initiation and peak velocity using your official
        %velocity/acceleration values
        vel_accel_around_sac = [velocity_diagonal(saccade_index - 50: saccade_index)';accelerations(saccade_index - 50: saccade_index)'];
         precise_sac_index = (size(vel_accel_around_sac,2) - (find(fliplr(vel_accel_around_sac(1,:)<velocityThreshold_spec) & fliplr(vel_accel_around_sac(2,:)<accelerationThreshold_spec), 1, 'first')-1))+1;
         if isempty(precise_sac_index)
             first_saccade_index = nan;
             first_saccade_time_spec = nan;
             peak_vel_spec = nan;
         else
             first_saccade_index = saccade_index - 30 + precise_sac_index -1;
             first_saccade_time_spec = (alltimes_index_sized(first_saccade_index))/30000;
             peak_vel_spec = max(velocity_diagonal(first_saccade_index:first_saccade_index+75));
         end 
    end     
    
    
        
%          if isempty(precise_sac_index)
%              Trial_3(i).first_saccade_index = nan;
%              Trial_3(i).first_saccade_time = nan;
%              Trial_3(i).PeakVelocity = nan;
%          else
%              Trial_3(i).first_saccade_index = saccade_index - 20 + precise_sac_index -1;
%              Trial_3(i).first_saccade_time = (Trial_3(i).alltimes_index_sized(Trial_3(i).first_saccade_index))/30000;
%              Trial_3(i).PeakVelocity = max(Trial_3(i).velocity_diagonal(Trial_3(i).first_saccade_index:Trial_3(i).first_saccade_index+75));
%          end 
    

done =1;


end

%%
% Changes from Liz from the Saccade_Detector_PV_2013MATLAB version
% - had the first saccade time taken from values from the Trial_3 struct
%- 8/15/19 put in 3rd condition of saccade detection, velocity must be
%above threshold for at least 10msec