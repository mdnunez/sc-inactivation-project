function analog_signal_plotter_STO(Trial_nsx,Trial_1,saccadeTimeFrame_0, w)
clc;
close all;

NumTrials = numel(Trial_1);






% ============================================================================= %
% ======= PART 2/2: Plotting Saccade Position, Velocity , Acceleration ======== %
% ============================================================================= %


w(find(w==0)) = [];
for n = w
    
    first_saccade_time_ms = Trial_1(n).first_saccade_time*1000;
    lineStyle =  '-';
    
    codes = Trial_1(n).eCodes;
    colindex0 = find(codes == 1001, 1, 'first');        %same timestamp as alltimesindex(n,1)
    colindex1 = size(codes, 2);
    colindex2 = find(codes == 3000, 1, 'first');        %struct index for FPOff code
    colindex3 = find(codes == 5050, 1, 'first');
    
    
    %     col5 = Trial_1(n).timeindexes(colindex1)-Trial_1(n).timeindexes(colindex0);
    %     col6 = double(Trial_1(n).timeindexes(colindex2)-Trial_1(n).timeindexes(colindex0))/30;      %FPOFF code ms time for plot
    %     col7 = double(first_saccade_time_ms(n) - (Trial_1(n).timeindexes(colindex0)/30));
    
    col7 = double(first_saccade_time_ms - (Trial_1(n).timeindexes(colindex2)/30) - saccadeTimeFrame_0);
    
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    % Eye Signal Position
    
    plotted_vertical_position = Trial_1(n).all_vertical_positions_sized(1:end);     
    plotted_horizontal_position = Trial_1(n).all_horizontal_positions_sized(1:end);
    alltimes_index_sized_1 = Trial_1(n).alltimes_index_sized(1:end) - Trial_1(n).alltimes_index_sized(1); % Position time index in ms
    
    subplot(3,1,1)
    plot(alltimes_index_sized_1/30, plotted_vertical_position, 'r', alltimes_index_sized_1/30, plotted_horizontal_position, 'b', 'linestyle', lineStyle,'linewidth',2);
    set(gca,'FontSize',12)
stop = 1;
    xlim([0, max(alltimes_index_sized_1)/30]);
    ylim([-40 40]);
    hold on;
    if ~isempty(col7)
        if ~isnan(col7)
            plot([col7,col7], [40, -40])
        end
    end
    %%
    % Eye Velocity
    
    plotted_velocity = Trial_1(n).velocity_diagonal(1:end)';
    
    subplot(3,1,2)
    plot(alltimes_index_sized_1/30, plotted_velocity, 'linestyle', lineStyle,'linewidth',2);
    set(gca,'FontSize',12)
    xlim([0, max(alltimes_index_sized_1)/30]);
    ylim([0 5000]);
    hold on;
    if ~isempty(col7)
        if ~isnan(col7)
            plot([col7,col7], [0, 5000])
        end
    end
    
    % Eye Acceleration
    
    plotted_acceleration = Trial_1(n).accelerations(1:end)';
    
    subplot(3,1,3)
    plot(alltimes_index_sized_1/30, plotted_acceleration, 'linestyle', lineStyle,'linewidth',2);
    set(gca,'FontSize',12)
    xlim([0, max(alltimes_index_sized_1)/30]);
    ylim([0 100000]);
    hold on;
    if ~isempty(col7)
        if ~isnan(col7)
            plot([col7,col7], [0, 100000])
        end
    end
    
    dim = [0.3 0.3 0.3 0.3];
    str = {'Trial' num2str(n), ', PeakVelocity=' num2str(Trial_1(n).PeakVelocity) ', RT=' num2str(Trial_1(n).RT)};
    annotation('textbox', dim, 'String',str,'FitBoxToText','on', 'Linestyle', 'none');
    
    done = 1;
    
    
end % end of the loop that goes through each of the trials
done =1;
end