function [Trial_2, remove_3] = SacDetectManualCheck_STO(Trial_nsx, Trial_1, saccadeTimeFrame_0,saccadeTimeFrame_1, velocityThreshold,accelerationThreshold,allh, allv)
Trial_2 = Trial_1;

N = numel(Trial_2);

remove_3 = [];
FirstPassThresholds = [300, 200, 100];
nPlotsBunch = 50;
t = 1:nPlotsBunch:numel(Trial_2);
input_plot = 'y';

input_plot1 = 'z';
while input_plot1 ~= 'y' && input_plot1 ~= 'n'
    input_plot1 = input(['Plot all trials to check? y/n: '], 's');
    if input_plot1 == 'y'
        for i = 1:ceil(numel(Trial_2)/nPlotsBunch)
            if input_plot == 'y'
                if i == ceil(numel(Trial_2)/nPlotsBunch)
                    w = t(i):numel(Trial_2);
                else
                    w = t(i):t(i+1)-1;                %vector of the trial number needed to be plotted, 50 at a time
                end
                analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);
                input_plot = input(['Plot another ' num2str(nPlotsBunch) ' more? y/n: '],'s');
            end
        end
    end
end

clc;
input_doyou = 'z';
while input_doyou ~= 'y' && input_doyou ~= 'n'
    input_doyou = input('Do you want look at or fix some trials? y/n: ', 's');
end
if input_doyou == 'y'
    input_trials = input('Which trials do you want to look at or fix? Define a vector of trials, or write 0 for none: \n');
    lookat = input_trials;
    
    for c = 1:numel(lookat)       %for each out of range trial,
        tryagain = 1;
        while tryagain == 1
            w = lookat(c);
            analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);    %graph the specific trial out
            discard_input = 'z';
            while discard_input ~= 'a' && discard_input ~= 'b' && discard_input ~= 'c' && discard_input ~= 'd' && discard_input ~= 'e' && discard_input ~= 'f' && discard_input ~= 'g' && discard_input ~= 'h'
                clc;
                prompt_2a = ['For this trial, type ...'...
                    '\n "a" to discard trial from analysis'...
                    '\n "b" to alter saccade detection window'...
                    '\n "c" to alter threshold'...
                    '\n "d" to continue'...
                    '\n "e" to escape and return'...
                    '\n "f" to alter BOTH window and thresholds'...
                    '\n "g" to input the RT manually'...
                    '\n "h" to adjust First Pass Threshold: '];
                discard_input = input(prompt_2a, 's');
            end
            if discard_input == 'a'
                remove_3 = [remove_3, w];
                tryagain = 0;
            elseif discard_input == 'b'
                saccadeTimeFrame_0_spec = saccadeTimeFrame_0;
                saccadeTimeFrame_1_spec = saccadeTimeFrame_1;
                while saccadeTimeFrame_0_spec > 2500 || saccadeTimeFrame_1_spec > 2500
                    saccadeTimeFrame_0_spec = str2num(input('Input saccade detection 1st bound (in ms):', 's'));
                    saccadeTimeFrame_1_spec = str2num(input('Input saccade detection 2nd bound (in ms):', 's'));
                end
                velocityThreshold_spec = velocityThreshold;
                accelerationThreshold_spec = accelerationThreshold;
                FirstPassThresholds_spec = FirstPassThresholds;
                [first_saccade_time_spec, peak_vel_spec, Trial_2] = Saccade_Detector_STO_spec(Trial_nsx, Trial_2, saccadeTimeFrame_0_spec,saccadeTimeFrame_1_spec, velocityThreshold_spec, accelerationThreshold_spec, w, allh, allv, FirstPassThresholds_spec);
                for n = w
                    Trial_2(n).first_saccade_time = first_saccade_time_spec;
                    Trial_2(n).RT =  first_saccade_time_spec -  Trial_2(n).FPOff;
                    Trial_2(n).PeakVelocity = peak_vel_spec;
                end
               analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);    %graph the specific trial out
                prompt_2f_input = 'z';
                while prompt_2f_input ~= 'd' && prompt_2f_input ~= 'a' && prompt_2f_input ~= 't'
                    clc
                    prompt_2f = ['Keep or discard adjusted trial? "d" for keep, "a" for discard, or "t" for try again: '];
                    prompt_2f_input = input(prompt_2f, 's');
                end
                if prompt_2f_input == 'a'
                    remove_3 = [remove_3, w];
                    tryagain = 0;
                elseif prompt_2f_input == 't'
                    tryagain = 1;
                elseif prompt_2f_input == 'd'
                    tryagain = 0;
                end
            elseif discard_input == 'c'
                saccadeTimeFrame_0_spec = saccadeTimeFrame_0;
                saccadeTimeFrame_1_spec = saccadeTimeFrame_1;
                velocityThreshold_spec = str2num(input('Input saccade detection velocity threshold (degrees per second):', 's'));
                accelerationThreshold_spec = str2num(input('Input saccade detection acceleration threshold (degrees per second^2):', 's'));
                FirstPassThresholds_spec = FirstPassThresholds;
                [first_saccade_time_spec, peak_vel_spec, Trial_3] = Saccade_Detector_STO_spec(Trial_nsx, Trial_2, saccadeTimeFrame_0_spec,saccadeTimeFrame_1_spec, velocityThreshold_spec, accelerationThreshold_spec, w, allh, allv, FirstPassThresholds_spec);
                for n = w
                    Trial_2(n).first_saccade_time = first_saccade_time_spec;
                    Trial_2(n).RT =  first_saccade_time_spec -  Trial_2(n).FPOff;
                    Trial_2(n).PeakVelocity = peak_vel_spec;
                end
                analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);    %graph the specific trial out
                prompt_2f_input = 'z';
                while prompt_2f_input ~= 'd' && prompt_2f_input ~= 'a' && prompt_2f_input ~= 't'
                    clc
                    prompt_2f = ['Keep or discard adjusted trial? "d" for keep, "a" for discard, "t" for try again: '];
                    prompt_2f_input = input(prompt_2f, 's');
                end
                if prompt_2f_input  == 'a'
                    remove_3 = [remove_3, w];
                    tryagain = 0;
                elseif prompt_2f_input == 'd'
                    tryagain = 0;
                elseif prompt_2f_input == 't'
                    tryagain = 1;
                end
            elseif discard_input == 'e'
                dbstop in 'SacDetectManualCheck_STO' at '121'
                tryagain = 0;
            elseif discard_input == 'f'
                saccadeTimeFrame_0_spec = saccadeTimeFrame_0;
                saccadeTimeFrame_1_spec = saccadeTimeFrame_1;
                while saccadeTimeFrame_0_spec > 2500 || saccadeTimeFrame_1_spec > 2500
                    saccadeTimeFrame_0_spec = str2num(input('Input saccade detection 1st bound (in ms):', 's'));
                    saccadeTimeFrame_1_spec = str2num(input('Input saccade detection 2nd bound (in ms):', 's'));
                    velocityThreshold_spec = str2num(input('Input saccade detection velocity threshold (degrees per second):', 's'));
                    accelerationThreshold_spec = str2num(input('Input saccade detection acceleration threshold (degrees per second^2):', 's'));
                end
                FirstPassThresholds_spec = FirstPassThresholds;
                [first_saccade_time_spec, peak_vel_spec, Trial_3] = Saccade_Detector_STO_spec(Trial_nsx, Trial_2, saccadeTimeFrame_0_spec,saccadeTimeFrame_1_spec, velocityThreshold_spec, accelerationThreshold_spec, w, allh, allv, FirstPassThresholds_spec);
                for n = w
                    Trial_2(n).first_saccade_time = first_saccade_time_spec;
                    Trial_2(n).RT =  first_saccade_time_spec -  Trial_2(n).FPOff;
                    Trial_2(n).PeakVelocity = peak_vel_spec;
                end
                analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);    %graph the specific trial out
                prompt_2f_input = 'z';
                while prompt_2f_input ~= 'd' && prompt_2f_input ~= 'a' && prompt_2f_input ~= 't'
                    clc
                    prompt_2f = ['Keep or discard adjusted trial? "d" for keep, "a" for discard, "t" for tryagain: '];
                    prompt_2f_input = input(prompt_2f, 's');
                end
                if prompt_2f_input == 'a'
                    remove_3 = [remove_3, w];
                    tryagain = 0;
                elseif prompt_2f_input == 'd'
                    tryagain = 0;
                elseif prompt_2f_input == 't'
                    tryagain = 1;
                end
            elseif discard_input == 'g'
                Trial_2(n).RT = input(['Manually input the RT here: '])/1000;
                Trial_2(n).PeakVelocity = input(['Manually input the PV here: ']);
                prompt_2f_input = 'z';
                while prompt_2f_input ~= 'd' && prompt_2f_input ~= 'a' && prompt_2f_input ~= 't'
                    clc
                    prompt_2f = ['Keep or discard adjusted trial? "d" for keep, "a" for discard, "t" for tryagain: '];
                    prompt_2f_input = input(prompt_2f, 's');
                end
                if prompt_2f_input == 'a'
                    remove_3 = [remove_3, w];
                    tryagain = 0;
                elseif prompt_2f_input == 'd'
                    tryagain = 0;
                elseif prompt_2f_input == 't'
                    tryagain = 1;
                end
            elseif discard_input == 'h'
                saccadeTimeFrame_0_spec = saccadeTimeFrame_0;
                saccadeTimeFrame_1_spec = saccadeTimeFrame_1;
                velocityThreshold_spec = velocityThreshold;
                accelerationThreshold_spec = accelerationThreshold;
                FirstPassThresholds_spec = input(['Enter First Pass Threshold in vector form: ']);
                [first_saccade_time_spec, peak_vel_spec, Trial_3] = Saccade_Detector_STO_spec(Trial_nsx, Trial_2, saccadeTimeFrame_0_spec,saccadeTimeFrame_1_spec, velocityThreshold_spec, accelerationThreshold_spec, w, allh, allv, FirstPassThresholds_spec);
                for n = w
                    Trial_2(n).first_saccade_time = first_saccade_time_spec;
                    Trial_2(n).RT =  first_saccade_time_spec -  Trial_2(n).FPOff;
                    Trial_2(n).PeakVelocity = peak_vel_spec;
                end
                analog_signal_plotter_STO(Trial_nsx,Trial_2,saccadeTimeFrame_0, w);    %graph the specific trial out
                prompt_2f_input = 'z';
                while prompt_2f_input ~= 'd' && prompt_2f_input ~= 'a' && prompt_2f_input ~= 't'
                    clc
                    prompt_2f = ['Keep or discard adjusted trial? "d" for keep, "a" for discard, "t" for tryagain: '];
                    prompt_2f_input = input(prompt_2f, 's');
                end
                if prompt_2f_input == 'a'
                    remove_3 = [remove_3, w];
                    tryagain = 0;
                elseif prompt_2f_input == 'd'
                    tryagain = 0;
                elseif prompt_2f_input == 't'
                    tryagain = 1;
                end
            elseif discard_input == 'd'
                tryagain = 0;
            end
            
            
            
            
            
        end
    end
    
end


end