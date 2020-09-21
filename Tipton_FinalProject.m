%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Final Project
% Filename: Tipton_FinalProject.m 
% Author: Natalie Tipton
% Class: EGR 635/534
% Date: 12/5/19
% Instructor: Dr. Rhodes
% Description: This algorithm detects the presence of Atrial Fibrillation
%   or normal sinus rhythm in an ECG signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

path = cd;

files = dir('Training\Normal\*.mat');   %open all .csv files in current folder w/ CS0 in name
%files = dir('Test\Normal\*.mat');   %open all .csv files in current folder w/ CS0 in name
%files = dir('Training\AF\*.mat');   %open all .csv files in current folder w/ CS0 in name
%files = dir('Test\AF\*.mat');   %open all .csv files in current folder w/ CS0 in name
%files = dir('Training\Other\*.mat');   %open all .csv files in current folder w/ CS0 in name
%files = dir('Test\Other\*.mat');   %open all .csv files in current folder w/ CS0 in name

[num_files,z] = size(files);    %determine number of files read

for x = 1:num_files         % run code for each data file
    if x > 1                % after first file, delete all variable values except cumulative ones
        clearvars -except ent_atrial_activity ent_qrs_removed ent_raw_data ent_qq_interval ent_wavelet...
            x path num_files files classification;
    end
    
    cd('Training\Normal');
    %cd('Test\Normal');
    %cd('Training\AF');
    %cd('Test\AF');
    %cd('Training\Other');
    %cd('Test\Other');
    
    load(files(x).name);    %read in all images in directory
    cd(path); 

    val = (val - min(val)) ./ (max(val) - min(val));
    
    entropy_data = shannon_ent(val, 0.001);
    
    fs = 300;                   % sampling frequency
    len = length(val);          % length of data
    t = 0:1/fs:len/fs - 1/fs;   % time vector
    
    % plot original data
    figure
    plot(t,val);
    title('Original Signal'); xlabel('Time (s)'); ylabel('Amplitude');
    
    %wavelet entropy calculations on original data
    [aa, coeffs, entropy_wavelet] = wavelets(val, t, fs);
    
    % revectorize data into overlapping vecotrs of size 300 with 50% overlap
    aa_windowed = buffer(1:numel(aa),300,150);
    aa_overlap = aa(aa_windowed(:,all(aa_windowed)))';
    [row,col] = size(aa_overlap);   % find size of vecotrized atrial activity
    
    % find entropy of each second of data
    for i = 1:row
        entropy_moving_aa(i) = shannon_ent(aa_overlap(i,:), 0.001);
    end
    
    % average moving entropy values
    entropy_aa = mean(entropy_moving_aa);
    
    % plot atrial activity from wavelet decomposition
    figure
    plot(t,aa);
    title('Atrial Activity'); xlabel('Time (s)'); ylabel('Amplitude');
    
    % bandpass filter the original data for time domain calculations
    ecg = removeNoise(5, 15, fs, t, val);
    
    % find peaks in data
    [pks, locs] = findpeaks(ecg);
    
    % plot peaks
    figure
    findpeaks(ecg);
    title('Peaks in ECG'); xlabel('Time (s)'); ylabel('Amplitude');
    
    % pick out the QRS peaks and their location index in data
    qrs_pks = pks(pks > 0.6);
    qrs_locs = locs(pks > 0.6);
    
    % plot QRS peaks over orignal data
    figure
    plot(t, val)
    hold on
    plot(t(qrs_locs), val(qrs_locs), 'o');
    title('QRS Peaks'); xlabel('Time (s)'); ylabel('Amplitude');
    legend('ECG', 'QRS Peak');
    
    % find outliers due to noise that are much higher than avg QRS peak
    outlier = find(qrs_pks > 1.5 * median(qrs_pks));
    
    % % plot the outlier peaks
    % figure
    % plot(t,val)
    % hold on
    % plot(t(qrs_locs(outlier)), val(qrs_locs(outlier)), 'o');
    % title('Outliers Detected'); xlabel('Time (s)'); ylabel('Amplitude');
    % legend('ECG', 'Outlier');
    
    % remove the outlier data up to the preceding and following QRS peak
    for i = length(outlier) : -1:  1
        ecg(qrs_locs(outlier(i)-1) : qrs_locs(outlier(i)+1)) = mean(ecg);
    end
    
    % % plot data with outliers removed
    % figure
    % plot(t,ecg)
    % title('ECG With Outliers Removed'); xlabel('Time (s)'); ylabel('Amplitude');
    
    % find Q and S points
    [q_points, s_points, t_2, v2] = QS_pointDetect(ecg, t);
    
    qq_int = diff(q_points(1:end));         % find Q to Q intervals
    entropy_qq = shannon_ent(qq_int, 0.001);   % calculate entropy of Q to Q intervals
    
    % Plot Q and S points
    figure
    plot(t_2,v2)
    hold on
    plot(t_2(q_points(:)), v2(q_points(:)), 'o');
    plot(t_2(s_points(:)), v2(s_points(:)), 'o');
    title('Pan-Tompkins with Q and S detected');
    legend('ECG', 'q', 's'); xlabel('Time (s)'); ylabel('Amplitude');
    
    % smooth data for slope analysis
    smoothed = smoothdata(ecg);
    smoothed = smoothdata(smoothed);
    
    % plot smoothed data
    figure
    plot(t,smoothed);
    title('Smoothed data'); xlabel('Time (s)'); ylabel('Amplitude');
    
    % refine location of the S point by finding where slope goes positive
    for i = 1:length(q_points) - 1
        s_points(i) = s_points(i) + find(diff(val(s_points(i):q_points(i+1))) > 0,1,'first');
    end
    
    % based on S-point location, find the peak of the t-wave
    for i = 1:length(q_points) - 1
        done = 0;
        start = s_points(i);    % start looking for slope change at s point
        while done == 0
            % find where slope goes negative (peak of T-wave)
            top_peak = start + find(diff(val(start:q_points(i+1))) < 0,1,'first');
            
            % check if next 5 differences are also negative
            if (all(diff(val(top_peak:top_peak + 10)) < 0))
                t_top(i) = top_peak;    % if so, it is really the peak
                done = 1;
            else
                start = top_peak;       % if not, it was from noise, recheck starting from that point
            end
        end
    end
    
    % based on peak of the t-wave, find the end of the t-wave
    for i = 1:length(q_points) - 1
        t_end(i) = t_top(i) + find(diff(val(t_top(i):q_points(i+1))) > 0,1,'first');
    end
    
    % subtract an average QRS complex from each pulse to remove ventricular activity
    [qrs, qrs_avg, qrs_removed] = QRSaveraging(val, q_points, t_end ,fs);
    
    % revectorize data into overlapping vecotrs of one second length
    qrs_rem_windowed = buffer(1:numel(qrs_removed),300,150);
    qrs_rem_overlap = qrs_removed(qrs_rem_windowed(:,all(qrs_rem_windowed)))';
    [row,col] = size(qrs_rem_overlap);
    
    % find entropy of each second of data with QRS-T removed
    for i = 1:row
        entropy_moving_qrs(i) = shannon_ent(qrs_rem_overlap(i,:), 0.001);
    end
    
    % find average entropy of all thirty seconds
    entropy_qrs_rem = mean(entropy_moving_qrs);
    
    % plot Q, S, and T points
    figure
    plot(t, val)
    hold on
    plot(t(q_points(:)), val(q_points(:)), 'o');
    plot(t(s_points(:)), val(s_points(:)), 'o');
    plot(t(t_end(:)), val(t_end(:)), 'o');
    title('Q, S, and T points detected on data');
    legend('ECG', 'q', 's', 't');
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    
    % Save all entropy values to an array
    ent_atrial_activity(x) = entropy_aa;
    ent_qrs_removed(x) = entropy_qrs_rem;
    ent_raw_data(x) = entropy_data;
    ent_qq_interval(x) = entropy_qq;
    ent_wavelet(x) = entropy_wavelet;
    
    % classify rhythm based on Q to Q interval entropy
    if(entropy_qq < 4)
        classification(x) = 0;      % normal sinus rhythm
    else
        classification(x) = 1;      % atrial fibrillation
    end
    
end

% writes data to an excel file
% myfile = 'Project_Results.xlsx' ;
% xlswrite(myfile,classification,'A42:J42')
% xlswrite(myfile,ent_qrs_removed,'A38:J38')
% xlswrite(myfile,ent_qq_interval,'A39:J39')
% xlswrite(myfile,ent_raw_data,'A40:J40')
% xlswrite(myfile,ent_wavelet,'A41:J41')
 





