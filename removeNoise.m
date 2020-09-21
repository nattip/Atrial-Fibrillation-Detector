function ecg = removeNoise(f1, f2, fs, t, val)

% This function creates a bandpass filter for given corner frequencies
% Inputs:
%       f1 = low end corner frequency
%       f2 = high end corner frequency
%       fs = samlping frequency of signal
%       t = time vector for signal
%       val = signal being filtered
% Outputs:
%       ecg = filtered signal
%       plots filtered signal

% filter to remove baseline wander and noise from muscles
cutoff =[f1 f2]*2/fs;                 % cutoff based on fs
order = 3;                            % order of 3 less processing
[a,b] = butter(order,cutoff);         % create filter coefficients
ecg = filtfilt(a,b,val);              % filter data
ecg = ecg/ max( abs(ecg));            % zero mean the data

% plot filtered data
figure
plot(t,ecg);
title('Baseline wander and noise removed'); xlabel('Time (s)'); ylabel('Amplitude');