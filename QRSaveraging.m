function [qrs, qrs_avg, qrs_removed] = QRSaveraging(val, q_points, t_end,fs)

% This function finds the average QRS complex for a given input signal
% Inputs:
%       val = input ECG signal
%       q_points = Detected locations of Q points
%       t_end = Detected locations of end of T wave
% Outputs:
%       qrs = the matrix with all QRS-T complexes, 1 in each row
%       qrs_avg = the average of all QRS-T complexes
%       qrs_removed = signal with GRS-T segments cancelled out

% find length of each QRS-T segment
for i = 1 : length(t_end)
    qrs_size(i) = t_end(i) - q_points(i);
end
% find longest QRS-T portion
max_size = max(qrs_size);

% create a vector for average QRS-T portion the length of largest QRS-T
qrs = zeros(length(t_end), max_size);

% separate out each QRS-T segment
for i = 1 : length(t_end)
    cur_len = t_end(i) - q_points(i) + 1;   % find length of current QRS-T
    qrs(i,1:cur_len) = val(q_points(i):t_end(i));   % place QRS-T in new row
end

% find average of all QRS-T segments
qrs_avg = mean(qrs(1:length(t_end),:));

% find size of QRS-T matrix
[m,n] = size(qrs);

% duplicate input signal for removal of QRS-T segments
qrs_removed = val;

% replaces each QRS-T segment in input signal with the difference between
% that segment and the average QRS-T segment
for i = 1 : m
    qrs_removed(q_points(i):q_points(i) + n - 1) = qrs(i,:) - qrs_avg;
end

% removes any partial QRS-T segments that were not cancelled out
%qrs_removed(end-300:end) = 0;
%qrs_removed(1:300) = 0;

% histo = hist(qrs_removed, min(qrs_removed):max(qrs_removed)); % Form a histogram of the input
% qrs_p = histo/sum(histo); % Estimate the probabilities (Pi)
% %aa_p = abs(qrs_p(:) / sum(qrs_p));
% qrs_p = qrs_p(qrs_p ~= 0);
% qrs_ent = -sum(qrs_p.*log2(qrs_p));

% time vector for average QRS-T segment
t_avg = 0:1/fs:length(qrs_avg)/fs-1/fs;

% plot average QRS-T segment
figure
plot(t_avg,qrs_avg);
title('Average QRS Complex'); xlabel('Time (s)'); ylabel('Amplitude');

% time vector for signal with QRS-T segments removed
t = 0:1/fs:length(qrs_removed)/fs - 1/fs;

% plot signal wit QRS-T segments removed
figure
plot(t,qrs_removed);
title('Data with QRS Complexes cancelled out'); xlabel('Time (s)'); ylabel('Amplitude');
    