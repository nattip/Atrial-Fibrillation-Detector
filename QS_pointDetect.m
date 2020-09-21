function [q_points, s_points, t_2, v2] = QS_pointDetect(ecg, t)

% This function detects the location of Q and S points in an ECG signal
% Inputs:
%       ecg = ECG signal on which Q and S points are being detected
%       t = time vector for input ECG signal
% Outputs:
%       q_points = array of index locations of input signal where Q points are
%       s_points = array of index location of input signal where S points are
%       t_2 = time vector for Pan-Tompkins double derivative result
%       v2 = Pan-Tomplins double derivative results

% for completion of pan-tompkins method find derivative of data
dx = ecg(2:end) - ecg(1:end-1);
dt = t(2:end) - t(1:end-1);
v1 = dx./dt;

% create time vector for derivative
t_2 = t;
t_2(end) = [];

% find derivative of derivative of data
dx = v1(2:end) - v1(1:end-1);
dt = t_2(2:end) - t_2(1:end-1);
v2 = dx./dt;

% adjust time vector again
t_2(end) = [];

% % plot the pan-tompkins result
% figure
% plot(t_2, v2);
% title('Pan-Tompkins'); xlabel('Time (s)'); ylabel('Amplitude');

% % plot all peaks
% figure
% findpeaks(v2);
% title('Peaks in Pan-Tompkins result'); xlabel('Time (s)'); ylabel('Amplitude');

% find all large peaks in pan-tompkins and the location of them
[pks_pan, locs_pan] = findpeaks(v2);    
qs_peaks = pks_pan(pks_pan > (0.3 * max(pks_pan)));
qs_locs = locs_pan(pks_pan > (0.3 * max(pks_pan)));

% find number of Q and S peaks in pean-tompkins
len_pan = length(qs_peaks);

% every other pan-thomkins peak is a Q or S point
q_points = ceil(qs_locs(1:2:end));
s_points = ceil(qs_locs(2:2:end));

% if there is more than 100 points of distance between first Q and S
% points, then the first point is actually an S point, not a Q point
if(abs(q_points(1) - s_points(1)) > 100)
    q_points(1) = [];       % get rid of first point because an S without a Q doesn't help us
    temp = s_points;        % temporarioy save the s_points
    s_points = q_points;    % make s_points equal to q_points
    q_points = temp;        % make q_points equal to temporarily saved s_points
end

qtos_dist = 30;             % an average Q to S distance is 30 for initial estimation

% finding if any Q or S points were not detected originally
for i = 1 : length(q_points)-1
    if((s_points(i) - q_points(i)) < 100)   % if less than 100 dist btwn same index Q and S
        qtos_dist(i) = s_points(i) - q_points(i);   % both correctly detected, record the distance
    else    % if distance more than 100 btwen same index Q and S, something was missed
        avg_dist = ceil(mean(nonzeros(qtos_dist))); % calculate the average Q to S distance
        qtos_dist(i) = avg_dist;
        s_points = [s_points(1:i-1), q_points(i) + avg_dist, s_points(i:end)]; % add in missing point at average distance away
        
        % we must switch following Q and S points to make room for added point
        temp = q_points(i+1:end);
        q_points = [q_points(1:i), s_points(i+1:end)];
        s_points = [s_points(1:i), q_points(i+1:end)];
    end
end 