function [aa, coeffs, we] = wavelets(val, t, fs)

% This function computes 9th order wavelets.
% Inputs:
%       val = input signal
%       t = time vector for input signal
%       fs = sampling frequency of input signal
% Outputs:
%       aa = atrial activity of input signal. All content less than 10 Hz
%       aa_ent = entropy of atrial activity
%       aa_went = wavelet entropy of atrial activity
%       plots the content of each of the 9 wavelets
% 

% obtain wavelet coefficients
[c,l] = wavedec(val, 9, 'dmey');

% determine the frequency ranges of each wavelet for labelling
range_start = fs/2;
r1 = range_start / 2;
r2 = r1 / 2;
r3 = r2 / 2;
r4 = r3 / 2;
r5 = r4 / 2;
r6 = r5 / 2;
r7 = r6 / 2;
r8 = r7 / 2;
r9 = r8 / 2;

% getting detailed and approximate coefficients
a9 = wrcoef('a', c, l, 'dmey', 9);
d1 = wrcoef('d', c, l, 'dmey', 1);
d2 = wrcoef('d', c, l, 'dmey', 2);
d3 = wrcoef('d', c, l, 'dmey', 3);
d4 = wrcoef('d', c, l, 'dmey', 4);
d5 = wrcoef('d', c, l, 'dmey', 5);
d6 = wrcoef('d', c, l, 'dmey', 6);
d7 = wrcoef('d', c, l, 'dmey', 7);
d8 = wrcoef('d', c, l, 'dmey', 8);
d9 = wrcoef('d', c, l, 'dmey', 9);

% plot and label all wavelets
figure
subplot(11,1,1)
plot(t, val, 'r');
title('Wavelets');
subplot(11,1,2)
plot(t,a9); ylabel('a9');xlabel(['0 - ', num2str(r9), ' Hz']);
subplot(11,1,3)
plot(t,d9); ylabel('d9');xlabel([num2str(r9), ' - ', num2str(r8), ' Hz']);
subplot(11,1,4)
plot(t,d8); ylabel('d8');xlabel([num2str(r8), ' - ', num2str(r7), ' Hz']);
subplot(11,1,5)
plot(t,d7); ylabel('d7');xlabel([num2str(r7), ' - ', num2str(r6), ' Hz']);
subplot(11,1,6)
plot(t,d6); ylabel('d6');xlabel([num2str(r6), ' - ', num2str(r5), ' Hz']);
subplot(11,1,7)
plot(t,d5); ylabel('d5');xlabel([num2str(r5), ' - ', num2str(r4), ' Hz']);
subplot(11,1,8)
plot(t,d4); ylabel('d4');xlabel([num2str(r4), ' - ', num2str(r3), ' Hz']);
subplot(11,1,9)
plot(t,d3); ylabel('d3');xlabel([num2str(r3), ' - ', num2str(r2), ' Hz']);
subplot(11,1,10)
plot(t,d2); ylabel('d2');xlabel([num2str(r2), ' - ', num2str(r1), ' Hz']);
subplot(11,1,11)
plot(t,d1); ylabel('d1'); xlabel([num2str(r1), ' - ', num2str(range_start), 'Hz']);


coeffs = [a9; d9; d8];  % coefficients that contribute to atrial activity

[nj, nk] = size(coeffs);

bottom = 0;     % set denominator to zero

% cycle through all coefficients
% double summation of C(j,k)^2 for k=1:Pj and j=1:N
for j = 1:nj 
    for k = 1:nk
        bottom = bottom + coeffs(j,k)^2;        
    end
end

top = zeros(1,nj);  % set numerator to all zeros

%find the numerator for each value of j
for j = 1:nj
    for k = 1:nk
        top(j) = top(j) + coeffs(j,k)^2; % numerator = sum of all k for each j
        E(j) = top(j)/ bottom;  % find E(j) 
    end
end

% calculate wavelet entropy
we = -sum(E .* log2(E));

% recreate atrial activity by including only the low frequency bands
aa = d8+d9+a9+d7;
