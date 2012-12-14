function [ time, F ] = load_and_calibrate(K, datafilename, calibfilename)

% Set things
% t_start = 2.8E6;    % Time index to start at
% t_stop = 2.9E6;     % Time index to stop at
t_start = 6.05E6;    % Time index to start at
t_stop = 7.05E6;     % Time index to stop at
dsr = 10;            % Down sampling ratio - sets low pass filter cut-off (5 --> 30Hz)

% Constants
fs = 300;

% Derivative constants
fhigh = fs/(2*dsr);
flow = 0.2;

% Load data
load(datafilename);

% Time overflow
t_diff = diff(time);
jump_idx = find(t_diff<0,1);
t_diff(jump_idx) = mod(time(jump_idx+1)-time(jump_idx),2^23);
jump_idx = find(t_diff<0,1);
t_diff(jump_idx) = 1;
time = cumsum([time(1); t_diff]);
clear t_diff;

% Pick out a nice chunk to work on
t_idx = find((time>t_start)&(time<t_stop));
time = time(t_idx);
value1 = value1(t_idx);
value2 = value2(t_idx);

%%

% Low pass filter to remove 50Hz
% b = fir1(300, [flow fhigh]/(fs/2));
b = fir1(300, fhigh/(fs/2)) - ones(1,301)/301;
value1_filt = filtfilt(b, 1, value1);
value2_filt = filtfilt(b, 1, value2);

% Downsample
value1_filt = value1_filt(1:dsr:length(time));
value2_filt = value2_filt(1:dsr:length(time));
time = time(1:dsr:length(time));

% Calibrate
load(calibfilename)
% F1 = (value1_filt - c1).*m1;
% F2 = (value2_filt - c2).*m2;
F1 = (value1_filt).*m1;
F2 = (value2_filt).*m2;
F = F1+F2;
time = time - time(1);
time = time / fs;

% Select a little chunk
% offset = 315;       % Good clear section with 2.8-2.9E6
% offset = 6716;      % Good clear section with 2.8-2.9E6
% offset = 58;        % Difficult section with 6.05-7.05E6
% offset = 62025;     % Really difficult section with 6.05-7.05E6
% offset = 62010;     % Really difficult section with 6.05-7.05E6 not beat aligned
% offset = 9023;      % Clutter section with 6.05-7.05E6
offset = 27010;     % Clutter section with 6.05-7.05E6
F = F(offset+1:offset+K)';
F = F-mean(F);
F(1) = NaN;
time = time(offset+1:offset+K)';
time = time - time(1);
% figure, plot(F)


end
