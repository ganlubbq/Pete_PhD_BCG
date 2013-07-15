function [ time, bcg_out, beats1, beats2 ] = load_and_calibrate(K, datafilename, calibfilename, t_start, t_stop, offset)

% Set things
dsr = 8;            % Down sampling ratio - sets low pass filter cut-off (5 --> 30Hz)

% Constants
fs = 240;

% Derivative constants
fhigh = fs/(2*dsr);
flow = 0.2;

% Load data
load(datafilename);

% Pick out a nice chunk to work on
t_idx = find((t>t_start)&(t<t_stop));
bcg = F(t_idx,:);
ecg = E(t_idx,:);

% Extract beats from ECG
[beats1] = ecgdetect(ecg(:,1),fs,1.5,0);
[beats2] = ecgdetect(ecg(:,2),fs,1.5,0);

% Low pass filter to remove 50Hz
% b = fir1(300, [flow fhigh]/(fs/2));
b = fir1(300, fhigh/(fs/2)) - ones(1,301)/301;
bcg_filt = zeros(size(bcg));
ecg_filt = zeros(size(ecg));
for ss = 1:4
    bcg_filt(:,ss) = filtfilt(b, 1, bcg(:,ss));
end
for ss = 1:2
    ecg_filt(:,ss) = filtfilt(b, 1, ecg(:,ss));
end

% Downsample
bcg_filt = bcg_filt(1:dsr:length(t_idx),:);
ecg_filt = ecg_filt(1:dsr:length(t_idx),:);

% Calibrate
% load(calibfilename)


% Select a little chunk
bcg_out = bcg_filt(offset+1:offset+K,:);
ecg_out = ecg_filt(offset+1:offset+K,:);

% Remove DC offset from bcg
bcg_out = bsxfun(@minus, bcg_out, mean(bcg_out))';

% Time
time = (0:K-1)/(fs/dsr);

% Chop beats
beats1(beats1(:,1)>time(end),:) = [];
beats2(beats2(:,1)>time(end),:) = [];

end
