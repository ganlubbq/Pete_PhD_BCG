filename = 'section_tests'
wid = 10;
hei = 5;

np1_test_num = 5;
np2_test_num = 1;

%% two people
addpath('standard8_NP2_NF200')

% Load data
load([filename '_' num2str(np2_test_num)]);

% Plot data
h_np2_obs = figure; hold on
plot(time, observ)
xlabel('time (s)')
ylabel('force (N)')
ylim([-1 1]);
export_pdf(h_np2_obs, 'np2_data.pdf', wid, hei);

% Plot timing results
h_np2_timing = figure; hold on
for ii = 1:length(ps), plot(ps(ii).beat(1).time(1:end-1), diff(ps(ii).beat(1).time), 'b*:'); end
for ii = 1:length(ps), plot(ps(ii).beat(2).time(1:end-1), diff(ps(ii).beat(2).time), 'bo:'); end
plot(beats1(:,1), beats1(:,2),'g*-')
plot(beats2(:,1), beats2(:,2),'go-')
xlabel('time (s)')
xlabel('duration (s)')
xlim([0 model.K/model.fs]);
ylim([0.5 1.5]);
export_pdf(h_np2_timing, 'np2_timing.pdf', wid, hei);

% Plot ESSs
h_np2_ess = figure; hold on
plot(ess)
xlabel('frame number')
ylabel('ESS')
ylim([0 200])
export_pdf(h_np2_ess, 'np2_ess.pdf', wid, hei);

rmpath('standard8_NP2_NF200')

%% one person
addpath('standard8_NP1_NF200')

% Load data
load([filename '_' num2str(np2_test_num)]);

% Plot data
h_np1_obs = figure; hold on
plot(time, observ)
xlabel('time (s)')
ylabel('force (N)')
ylim([-1 1]);
export_pdf(h_np1_obs, 'np1_data.pdf', wid, hei);

% Plot timing results
h_np1_timing = figure; hold on
for ii = 1:length(ps), plot(ps(ii).beat(1).time(1:end-1), diff(ps(ii).beat(1).time), 'b*:'); end
plot(beats1(:,1), beats1(:,2),'g*-')
xlim([0 model.K/model.fs]);
xlabel('time (s)')
xlabel('duration (s)')
ylim([0.5 1.5]);
export_pdf(h_np1_timing, 'np1_timing.pdf', wid, hei);

% Plot ESSs
h_np1_ess = figure; hold on
plot(ess)
xlabel('frame number')
ylabel('ESS')
ylim([0 200])
export_pdf(h_np1_ess, 'np1_ess.pdf', wid, hei);

rmpath('standard8_NP1_NF200')
