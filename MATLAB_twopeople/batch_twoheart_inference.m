%% Setup

clup
% dbstop if error
% dbstop if warning

% Set flag to batch
flags.batch = true;


% Add pseudo-class directory to path
addpath('pseudoclasses/');

% Set model and algorithm parameters
set_model;
set_algo;

%%% SETTINGS %%%

% DEFINE RANDOM SEED
rand_seed = 0;

filename = 'section_tests'
t_length = 1E6;
offset = 110;

if model.np == 1
    data_file = '../data/data3.mat'
elseif model.np == 2
    data_file = '../data/data16.mat'
end

%%%%%%%%%%%%%%%%

% Set display options
display.text = true;
display.plot_during = false;
display.plot_after = false;

% TESTS %
test_t_start = [2.5E6, 3E6, 6E6, 8E6, 10E6, 12E6, 14E6, 16E6];
num_tests = length(test_t_start);

%% Test

for tt = 1:num_tests
    
    disp(['*** RUNNING TEST NUMBER ' num2str(tt) '. ***']);
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
    % Start time
    t_start = test_t_start(tt);
    
    % Run the algorithm
    run_twoheart_inference;
    
    % Save it all
    save([filename '_' num2str(tt)], 'data_file', 't_start', 'algo', 'model', 'time', 'observ', 'beats1', 'beats2', 'ps', 'ess');
    
end

%% Process

file_path = 'standard8_NF200';
addpath(file_path);

for tt = 1:num_tests%[3 6]%
    
    % Load data
    load([filename '_' num2str(tt)]);
    
    % Plot data
    figure, hold on
    plot(time, observ)
    plot(beats1(:,1), beats1(:,2),'g*')
    plot(beats2(:,1), beats2(:,2),'go')
    
    % Plot timing results
    figure, hold on
    plot(beats1(:,1), beats1(:,2),'g*')
    for ii = 1:length(ps), plot(ps(ii).beat(1).time(1:end-1), diff(ps(ii).beat(1).time), 'b*-'); end
    for ii = 1:length(ps), plot(ps(ii).beat(1).time, ps(ii).beat(1).param, 'c*:'); end
    if model.np == 2
        plot(beats2(:,1), beats2(:,2),'go')
        for ii = 1:length(ps), plot(ps(ii).beat(2).time(1:end-1), diff(ps(ii).beat(2).time), 'bo-'); end
        for ii = 1:length(ps), plot(ps(ii).beat(2).time, ps(ii).beat(2).param, 'co:'); end
    end
    
    % Plot ESSs
    figure, hold on
    plot(ess)
    ylim([0 100])
    
end