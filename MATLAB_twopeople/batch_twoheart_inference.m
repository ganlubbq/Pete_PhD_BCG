clup
% dbstop if error
% dbstop if warning

% Set flag to batch
flags.batch = true;
display.plotting = false;

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
test_t_start = [3E6, 5E6, 8E6, 12E6, 16E6, 19E6, 21E6, 24E6];
num_tests = length(test_t_start);

for tt = 1:num_tests
    
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