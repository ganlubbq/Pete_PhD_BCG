% Base script for running hearbeat inference algorithm

% Structures
% display - display options while running
% algo - algorithm parameters
% model - model parameters

%% Set Up

if ~exist('flags.batch', 'var') || (~flags.batch)
    
    clup
    dbstop if error
    % dbstop if warning
    
    % DEFINE RANDOM SEED
    rand_seed = 0;
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
    % Set flag to non-batch
    flags.batch = false;
    
    % Set model and algorithm parameters
    set_model;
    set_algo;
    
    % Set display options
    display.text = true;
    display.plot_during = true;
    if display.plot_during
        display.h_pf(1) = figure;
        display.h_pf(2) = figure;
        display.h_pf(3) = figure;
        display.h_pf(4) = figure;
%         display.h_pf(5) = figure;
    end
    display.plot_after = true;
    
    
end

%% Load some data
[time, observ] = load_and_calibrate(model.K, '../data/F_data1.mat', 'calibration.mat');

%% Run the particle filter
[pf] = hearbeat_vrpf(display, algo, model, time, observ);

%% Evaluation
[cp_list, pf_cp] = process_pf(algo, model, pf);

%% Plot graphs

if (~flags.batch) && display.plot_after
    
    figure, hold on
    plot(time, observ)
    for ii = 1:algo.Nf
        plot(pf_cp{ii}, zeros(size(pf_cp{ii})), 'r*')
    end
    
end
