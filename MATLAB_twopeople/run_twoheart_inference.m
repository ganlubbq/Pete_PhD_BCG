% Base script for heartbeat inference algorithm testing for two people

% Structures
% display - display options while running
% algo - algorithm parameters
% model - model parameters

%% Set Up

if ~exist('flags.batch', 'var') || (~flags.batch)
    
    clup
    dbstop if error
    % dbstop if warning
    
    %%% SETTINGS %%%
    
    % DEFINE RANDOM SEED
    rand_seed = 1;
    
    t_start = 6E6;
    t_length = 1E6;
    offset = 300;
    
    %%%%%%%%%%%%%%%%
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
    % Set flag to non-batch
    flags.batch = false;
    
    % Add pseudo-class directory to path
    addpath('pseudoclasses/');
    
    % Set model and algorithm parameters
    set_model;
    set_algo;
    
    % Set display options
    display.text = true;
    display.plot_during = false;
    display.plot_after = true;
    display.plot_after_frame = 50;
    if display.plot_during
        num_figs = 5;
        for ff = 1:num_figs
            display.h_pf(ff) = figure;
        end
    end
    
    
end

%% Load some data
[time, observ, beats1, beats2] = load_and_calibrate(model.K, '../data/data16.mat', 'calibration.mat', t_start+offset, t_start+t_length+offset, 0);

% %% Plot data
% figure, hold on
% plot(time, observ)
% plot(beats1(:,1), beats1(:,2),'g*')
% plot(beats2(:,1), beats2(:,2),'bo')

%% Testing
heartbeat_inference(display, algo, model, time, observ);