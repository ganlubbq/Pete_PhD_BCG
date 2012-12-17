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
    rand_seed = 1;
    
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
    display.plot_during = false;
    display.plot_after_frame = 50;
    if display.plot_during
        display.h_pf(1) = figure;
        display.h_pf(2) = figure;
        display.h_pf(3) = figure;
        display.h_pf(4) = figure;
        display.h_pf(5) = figure;
%         display.h_pf(6) = figure;
    end
    display.plot_after = true;
    
    
end

%% Load some data
[time, observ] = load_and_calibrate(model.K, '../data/F_data1.mat', 'calibration.mat');
% offset = 9002;
% time(1:offset) = [];
% time = time - time(1);
% observ(1:offset) = [];
% model.K = model.K - offset;

%% Run the particle filter
[pf, ps] = hearbeat_vrpf(display, algo, model, time, observ);

%% Evaluation

%% Plot graphs

if (~exist('flags.batch', 'var')||(~flags.batch)) && display.plot_after
    
    figure, surf(mean(cat(3,ps.rb_mn),3)), shading interp
    figure, hold on, arrayfun(@(x) plot(x.cp_time,x.cp_param(1,:)), ps); 
                     arrayfun(@(x) plot(x.cp_time(2:end), diff(x.cp_time), 'r'), ps);
    figure, hold on, arrayfun(@(x) plot(x.cp_time,x.cp_param(2,:)), ps);
    figure, hold on, plot(time, observ, 'r'),
                     plot(time, mean(cat(1,ps.signal_mn),1), 'b'),
                     plot(time, mean(cat(1,ps.signal_mn),1)+2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
                     plot(time, mean(cat(1,ps.signal_mn),1)-2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
                     plot(unique(cat(2,ps.cp_time)), zeros(size(unique(cat(2,ps.cp_time)))), 'g*');
    figure, hold on, plot(time, observ-mean(cat(1,ps.signal_mn),1), 'm');
    figure, hold on, plot(time, mean(cat(1,ps.clut),1));
    
end
