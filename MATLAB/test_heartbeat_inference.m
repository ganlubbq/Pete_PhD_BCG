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
% t_start = 2.8E6;    % Time index to start at
% t_stop = 2.9E6;     % Time index to stop at
t_start = 6.05E6;    % Time index to start at
t_stop = 7.05E6;     % Time index to stop at
% offset = 315;       % Good clear section with 2.8-2.9E6
% offset = 6716;      % Good clear section with 2.8-2.9E6
% offset = 58;        % Difficult section with 6.05-7.05E6
% offset = 62025;     % Really difficult section with 6.05-7.05E6
offset = 62010;     % Really difficult section with 6.05-7.05E6 not beat aligned
% offset = 9023;      % Clutter section with 6.05-7.05E6
% offset = 27010;     % Clutter section with 6.05-7.05E6
[time, observ] = load_and_calibrate(model.K, '../data/F_data1.mat', 'calibration.mat', t_start, t_stop, offset);

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
