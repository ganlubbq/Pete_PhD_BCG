% Base script for heartbeat inference algorithm testing for two people

% Structures
% display - display options while running
% algo - algorithm parameters
% model - model parameters

%% Set Up

if ~exist('flags', 'var') || ~isfield(flags,'batch') || (~flags.batch)
    
    clup
    dbstop if error
    % dbstop if warning
    
    % Set flag to non-batch
    flags.batch = false;
    
    % Add pseudo-class directory to path
    addpath('pseudoclasses/');
    
    % Set model and algorithm parameters
    set_model;
    set_algo;
    
    %%% SETTINGS %%%
    
    % DEFINE RANDOM SEED
    rand_seed = 0;
    
    t_length = 1E6;
    offset = 100;
    
    if model.np == 1
        data_file = '../data/one_person/data12.mat' % Default 3
        t_start = 3E6;
    elseif model.np == 2
        data_file = '../data/two_person/ash_phil.mat' %'../data/two_person/data16.mat' % Default 16
        t_start = 1E6;
    end
    
    %%%%%%%%%%%%%%%%
    
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
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
else
    
    disp('_____RUNNING IN BATCH MODE_____');
    
end

%% Load some data
[time, observ, beats1, beats2] = load_and_calibrate(model.K, data_file, 'calibration.mat', t_start, t_start+t_length, offset);

%% Testing
[ps, ess] = heartbeat_inference(display, algo, model, time, observ);

%% Plot data
if display.plot_after
    figure, hold on
    plot(time, observ)
    plot(beats1(:,1), beats1(:,2),'g*')
    plot(beats2(:,1), beats2(:,2),'go')
end

%% Plot Results

if display.plot_after
    ii=1;
    [H, Y] = heartbeat_obsmat(algo, model, time, observ, ps(ii).beat);
    [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time, observ, ps(ii).beat );
    x_vec = H*wf_mn;
    x = reshape(x_vec, model.K, 4)';
    
    w_av{1} = reshape(wf_mn(1:4*model.dw), model.dw, 4);
    figure, plot(w_av{1});
    if model.np == 2
        w_av{2} = reshape(wf_mn(4*model.dw+1:end), model.dw, 4);
        figure, plot(w_av{2});
    end
    
    
    figure, hold on
    plot(time, x);
    for ii = 1:length(ps), plot(ps(ii).beat(1).time, [diff(ps(ii).beat(1).time) 0], 'b*'); end
    for ii = 1:length(ps), plot(ps(ii).beat(1).pre_time, ps(ii).beat(1).time(1)-ps(ii).beat(1).pre_time, 'b*'); end
    if model.np == 2
        for ii = 1:length(ps), plot(ps(ii).beat(2).time, [diff(ps(ii).beat(2).time) 0], 'bo'); end
        for ii = 1:length(ps), plot(ps(ii).beat(2).pre_time, ps(ii).beat(2).time(1)-ps(ii).beat(2).pre_time, 'bo'); end
    end
    
    figure, hold on
    plot(time, x-observ);
    
end

%% How good it it?

if display.plot_after
    
    figure, hold on
    plot(beats1(:,1), beats1(:,2),'g*')
    for ii = 1:length(ps), plot(ps(ii).beat(1).time(1:end-1), diff(ps(ii).beat(1).time), 'b*'); end
    % for ii = 1:length(ps), plot(ps(ii).beat(1).time, ps(ii).beat(1).param, 'b*-'); end
    if model.np == 2
        plot(beats2(:,1), beats2(:,2),'go')
        for ii = 1:length(ps), plot(ps(ii).beat(2).time(1:end-1), diff(ps(ii).beat(2).time), 'bo'); end
        %     for ii = 1:length(ps), plot(ps(ii).beat(2).time, ps(ii).beat(2).param, 'bo-'); end
    end
    ylim([0 2])
    
end