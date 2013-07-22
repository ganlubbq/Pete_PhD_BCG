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
    
    % Set flag to non-batch
    flags.batch = false;
    
    % Add pseudo-class directory to path
    addpath('pseudoclasses/');
    
    % Set model and algorithm parameters
    set_model;
    set_algo;
    
    %%% SETTINGS %%%
    
    % DEFINE RANDOM SEED
    rand_seed = 1;
    
    t_start = 6E6;%3E6;%
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
    
end

%% Load some data
[time, observ, beats1, beats2] = load_and_calibrate(model.K, data_file, 'calibration.mat', t_start, t_start+t_length, offset);

%% Separate tasks
if 1
    ecg_beat(1) = beat_init(model, -Inf, 0, [], beats1(:,1)', zeros(size(beats1(:,1)')) );
    if model.np > 1
        ecg_beat(2) = beat_init(model, -Inf, 0, [], beats2(:,1)', zeros(size(beats2(:,1)')) );
    end
    [wf_mn, wf_vr, H] = heartbeat_separation( display, algo, model, time, observ, ecg_beat);
    w_av{1} = reshape(wf_mn(1:4*model.dw), model.dw, 4);
    w_sd{1} = reshape(sqrt(diag(wf_vr(1:4*model.dw,1:4*model.dw))), model.dw, 4);
    figure, hold on
    plot(w_av{1}); plot(w_av{1}+2*w_sd{1}, '--'); plot(w_av{1}-2*w_sd{1}, '--');
    if model.np == 2
        w_av{2} = reshape(wf_mn(4*model.dw+1:end), model.dw, 4);
        w_sd{2} = reshape(sqrt(diag(wf_vr(4*model.dw+1:end,4*model.dw+1:end))), model.dw, 4);
        figure, hold on
        plot(w_av{2}); plot(w_av{2}+2*w_sd{2}, '--'); plot(w_av{2}-2*w_sd{2}, '--');
    end
%     [ps] = heartbeat_timeonlyinference(display, algo, model, time, observ, wf_mn, wf_vr);    
    
    model.w_prior_mn = wf_mn;
    model.w_prior_vr = wf_vr;
    
end

%% Testing
if 1
    [ps] = heartbeat_inference(display, algo, model, time, observ);
end

%% Plot data
figure, hold on
plot(time, observ)
plot(beats1(:,1), beats1(:,2),'g*')
plot(beats2(:,1), beats2(:,2),'go')

%% Plot Results

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
if model.np == 2
    for ii = 1:length(ps), plot(ps(ii).beat(2).time, [diff(ps(ii).beat(2).time) 0], 'bo'); end
end

figure, hold on
plot(time, x-observ);

%% How good it it?
figure, hold on
plot(beats1(:,1), beats1(:,2),'g*')
for ii = 1:length(ps), plot(ps(ii).beat(1).time, [diff(ps(ii).beat(1).time) 0], 'b*'); end
if model.np == 2
    plot(beats2(:,1), beats2(:,2),'go')
    for ii = 1:length(ps), plot(ps(ii).beat(2).time, [diff(ps(ii).beat(2).time) 0], 'bo'); end
end