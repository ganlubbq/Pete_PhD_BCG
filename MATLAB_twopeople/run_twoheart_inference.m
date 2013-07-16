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
    
    t_start = 6E6;%3E6;%
    t_length = 1E6;
    offset = 110;
    
    data_file = '../data/data16.mat'
    
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
[time, observ, beats1, beats2] = load_and_calibrate(model.K, data_file, 'calibration.mat', t_start, t_start+t_length, offset);

%% Testing
[ps] = heartbeat_inference(display, algo, model, time, observ);

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
for ii = 1:algo.Nf, plot(ps(ii).beat(1).time, [diff(ps(ii).beat(1).time) 0], 'b*'); end
for ii = 1:algo.Nf, plot(ps(ii).beat(2).time, [diff(ps(ii).beat(2).time) 0], 'bo'); end

%% How good it it?
figure, hold on
plot(beats1(:,1), beats1(:,2),'g*')
plot(beats2(:,1), beats2(:,2),'go')
for ii = 1:algo.Nf, plot(ps(ii).beat(1).time, [diff(ps(ii).beat(1).time) 0], 'b*'); end
for ii = 1:algo.Nf, plot(ps(ii).beat(2).time, [diff(ps(ii).beat(2).time) 0], 'bo'); end