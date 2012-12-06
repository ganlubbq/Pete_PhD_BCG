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
    display.plot_during = false;
    display.plot_after_frame = 0;
    if display.plot_during
        display.h_pf(1) = figure;
        display.h_pf(2) = figure;
        display.h_pf(3) = figure;
        display.h_pf(4) = figure;
%         display.h_pf(5) = figure;
%         display.h_pf(6) = figure;
    end
    display.plot_after = true;
    
    
end

%% Load some data
[time, observ] = load_and_calibrate(model.K, '../data/F_data1.mat', 'calibration.mat');

%% Run the particle filter
[pf, ps] = hearbeat_vrpf(display, algo, model, time, observ);

%% Evaluation

%% Plot graphs

if (~exist('flags.batch', 'var')||(~flags.batch)) && display.plot_after
    
    figure, hold on, plot(time, observ), for ii = 1:algo.Nf, plot(ps(ii).cp_time, zeros(size(ps(ii).cp_time)), 'g*'), end
    figure, surf(mean(cat(3,ps.rb_mn),3)), shading interp
%     figure, hold on, cellfun(@(x,y) plot(x,y), pf_cp, pf_p);
%     figure, hold on, cellfun(@(x,y) plot(x,y), pf_cp, pf_a);
%     figure, hold on, cellfun(@(x,y) plot(x,y), pf_cp, pf_b);
%     figure, hold on, cellfun(@(x) plot(x(2:end),diff(x)), pf_cp);
%     figure, hold on, cellfun(@(x,y) plot(x,y), pf_cp, pf_p), cellfun(@(x,y,z) plot(x,y+z,'m'), pf_cp, pf_p, pf_b), cellfun(@(x) plot(x(1:end-1),diff(x),'r'), pf_cp)
%     figure, hold on, surf(time, 1:model.dw, rb_est), shading interp
%     figure, hold on, plot(time, reconstructed,'b'), plot(time, observ, 'r'), for ii = 1:algo.Nf, plot(pf_cp{ii}, zeros(size(pf_cp{ii})), 'g*'); end
%     figure, hold on, plot(time, observ-reconstructed)
%     figure, hold on, plot(time, mean(pf_clut));
%     figure, hold on, plot(time, lhood_est);
%     
    
end
