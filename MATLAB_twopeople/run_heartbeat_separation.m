% Base script for heartbeat separation algorithm, i.e. with known heartbeat
% times from ECG.

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
t_start = 6E6;
t_stop = t_start+1E6;
offset = 300;

[time, observ, beats1, beats2] = load_and_calibrate(model.K, '../data/data16.mat', 'calibration.mat', t_start+offset, t_stop+offset, 0);

%% Plot data
figure, hold on
plot(time, observ)
plot(beats1(:,1), beats1(:,2),'g*')
plot(beats2(:,1), beats2(:,2),'bo')

%% Testing
heartbeat_separation(display, algo, model, time, observ, {beats1, beats2});
heartbeat_preliminary_timing(display, algo, model, time, observ);

%% Evaluation

%% Plot graphs

% if (~exist('flags.batch', 'var')||(~flags.batch)) && display.plot_after
%     
%     figure, surf(mean(cat(3,ps.rb_mn),3)), shading interp
%     figure, hold on, arrayfun(@(x) plot(x.cp_time,x.cp_param(1,:)), ps); 
%                      arrayfun(@(x) plot(x.cp_time(2:end), diff(x.cp_time), 'r'), ps);
%     figure, hold on, arrayfun(@(x) plot(x.cp_time,x.cp_param(2,:)), ps);
%     figure, hold on, plot(time, observ, 'r'),
%                      plot(time, mean(cat(1,ps.signal_mn),1), 'b'),
%                      plot(time, mean(cat(1,ps.signal_mn),1)+2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
%                      plot(time, mean(cat(1,ps.signal_mn),1)-2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
%                      plot(unique(cat(2,ps.cp_time)), zeros(size(unique(cat(2,ps.cp_time)))), 'g*');
%     figure, hold on, plot(time, observ-mean(cat(1,ps.signal_mn),1), 'm');
%     figure, hold on, plot(time, mean(cat(1,ps.clut),1));
%     
% end
