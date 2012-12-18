% Plot graphs for preliminary report - December 2012

%% Set Up

clup
dbstop if error

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

display.text = true;
display.plot_during = false;

% Set model and algorithm parameters
set_model;
set_algo;
model.K = 500;

% Path
graph_path = '../graphs/';

% Chunk of the night
t_start = 6.05E6;    % Time index to start at
t_stop = 7.05E6;     % Time index to stop at

suffixes = {'easy', 'ok', 'hard', 'clutter'};
offsets = {57110, 8903, 62010, 9250};

for ii = 1:4
    
    offset = offsets{ii};
    suffix = suffixes{ii};
    
    if ii == 4
        model.K = 1500;
    end
    
    % Load data
    [time, observ] = load_and_calibrate(model.K, '../data/F_data1.mat', 'calibration.mat', t_start, t_stop, offset);
    
    % Plot signal
    figure(1), clf, hold on
    plot(time, observ, 'r');
    xlabel('time / s'); ylabel('Force / N');
    ylim([-3, 3]);
    set(1, 'units', 'inches');
    set(1, 'position', [2 2 9 3]);
    
    % Save plot
    filename = [graph_path 'BCG_signal_' suffix '.pdf'];
    export_pdf(1, filename, 6,2,'inches');
    close(1);
    
    % Run algorithm
    [pf, ps] = hearbeat_vrpf(display, algo, model, time, observ);
    
    % Plot reconstructed signal and timings
    figure(1), clf, hold on
    plot(time, observ, 'r'),
    plot(time, mean(cat(1,ps.signal_mn),1), 'b'),
    plot(time, mean(cat(1,ps.signal_mn),1)+2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
    plot(time, mean(cat(1,ps.signal_mn),1)-2*sqrt(mean(cat(1,ps.signal_vr),1)), ':b'),
    plot(unique(cat(2,ps.cp_time)), zeros(size(unique(cat(2,ps.cp_time)))), 'g*');
    xlabel('time / s'); ylabel('Force / N');
    ylim([-3, 3]);
    set(1, 'units', 'inches');
    set(1, 'position', [2 2 9 3]);
    
    % Save plot
    filename = [graph_path 'BCG_timing_reconstruction_' suffix '.pdf'];
    export_pdf(1, filename, 6,2,'inches');
    close(1);
    
    % Plot periods
    figure(1), clf, hold on
    arrayfun(@(x) plot(x.cp_time,x.cp_param(1,:)), ps); 
    arrayfun(@(x) plot(x.cp_time(2:end), diff(x.cp_time), 'r'), ps);
    xlabel('time / s'); ylabel('Heartbeat Period / s');
    ylim([0, 2]);
    set(1, 'units', 'inches');
    set(1, 'position', [2 2 9 3]);
    
    % Save plot
    filename = [graph_path 'BCG_periods_' suffix '.pdf'];
    export_pdf(1, filename, 6,2,'inches');
    close(1);
    
    if ii == 2
        % Template evolution
        figure(1), clf, hold on
        surf(time, 1:model.dw, mean(cat(3,ps.rb_mn),3)), shading interp
        set(gca, 'cameraposition',[30,0,4])
        saveas(1,[graph_path 'template_surface_' suffix '.png'], 'png');
        close(1);
        
        % Initial template
        figure(1), clf, hold on
        plot(1:model.dw, model.w_prior_mn);
        filename = [graph_path 'initial_template_' suffix '.pdf'];
        export_pdf(1, filename, 4,3,'inches');
        close(1);
    end
    
end