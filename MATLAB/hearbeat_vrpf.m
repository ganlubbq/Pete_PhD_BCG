function [ pf, ps ] = hearbeat_vrpf( display, algo, model, time, observ )
%HEARBEAT_PF Run a particle filter to infer heartbeats from a BCG or ECG
%signal.

% The changepoint paramters for estimation are the beat periods and
% adjustments. An additional changepoint parameter, the beat waveform
% is Rao-blackwellised.

% Make local copies of useful numbers
Nf = algo.Nf;
K = model.K;

% Number of particle filter steps we need to run
M = ceil(K/algo.S)+1;

% Create particle filter structure arrays
[ps] = init_ps(algo, model);

if display.text
    fprintf(1, 'Particle filter iteration %u.\n', 1);
end

%%% Initial particle filter step

% Time indexes
L = algo.L;
S = algo.S;
kk = 1 - S;

% Inialise strucutre
pf = init_pf(algo, model, algo.L);

% Loop through particles
for ii = 1:Nf
    
    % Sample parameter from prior
    [cp_param, ~] = heartbeat_cpparamprior(model);
    pf(ii).pre_cp_param(:) = cp_param;
    
    % Linear state and likelihood
    pf(ii).win_rb_mn(:,S) = model.w_prior_mn;
    pf(ii).win_rb_vr(:,:,S) = model.w_prior_vr;
    rb_mn = model.w_prior_mn;
    rb_vr = model.w_prior_vr;
    for ll = S+1:L
        
        H = heartbeat_interpolation(algo, model, time(kk+ll), 0);
        [rb_mn, rb_vr, ~,~,~, lhood] = log_kf_update(rb_mn, rb_vr, observ(kk+ll), H, model.y_obs_vr);
        
        pf(ii).win_rb_mn(:,ll) = rb_mn;
        pf(ii).win_rb_vr(:,:,ll) = rb_vr;
        pf(ii).win_obslhood(ll) = lhood;
        
    end
    
end

% Loop through time
for mm = 2:M
    
    last_pf = pf;
    last_ps = ps;
    
    % Time indexes
    kk = kk + algo.S;
    L = min(algo.L, K-kk);
    S = min(algo.S, K-kk);
    
    if display.text
        fprintf(1, 'Particle filter iteration %u of %u, using a window of %f to %f (indexes %u to %u)\n', mm, M, time(kk), time(kk+L), kk, kk+L);
    end
    
%     % Diagnostics
    diagnostic_lastest_cp_time = zeros(1,Nf);
    diagnostic_lastest_cp_param = zeros(model.dp,Nf);
%     diagnostic_last_clut = zeros(1,Nf);
    
    % Initialise pf frame
    pf = init_pf(algo, model, L);
    
    % Sample ancestors
    sampling_weight = [last_pf.weight];
    ancestor = sample_weights(algo, sampling_weight, Nf);
    
    % Loop through particles
    for ii = 1:Nf
        
        % Get ancestor index
        a_idx = ancestor(ii);
        pf(ii).ancestor = a_idx;
        
        % Copy forward
        if ~isempty(last_pf(a_idx).win_cp_time) && (last_pf(a_idx).win_cp_time<time(kk))
            pf(ii).pre_cp_time = last_pf(a_idx).win_cp_time;
            pf(ii).pre_cp_param = last_pf(a_idx).win_cp_param;
        else
            pf(ii).pre_cp_time = last_pf(a_idx).pre_cp_time;
            pf(ii).pre_cp_param = last_pf(a_idx).pre_cp_param;
        end
        pf(ii).pre_rb_mn = last_pf(a_idx).win_rb_mn(:,S);
        pf(ii).pre_rb_vr = last_pf(a_idx).win_rb_vr(:,:,S);
        pf(ii).pre_clut = last_pf(a_idx).win_clut(S);
        
        % Sample changepoints in the window
        [cp_time, cp_param, new_trans_prob]  = heartbeat_cptransition(model, pf(ii).pre_cp_time, pf(ii).pre_cp_param, time(kk), time(kk+L));
        pf(ii).win_cp_time = cp_time;
        pf(ii).win_cp_param = cp_param;
        
        % Find old transition prob
        last_cp_time = last_pf(a_idx).pre_cp_time;
        last_cp_param = last_pf(a_idx).pre_cp_param;
        if last_pf(a_idx).win_cp_time < time(kk)
            last_cp_time = last_pf(a_idx).win_cp_time;
            last_cp_param = last_pf(a_idx).win_cp_param;
        end
        [~, ~, old_trans_prob]  = heartbeat_cptransition(model, last_cp_time, last_cp_param, time(kk), time(kk+L));
        
        % Loop through observations
        last_cp_time = pf(ii).pre_cp_time;
        rb_mn = pf(ii).pre_rb_mn;
        rb_vr = pf(ii).pre_rb_vr;
        for ll = 1:L
            
            % Update changepoint if we've past one
            if time(kk+ll) > pf(ii).win_cp_time
                last_cp_time = pf(ii).win_cp_time;
                last_cp_param = pf(ii).win_cp_param;
            end
            
            % interpolation and Kalman filtering
            H = heartbeat_interpolation(algo, model, time(kk+ll), last_cp_time);
            [rb_mn, rb_vr, ~, s_mn, s_vr, lhood] = log_kf_update(rb_mn, rb_vr, observ(kk+ll), H, model.y_obs_vr);
            
            pf(ii).win_rb_mn(:,ll) = rb_mn;
            pf(ii).win_rb_vr(:,:,ll) = rb_vr;
            pf(ii).win_obslhood(ll) = lhood;
            pf(ii).win_signal_mn(ll) = s_mn;
            pf(ii).win_signal_vr(ll) = s_vr;
        end
        
        % Likelihoods
        new_lhood = sum(pf(ii).win_obslhood);
        old_lhood = sum(last_pf(a_idx).win_obslhood(S+1:end));
        
        % Weight
        pf(ii).weight = (new_lhood - old_lhood);
        
        % Diagnostics
        diagnostic_lastest_cp_time(ii) = last_cp_time;
        diagnostic_lastest_cp_param(:,ii) = last_cp_param;
    end
    
    assert(~any(isnan(pf(mm).weight)));
    assert(~all(isinf(pf(mm).weight)));
    
    % Particle smoother
    for ii = 1:Nf
        
        a_idx = pf(ii).ancestor;
        
        ps(ii) = last_ps(a_idx);
        if pf(ii).win_cp_time < time(kk+S);
            ps(ii).Ncp = ps(ii).Ncp + 1;
            ps(ii).cp_time(ps(ii).Ncp) = pf(ii).win_cp_time;
            ps(ii).cp_param(:,ps(ii).Ncp) = pf(ii).win_cp_param;
        end
        ps(ii).rb_mn(:,kk+1:kk+S) = pf(ii).win_rb_mn(:,1:S);
        ps(ii).signal_mn(kk+1:kk+S) = pf(ii).win_signal_mn(1:S);
        ps(ii).signal_vr(kk+1:kk+S) = pf(ii).win_signal_vr(1:S);
        
    end
    
    % Diagnostics
    if display.plot_during && (M>display.plot_after_frame)
        figure(display.h_pf(1)); clf; hold on; hist(diagnostic_lastest_cp_time, 100);
        figure(display.h_pf(2)); clf; hold on; hist(diagnostic_lastest_cp_param(1,:), 100);
        figure(display.h_pf(3)); clf; hold on; hist(diagnostic_lastest_cp_param(2,:), 100);
%         figure(display.h_pf(4)); clf; hold on; hist(diagnostic_lastest_cp_param(3,:), 100);
        figure(display.h_pf(4)); clf; hold on; plot( cell2mat(arrayfun(@(x) {x.win_rb_mn(:,end)}, pf')) );
%         figure(display.h_pf(5)); clf; hold on; hist(diagnostic_last_clut, 100);
        pause
    end
    
end

% Resample smoother samples
sampling_weight = [pf.weight];
ancestor = sample_weights(algo, sampling_weight, Nf);
ps = ps(ancestor);

end






function [ps] = init_ps(algo, model)
% Initialise particle smoother structure

Nf = algo.Nf;
K = model.K;

ps = struct('Ncp', cell(Nf,1), ...                  Number of changepoints over all time
            'cp_time', cell(Nf,1), ...              Changepoint times ...
            'cp_param', cell(Nf,1), ...             ... and the associated parameters
            'rb_mn', cell(Nf,1), ...                Mean of the Rao-Blackwellised bit over all time
            'signal_mn', cell(Nf,1), ...            Predicted signal mean over all time
            'signal_vr', cell(Nf,1));%              Predicted signal variance over all time

[ps.Ncp] = deal(0);
[ps.cp_time] = deal(zeros(1,0));
[ps.cp_param] = deal(zeros(model.dp,0));
[ps.rb_mn] = deal(zeros(model.dw,model.K));

end

function [pf] = init_pf(algo, model, L)
% Initialise particle smoother structure

Nf = algo.Nf;
K = model.K;
M = ceil(K/algo.S)+1;

pf = struct('pre_cp_time', cell(Nf,1), ...          Most recent changepoint to occur before the window ...
            'pre_cp_param', cell(Nf,1), ...         ... and the corresponding parameters
            'pre_rb_mn', cell(Nf,1), ...            Mean of the Rao-Blackwellised bit before the window
            'pre_rb_vr', cell(Nf,1), ...            (Co)variance of the Rao-Blackwellised bit before the window
            'pre_clut', cell(Nf,1), ...             Clutter indicator variable for the most recent observation before the window
            'win_Ncp', cell(Nf,1), ...              Number of changepoints in the window
            'win_cp_time', cell(Nf,1), ...          Changepoints occuring within the window ...
            'win_cp_param', cell(Nf,1), ...         ... and the corresponding parameters.
            'win_rb_mn', cell(Nf,1), ...            Mean of the Rao-Blackwellised bit during the window
            'win_rb_vr', cell(Nf,1), ...            (Co)variance of the Rao-Blackwellised bit during the window
            'win_clut', cell(Nf,1), ...             Clutter indicators for the observations in the window
            'win_signal_mn', cell(Nf,1),...         Predicted mean of the signal over the window
            'win_signal_vr', cell(Nf,1),...         Predicted variance of the signal over the window
            'win_obslhood', cell(Nf,1), ...         Observation likelihoods during the window
            'ancestor', cell(Nf,1), ...             Ancestor particle
            'weight', cell(Nf,1));%                 Particle weight

[pf.pre_cp_time] = deal(0);
[pf.pre_cp_param] = deal(zeros(model.dp,1));
[pf.pre_rb_mn] = deal(zeros(model.dw,1));
[pf.pre_rb_vr] = deal(zeros(model.dw,model.dw,1));
[pf.pre_clut] = deal(0);
[pf.win_cp_time] = deal(zeros(1,0));
[pf.win_cp_param] = deal(zeros(2,0));
[pf.win_rb_mn] = deal(zeros(model.dw,L));
[pf.win_rb_vr] = deal(zeros(model.dw,model.dw,L));
[pf.win_clut] = deal(zeros(1,L));
[pf.win_siganl_mn] = deal(zeros(1,L));
[pf.win_siganl_vr] = deal(zeros(1,L));
[pf.win_obslhood] = deal(zeros(1,L));
[pf.ancestor] = deal(0);
[pf.weight] = deal(0);

end


