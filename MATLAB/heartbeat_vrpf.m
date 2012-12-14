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
M = ceil(K/algo.S);

% Create particle filter structure arrays
ps = init_ps(algo, model);
pf = init_pf(algo, model, algo.L);

% Initialise
for ii = 1:Nf
    pf(ii).win_rb_mn(:,algo.S) = model.w_prior_mn;
    pf(ii).win_rb_vr(:,:,algo.S) = model.w_prior_vr;
end

%%% Particle filter %%%

% Initialise time
kk = 1 - algo.S;

% Loop through time
for mm = 1:M
    
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
%     length(unique(ancestor))
    
    % Loop through particles
    for ii = 1:Nf
        
        % Get ancestor index
        a_idx = ancestor(ii);
        pf(ii).ancestor = a_idx;
        
        % Copy forward
        if ~isempty(last_pf(a_idx).win_cp_time) && (last_pf(a_idx).win_cp_time<time(kk))
            pf(ii).pre_cp_time = last_pf(a_idx).win_cp_time;
            pf(ii).pre_cp_param = last_pf(a_idx).win_cp_param;
            pf(ii).ante_cp_time = last_pf(a_idx).pre_cp_time;
            pf(ii).ante_cp_param = last_pf(a_idx).pre_cp_param;
        else
            pf(ii).pre_cp_time = last_pf(a_idx).pre_cp_time;
            pf(ii).pre_cp_param = last_pf(a_idx).pre_cp_param;
            pf(ii).ante_cp_time = last_pf(a_idx).ante_cp_time;
            pf(ii).ante_cp_param = last_pf(a_idx).ante_cp_param;
        end
        pf(ii).pre_rb_mn = last_pf(a_idx).win_rb_mn(:,S);
        pf(ii).pre_rb_vr = last_pf(a_idx).win_rb_vr(:,:,S);
        pf(ii).pre_clut = last_pf(a_idx).win_clut(S);
        last_clut = kk + find(last_pf(a_idx).win_clut(1:S)==1,1,'last');
        last_noclut = kk + find(last_pf(a_idx).win_clut(1:S)==0,1,'last');
        if isempty(last_clut)
            last_clut = last_pf(a_idx).pre_last_clut;
        end
        if isempty(last_noclut)
            last_noclut = last_pf(a_idx).pre_last_noclut;
        end
        pf(ii).pre_last_clut = last_clut;
        pf(ii).pre_last_noclut = last_noclut;
        
        %%% CHANGEPOINT SAMPLING %%%
        
        % Sample a new parameter for the preceeding changepoint (resample move)
        if isempty(pf(ii).ante_cp_time)
            [~, ppsl_pre_cp_param] = heartbeat_cpprior(model, 0, inf);
        else
            [~, ppsl_pre_cp_param] = heartbeat_cptransition(model, pf(ii).ante_cp_time, pf(ii).ante_cp_param, 0, inf);
        end
        old_tp = log(1 - invgamcdf(time(kk)-pf(ii).pre_cp_time-pf(ii).pre_cp_param, model.tau_trans_shape, model.tau_trans_scale));
        new_tp = log(1 - invgamcdf(time(kk)-pf(ii).pre_cp_time-ppsl_pre_cp_param,   model.tau_trans_shape, model.tau_trans_scale));
        ap = new_tp-old_tp;
        if log(rand)<ap
            pf(ii).pre_cp_param = ppsl_pre_cp_param;
        end
        
        % Sample changepoints in the window and find transition prob
        if isempty(pf(ii).pre_cp_time)
            [cp_time, cp_param, new_trans_prob]  = heartbeat_cpprior(model, time(kk), time(kk+L));
        else
            [cp_time, cp_param, new_trans_prob]  = heartbeat_cptransition(model, pf(ii).pre_cp_time, pf(ii).pre_cp_param, time(kk), time(kk+L));
        end
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
        
        %%% Clutter and likelihood %%%

        % Loop through observations
        last_cp_time = pf(ii).pre_cp_time;
        last_cp_param = pf(ii).pre_cp_time;
        if isempty(last_cp_time)
            last_cp_time = -inf;
        end
        rb_mn = pf(ii).pre_rb_mn;
        rb_vr = pf(ii).pre_rb_vr;
        clut_indic = pf(ii).pre_clut;
        for ll = 1:L
            
            % Update changepoint if we've past one
            if ~isempty(pf(ii).win_cp_time) && (time(kk+ll) > pf(ii).win_cp_time) && (last_cp_time < pf(ii).win_cp_time)
                last_cp_time = pf(ii).win_cp_time;
                last_cp_param = pf(ii).win_cp_param;
                rb_vr = rb_vr + model.w_trans_vr;
            end
            
            % Interpolation and Kalman filtering
            H = heartbeat_interpolation(algo, model, time(kk+ll), last_cp_time);
            [noclut_rb_mn, noclut_rb_vr, ~, s_mn, s_vr] = kf_update(rb_mn, rb_vr, observ(kk+ll), H, model.y_obs_vr);
            noclut_lhood = loggausspdf(observ(kk+ll), s_mn, s_vr);
            
            % Clutter sampling
            clut_rb_mn = rb_mn; clut_rb_vr = rb_vr;
            clut_lhood = loggausspdf(observ(kk+ll), 0, model.y_clut_vr);
            if (mm==1)
                clut_prior = log([0; 1]);
            elseif (clut_indic==0) && ((kk+ll)<(last_clut+algo.min_noclut))
                clut_prior = log([0; 1]);
            elseif (clut_indic==1) && ((kk+ll)<(last_noclut+algo.min_clut))
                clut_prior = log([1; 0]);
            else
                clut_prior = log(model.clut_trans(:,clut_indic+1));
            end
            clut_prob = clut_prior + [clut_lhood; noclut_lhood];
            lhood = logsumexp(clut_prob);
            clut_prob = clut_prob - lhood;
            clut_indic = log(rand)<clut_prob(1);
            if clut_indic
                rb_mn = clut_rb_mn;
                rb_vr = clut_rb_vr;
                last_clut = kk+ll;
            else
                rb_mn = noclut_rb_mn;
                rb_vr = noclut_rb_vr;
                last_noclut = kk+ll;
            end

            % Store everything
            pf(ii).win_rb_mn(:,ll) = rb_mn;
            pf(ii).win_rb_vr(:,:,ll) = rb_vr;
            pf(ii).win_obslhood(ll) = lhood;
            pf(ii).win_signal_mn(ll) = s_mn;
            pf(ii).win_signal_vr(ll) = s_vr;
            pf(ii).win_clut(ll) = clut_indic;

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
    
    assert(~any(isnan([pf.weight])));
    assert(~all(isinf([pf.weight])));
    
    % Particle smoother
    for ii = 1:Nf
        
        a_idx = pf(ii).ancestor;
        
        ps(ii) = last_ps(a_idx);
        if ps(ii).Ncp > 0
            ps(ii).cp_time(ps(ii).Ncp) = pf(ii).pre_cp_time;
        end
        if pf(ii).win_cp_time < time(kk+S);
            ps(ii).Ncp = ps(ii).Ncp + 1;
            ps(ii).cp_time(ps(ii).Ncp) = pf(ii).win_cp_time;
            ps(ii).cp_param(:,ps(ii).Ncp) = pf(ii).win_cp_param;
        end
        ps(ii).rb_mn(:,kk+1:kk+S) = pf(ii).win_rb_mn(:,1:S);
        ps(ii).clut(kk+1:kk+S) = pf(ii).win_clut(1:S);
        ps(ii).signal_mn(kk+1:kk+S) = pf(ii).win_signal_mn(1:S);
        ps(ii).signal_vr(kk+1:kk+S) = pf(ii).win_signal_vr(1:S);
        
    end
    
    % Diagnostics
    if display.plot_during && (mm>display.plot_after_frame)
        figure(display.h_pf(1)); clf; hold on; hist(diagnostic_lastest_cp_time, 100);
        figure(display.h_pf(2)); clf; hold on; hist(diagnostic_lastest_cp_param(1,:), 100);
%         figure(display.h_pf(3)); clf; hold on; hist(diagnostic_lastest_cp_param(2,:), 100);
%         figure(display.h_pf(4)); clf; hold on; hist(diagnostic_lastest_cp_param(3,:), 100);
        figure(display.h_pf(3)); clf; hold on; plot( cell2mat(arrayfun(@(x) {x.win_rb_mn(:,end)}, pf')) );
        figure(display.h_pf(4)); clf; hold on; plot( cat(1,pf.win_signal_mn)' ); plot(observ(kk+1:kk+L),'linewidth',3);
        figure(display.h_pf(5)); clf; hold on; plot( cat(2,pf.pre_rb_mn) );
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
            'clut', cell(Nf,1), ...                 Clutter over all time
            'signal_mn', cell(Nf,1), ...            Predicted signal mean over all time
            'signal_vr', cell(Nf,1));%              Predicted signal variance over all time

[ps.Ncp] = deal(0);
[ps.cp_time] = deal(zeros(1,0));
[ps.cp_param] = deal(zeros(model.dp,0));
[ps.rb_mn] = deal(zeros(model.dw,K));

end

function [pf] = init_pf(algo, model, L)
% Initialise particle smoother structure

Nf = algo.Nf;
K = model.K;

pf = struct('pre_cp_time', cell(Nf,1), ...          Most recent changepoint to occur before the window ...
            'pre_cp_param', cell(Nf,1), ...         ... and the corresponding parameters
            'ante_cp_time', cell(Nf,1), ...         Next most recent changepoint to occur before the window ...
            'ante_cp_param', cell(Nf,1), ...        ... and the corresponding parameters
            'pre_rb_mn', cell(Nf,1), ...            Mean of the Rao-Blackwellised bit before the window
            'pre_rb_vr', cell(Nf,1), ...            (Co)variance of the Rao-Blackwellised bit before the window
            'pre_clut', cell(Nf,1), ...             Clutter indicator variable for the most recent observation before the window
            'pre_last_clut', cell(Nf,1), ...        Time index of the last clutter observation
            'pre_last_noclut', cell(Nf,1), ...      Time index of the last no-clutter observation
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

[pf.pre_rb_mn] = deal(zeros(model.dw,1));
[pf.pre_rb_vr] = deal(zeros(model.dw,model.dw,1));
[pf.pre_clut] = deal(0);
[pf.pre_last_clut] = deal(-inf);
[pf.pre_last_clut] = deal(-inf);
[pf.win_rb_mn] = deal(zeros(model.dw,L));
[pf.win_rb_vr] = deal(zeros(model.dw,model.dw,L));
[pf.win_clut] = deal(zeros(1,L));
[pf.win_signal_mn] = deal(zeros(1,L));
[pf.win_signal_vr] = deal(zeros(1,L));
[pf.win_obslhood] = deal(zeros(1,L));
[pf.ancestor] = deal(0);
[pf.weight] = deal(0);

end


