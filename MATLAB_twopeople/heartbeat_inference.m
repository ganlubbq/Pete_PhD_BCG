function [ ps, ess ] = heartbeat_inference( display, algo, model, time, observ )
%HEARTBEAT_INFERENCE Separate BCG heartbeat signals. Joint estimation of
%beat timings and waveforms, using a variable rate particle filter on a
%rolling window.

Nf = algo.Nf;

% Initialise particle smoother structure
ps = ps_init(model, Nf);

% Initialise a blank particle filter structure
last_pf = pf_init(model, Nf);

pqratio = zeros(model.np,1);

if display.text
    
    fprintf(1, '\n\nInitialisation Loop:\n');
    tic;
end
% Particle loop
for ii = 1:Nf
    
    % How long is the initial stretch?
    Lstart = algo.L-algo.S;
    
    % Sample initial values for the parameters and changepoints
    start_time = time(1); end_time = time(Lstart);
    beat = heartbeat_initialise( algo, model, time(1:Lstart), observ(:,1:Lstart), start_time, end_time );
    
    % Propose sequence
    [beat] = heartbeat_sequenceproposal(algo, model, beat, time(1:Lstart), observ(:,1:Lstart), model.w_prior_mn, model.w_prior_vr);
    
    % Store
    last_pf(ii).beat = beat;
    ps(ii).beat = beat;
    
    % Probabilities
    [H, Y] = heartbeat_obsmat(algo, model, time(1:Lstart), observ(:,1:Lstart), ps(ii).beat);
    lhood = loggausspdf(Y, H*model.w_prior_mn, H*model.w_prior_vr*H'+model.y_obs_vr*eye(length(Y)));
    for pp = 1:model.np
        pqratio(pp) = sum(beat(pp).pqratio);
    end
    
    % Update weight
    last_pf(ii).weight = lhood + sum(pqratio);
    
end
if display.text
    % Text output
    fprintf(1, '     - Effective sample size: %f.\n', calc_ESS([last_pf.weight]));
    fprintf(1, '     - That took %f seconds.\n', toc);
end

% Particle filter loop
kk = -algo.S;
mm = 0;
while 1
    
    % Time, step size and window length
    mm = mm + 1;
    kk = kk + algo.S;
    if kk >= model.K
        break;
    end
    
    if mm == 1
        start_time = time(kk+1);
    else
        start_time = time(kk);
    end
    
    L = min(algo.L, model.K-kk);
    S = min(algo.S, model.K-kk);
    end_time = time(kk+L);
    last_L = algo.L;
    
%     if mm == 3
%         Nf = Nf/10;
%     end
        
    if display.text
        % Text output
        fprintf(1, 'Particle filter iteration %u:\n', mm);
        fprintf(1, '     - Window from %f to %f (indexes %u to %u of %u)\n',  start_time, end_time, kk+1, kk+L, model.K);
        tic;
    end
    
    % Initialise particle filter structure
    pf = pf_init(model, Nf);
    pqratio = zeros(model.np,1);
    old_pqratio = zeros(model.np,1);
    
    % Particle selection
    ancestor = sample_weights([last_pf.weight], Nf, algo.resam_type);
    ps = ps(ancestor);
    
    % Particle loop
    for ii = 1:Nf
        
        % Which ancestor?
        anc = ancestor(ii);
        
        % Initialise particle
        pt = pf_forwardparticle(model, anc, last_pf(anc), start_time, end_time);

        % Get the beat
        beat = pt.beat;
        
        % Find parameters of Gaussian posterior waveform distribution
        [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time(1:kk), observ(:,1:kk), ps(ii).beat );
        
        % Old probabilities
        [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:min(kk+last_L-S,model.K)), observ(:,kk+1:min(kk+last_L-S,model.K)), beat);
        old_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        for pp = 1:model.np
            old_pqratio(pp) = sum(beat(pp).pqratio);
        end
        
%         % MH to refresh pre_param
%         for pp = 1:model.np
%             
%             pre_time = beat(pp).pre_time;
%             pre_param = beat(pp).pre_param;
%             ante_param = beat(pp).ante_param;
%             
%             % Propose a new pre_param
%             ppsl_pre_param = heartbeat_paramtrans(model, ante_param, []);
%             
%             % Find new and old target probabilities
%             if isempty(beat(pp).time)
%                 % Its the last beat - use the survival probability
%                 beat_period_offset = start_time - (pre_time + pre_param);
%                 ppsl_beat_period_offset = start_time - (pre_time + ppsl_pre_param);
%                 param_prob = log(1 - gamcdf( beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
%                 ppsl_param_prob = log(1 - gamcdf( ppsl_beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
%             else
%                 % Its not - use the density of the next changepoint
%                 period_offset = beat(pp).time(1) - (pre_time + pre_param);
%                 ppsl_period_offset = beat(pp).time(1) - (pre_time + ppsl_pre_param);
%                 param_prob = log(gampdf(period_offset, model.tau_trans_shape, model.tau_trans_scale))...
%                     + log(gampdf(pre_param, ante_param/model.p_trans_scale, model.p_trans_scale));
%                 ppsl_param_prob = log(gampdf(ppsl_period_offset, model.tau_trans_shape, model.tau_trans_scale))...
%                     + log(gampdf(ppsl_pre_param, ante_param/model.p_trans_scale, model.p_trans_scale));
%             end
%             
%             if log(rand) < (ppsl_param_prob-param_prob)
%                 beat(pp).pre_param = ppsl_pre_param;
%                 
%                 % Update smoothed particle too
%                 idx = find( ps(ii).beat(pp).time==beat(pp).pre_time );
%                 if ~isempty(idx)
%                     ps(ii).beat(pp).param(:,idx) = ppsl_pre_param;
%                 elseif ps(ii).beat(pp).pre_time == beat(pp).pre_time
%                     ps(ii).beat(pp).pre_param = ppsl_pre_param;
%                 else
%                     error('Where''d it go?')
%                 end
%             end
%             
%         end

        % Heartbeat sequence proposal
        [beat] = heartbeat_sequenceproposal(algo, model, beat, time(kk+1:kk+L), observ(:,kk+1:kk+L), wf_mn, wf_vr);
        
        % Store it
        pt.beat = beat;
        for pp = 1:model.np
            % Update smoothed particles too
            idx = find( ps(ii).beat(pp).time>start_time );
            ps(ii).beat(pp).time(idx) = [];
            ps(ii).beat(pp).param(:,idx) = [];
            ps(ii).beat(pp).time = [ps(ii).beat(pp).time, beat(pp).time];
            ps(ii).beat(pp).param = [ps(ii).beat(pp).param, beat(pp).param];
        end

        % Probabilities
        [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), beat);
        new_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        
        % P/Q ratio
        for pp = 1:model.np
            pqratio(pp) = sum(beat(pp).pqratio);
        end
        
        % Update weight
        pt.weight =   new_lhood    - old_lhood ...
                    + sum(pqratio) - sum(old_pqratio);
        if isinf(pt.weight)||isnan(pt.weight)
            pt.weight = -inf;
        end
        
        % Store particle
        pf(ii) = pt;
        
        % Plot
        if display.plot_during
            x_vec = H*wf_mn;
            x = reshape(x_vec,4,L);
            figure, plot(time(kk+1:kk+L), x');
        end
        
    end
    
    % Effective sample size
    ess(mm) = calc_ESS([pf.weight]);
    
    if display.text
        % Text output
        fprintf(1, '     - Effective sample size: %f.\n', ess(mm));
        fprintf(1, '     - That took %f seconds.\n', toc);
    end
    
    % Keep things that we'll need in the next interation
    last_pf = pf;
    
end

% Final smoother particle selection step
[ancestor] = sample_weights([last_pf.weight], Nf, algo.resam_type);
ps = ps(ancestor);

end
