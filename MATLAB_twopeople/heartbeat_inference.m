function [ ps, ess ] = heartbeat_inference( display, algo, model, time, observ )
%HEARTBEAT_INFERENCE Separate BCG heartbeat signals. Joint estimation of
%beat timings and waveforms, using a variable rate particle filter on a
%rolling window.

Nf = algo.Nf;

% Initialise particle structures and arrays
ps = ps_init(model, Nf);
last_pf = pf_init(model, Nf);
ess = zeros(ceil(model.K/algo.S),1);

pqratio = zeros(model.np,1);

if display.text
    
    fprintf(1, 'Initialisation Loop:\n');
    tic;
end
% Particle loop
for ii = 1:Nf
    
    % How long is the initial stretch? (MAKE THIS BIG, OR 0)
    Lstart = 0;algo.L;%algo.L-algo.S;
    
    if Lstart > 0
        % Sample initial values for the parameters and changepoints
        start_time = time(1);
        end_time = time(Lstart);
        beat = heartbeat_initialise( algo, model, time(1:Lstart), observ(:,1:Lstart), start_time, end_time );
        
        % Propose sequence
        [beat] = heartbeat_sequenceproposal(algo, model, beat, time(1:Lstart), observ(:,1:Lstart), model.w_prior_mn, model.w_prior_vr);
        
    else
        beat = heartbeat_initialise( algo, model, time(1:Lstart), observ(:,1:Lstart), 0, 0 );
        
    end
        
    % Store
    last_pf(ii).beat = beat;
    last_pf(ii).sf_ratio = zeros(model.np,1);
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
    
    if mm == 3
        algo.Nf = algo.Nf/algo.drop_fact;
        Nf = algo.Nf;
    end
    
    % Time, step size and window length
    mm = mm + 1;
    kk = kk + algo.S;
    if kk >= model.K
        break;
    end
    
    if mm == 1
        start_time = time(kk+1);
        last_L = algo.S+Lstart;
    else
        start_time = time(kk);
        last_L = algo.L;
    end
    
    L = min(algo.L, model.K-kk);
    S = min(algo.S, model.K-kk);
    end_time = time(kk+L);
        
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
        old_sf_ratio = last_pf(anc).sf_ratio;

        % Heartbeat sequence proposal
        [beat, sf_ratio] = heartbeat_sequenceproposal(algo, model, beat, time(kk+1:kk+L), observ(:,kk+1:kk+L), wf_mn, wf_vr);
        
        % Store it
        pt.beat = beat;
        pt.sf_ratio = sf_ratio;
        for pp = 1:model.np
            % Update smoothed particles too
            mrc_idx = most_recent_changepoint(ps(ii).beat(pp).time, start_time);
            ps(ii).beat(pp).time(mrc_idx+1:end) = [];
            ps(ii).beat(pp).param(:,mrc_idx+1:end) = [];
            ps(ii).beat(pp).param(:,mrc_idx) = beat(pp).pre_param;
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
                    + sum(pqratio) - sum(old_pqratio) ...
                    + sum(sf_ratio) - sum(old_sf_ratio);
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
