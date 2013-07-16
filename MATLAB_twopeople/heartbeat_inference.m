function [ ps ] = heartbeat_inference( display, algo, model, time, observ )
%HEARTBEAT_INFERENCE Separate BCG heartbeat signals. Joint estimation of
%beat timings and waveforms, using a variable rate particle filter on a
%rolling window.

% Make local copies of useful numbers
K = model.K;        % Number of time steps
Nf = algo.Nf;       % Number of filter particles
S = algo.S;         % Time increment (number of time steps to move forward each processing frame)
L = algo.L;         % Window length

% Number of particle filter steps we need to run
M = ceil(K/S);

% Initialise particle smoother structure
ps = ps_init(model, Nf);

% Initialise a blank particle filter structure
last_pf = pf_init(model, Nf);

% Particle loop
for ii = 1:Nf
    % People loop
    for pp = 1:model.np
        
        % Sample initial values for the parameters and changepoints from the prior
        last_pf(ii).beat(pp) = heartbeat_beatprior(model, 0);
        ps(ii).beat(pp).pre_time = last_pf(ii).beat(pp).pre_time;
        ps(ii).beat(pp).pre_param = last_pf(ii).beat(pp).pre_param;
        
    end
end

% Particle filter loop
for mm = 1:M
    
    % Time, step size and window length
    if mm == 1, kk = 0; else kk = kk + S; end
    L = min(L, K-kk);
    S = min(S, K-kk);
    start_time = time(kk+1);
    end_time = time(kk+L);
    
    if display.text
        % Text output
        fprintf(1, 'Particle filter iteration %u of %u:\n', mm, M);
        fprintf(1, '     - Window from %f to %f (indexes %u to %u)\n',  start_time, end_time, kk+1, kk+L);
        tic;
    end
    
    % Initialise particle filter structure
    pf = pf_init(model, Nf);
    
    % Particle selection
    if mm > 2
        auxiliary_weights = [last_pf.weight];
    else
        auxiliary_weights = zeros(1,Nf);
    end
    [ancestor, selected_weights] = sample_weights(auxiliary_weights, Nf, algo.resam_type);
    ps = ps(ancestor);
    
    % Particle loop
    for ii = 1:Nf
        
        % Which ancestor?
        anc = ancestor(ii);
        
        % Initialise particle
        pt = pf_forwardparticle(model, anc, last_pf(anc), start_time, end_time);
        
%         % Propose new parameters for the preceeding heartbeats from a MH kernel
%         for pp = 1:model.np
%             
%             beat = pt.beat(pp);
%             ante_param = beat.ante_param;
%             
%             % Sample a new value for pre_param
%             new_pre_param = heartbeat_paramtrans(model, ante_param, []);
%             new_beat = beat;
%             new_beat.pre_param = new_pre_param;
%             
%             % Probabilities
%             new_pre_param_prob = heartbeat_preparamchange(model, new_beat, start_time);
%             old_pre_param_prob = heartbeat_preparamchange(model, beat,     start_time);
% 
%             % Accept?
%             if log(rand) < new_pre_param_prob-old_pre_param_prob
%                 pt.beat(pp) = new_beat;
%                 
%                 % Update smoothed particle too
%                 idx = find( ps(ii).beat(pp).time==pt.beat(pp).pre_time );
%                 if ~isempty(idx)
%                     ps(ii).beat(pp).param(:,idx) = new_pre_param;
%                 end
%             end
%             
%         end
        
        % Find parameters of Gaussian posterior waveform distribution
        [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time(1:kk), observ(:,1:kk), ps(ii).beat );
        
        % Old probabilities
        if mm > 1
            [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L-S), observ(:,kk+1:kk+L-S), ps(ii).beat);
            old_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        else
            old_lhood = 0;
        end
        
        % Propose new heartbeats for the window
        for pp = 1:model.np
            
            pre_time = pt.beat(pp).pre_time;
            pre_param = pt.beat(pp).pre_param;
            ante_param = pt.beat(pp).ante_param;
            pt.beat(pp) = heartbeat_beattrans(model, pre_time, pre_param, ante_param, start_time, end_time, []);
            
            % Update smoothed particles too
            idx = find( ps(ii).beat(pp).time>start_time );
            ps(ii).beat(pp).time(idx) = [];
            ps(ii).beat(pp).param(:,idx) = [];
            ps(ii).beat(pp).time = [ps(ii).beat(pp).time, pt.beat(pp).time];
            ps(ii).beat(pp).param = [ps(ii).beat(pp).param, pt.beat(pp).param];
            
        end
        
        % Probabilities
        [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), ps(ii).beat);
        new_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        
        % Update weight
        pt.weight = pt.weight - selected_weights(ii) + new_lhood - old_lhood;
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

    if display.text
        % Text output
        fprintf(1, '     - Effective sample size: %f.\n', calc_ESS([pf.weight]));
        fprintf(1, '     - That took %f seconds.\n', toc);
    end
    
    % Keep things that we'll need in the next interation
    last_pf = pf;
    
end

% Final smoother particle selection step

end

