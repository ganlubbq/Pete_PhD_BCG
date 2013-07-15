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
        
    end
end

% Particle filter loop
for mm = 1:M
    
    % Time, step size and window length
    if mm == 1, kk = 1; else kk = kk + S; end
    L = min(L, K-kk);
    S = min(S, K-kk);
    start_time = time(kk);
    end_time = time(kk+L);
    
    if display.text
        % Text output
        fprintf(1, 'Particle filter iteration %u of %u:\n', mm, M);
        fprintf(1, '     - Window from %f to %f (indexes %u to %u)\n',  start_time, end_time, kk, kk+L);
        tic;
    end
    
    % Initialise particle filter structure
    pf = pf_init(model, Nf);
    
    % Particle selection
    if mm > 1
        auxiliary_weights = [last_pf.weight];
    else
        auxiliary_weights = zeros(1,Nf);
    end
    ancestor = sample_weights(auxiliary_weights, Nf, algo.resam_type);
    ps = ps(ancestor);
    
    % Particle loop
    for ii = 1:Nf
        
        % Which ancestor?
        anc = ancestor(ii);
        
        % Initialise particle
        pt = pf_forwardparticle(model, anc, last_pf(anc), start_time, end_time);
        
        % Propose new parameters for the preceeding heartbeats from a MH kernel
        for pp = 1:model.np
            
            beat = pt.beat(pp);
            ante_param = beat.ante_param;
            
            % Sample a new value for pre_param
            new_pre_param = heartbeat_paramtrans(model, ante_param, []);
            new_beat = beat;
            beat.pre_param = new_pre_param;
            
            % Probabilities
            new_pre_param_prob = heartbeat_preparamchange(model, new_beat, end_time);
            old_pre_param_prob = heartbeat_preparamchange(model, beat,     end_time);

            % Accept?
            if log(rand) < new_pre_param_prob-old_pre_param_prob
                pt.beat(pp) = new_beat;
                
                % Update smoothed particle too
                idx = find( ps(ii).beat(pp).time==pt.beat(pp).pre_time );
                if ~isempty(idx)
                    ps(ii).beat(pp).param(:,idx) = new_pre_param;
                end
            end
            
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
        
        % Find parameters of Gaussian posterior waveform distribution
        [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time(1:kk+L), observ(:,1:kk+L), ps(ii).beat );
        
        % Probabilities
        
        % Update weight
        
        % Store particle
        pf(ii) = pt;
        
    end
    
    % Normalise weights
    
    % Store in smoother structure
    
    if display.text
        % Text output
        fprintf(1, '     - That took %f seconds.\n', toc);
        tic;
    end
    
    % Keep things that we'll need in the next interation
    last_pf = pf;
    
end

% Final smoother particle selection step

end

