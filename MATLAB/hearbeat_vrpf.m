function [ pf ] = hearbeat_vrpf( display, algo, model, time, observ )
%HEARBEAT_PF Run a particle filter to infer heartbeats from a BCG or ECG
%signal.

% This is a BOOTSTRAP variable rate particle filter
%
% The nonlinear states for estimation are the beat times, the pulse and the
% beat amplitude. An additional changepoint parameters, the beat waveform
% is Rao-blackwellised.

% pf(k).cp_time  is a cell array containing the changepoint time(s) sampled
%                for this particle between time(k-1) and time(k). time(0)=0.
% pf(k).cp_param contains the corresponding changepoint parameters in row 
%                vectors. So cp_param{ii}(:,c) corresponds to
%                cp_time{ii}(c)
% pf(k).rb_mn and pf(k).rb_vr contain the mean and covariance of the
%                Gaussian density for w given observations 1:k.
% Ancestors and weights as usual particle filter

% Make local copies of useful numbers
Nf = algo.Nf;
K = model.K;

% Initialise particle filter structure array
pf = struct('Ncp', cell(K,1), ...                       Number of changepoints
            'cp_time', cell(K,1), ...           Most recent changepoint...
            'cp_param', cell(K,1), ...          ... and the corresponding parameters
            'last_cp_time', cell(K,1), ...      Changepoint before that ...
            'last_cp_param', cell(K,1), ...     ... and the corresponding parameters
            'clut_indic', cell(K,1), ...        Clutter indicator variable
            'rb_mn', cell(K,1), ...             Mean of the Rao-Blackwellised bit
            'rb_vr', cell(K,1), ...             (Co)variance of the Rao-Blackwellised bit
            'ancestor', cell(K,1), ...          Ancestor particle
            'weight', cell(K,1));
pf(1).cp_time = repmat({0}, 1, Nf);
pf(1).cp_param = repmat({zeros(model.dp,1)}, 1, Nf);
pf(1).Ncp = ones(1, Nf);
pf(1).cp_time = zeros(1, Nf);
pf(1).cp_param = zeros(model.dp, Nf);
pf(1).last_cp_time = NaN(1, Nf);
pf(1).last_cp_param = NaN(model.dp, Nf);
pf(1).clut_indic = zeros(1, Nf);
pf(1).rb_mn = zeros(model.dw, Nf);
pf(1).rb_vr = zeros(model.dw, model.dw, Nf);
pf(1).cp_rb_mn = zeros(model.dw, Nf);
pf(1).cp_rb_vr = zeros(model.dw, model.dw, Nf);
pf(1).ancestor = zeros(1, Nf);
pf(1).weight = zeros(1, Nf);

if display.text
    fprintf(1, 'Particle filter time step %u.\n', 1);
end

% Sample initial changepoint from prior
for ii = 1:Nf
    [pf(1).cp_param(:,ii), ~] = heartbeat_cpparamprior(model);
    pf(1).rb_mn(:,ii) = model.w_prior_mn;
    pf(1).rb_vr(:,:,ii) = model.w_prior_vr;
end

% Loop through time
for kk = 2:K
    
    if display.text
        fprintf(1, 'Particle filter time step %u.\n', kk);
    end
    
    % Diagnostics
    diagnostic_lastest_cp_time = zeros(1,Nf);
    diagnostic_lastest_cp_param = zeros(model.dp,Nf);
    
    % Shall we resample this frame?
%     if rem(kk,1)==0
    if kk > algo.no_resamp_period
        flag_resampling = true;
    else
        flag_resampling = false;
    end
    
    % Sample ancestors
    if flag_resampling
        sampling_weights = pf(kk-1).weight;
        sampling_weights = sampling_weights-max(sampling_weights);
        sampling_weights = exp(sampling_weights);
        sampling_weights = sampling_weights/sum(sampling_weights);
        recent_jumps = pf(kk-1).cp_time > (time(kk)-10/model.fs);
        sampling_weights(recent_jumps) = 1/Nf;
        sampling_weights(~recent_jumps) = sampling_weights(~recent_jumps)*((Nf-sum(recent_jumps))/Nf)/sum(eps+sampling_weights(~recent_jumps));
%         sampling_weights = max(sampling_weights, 1);
%         sampling_weights = sampling_weights/sum(sampling_weights);
%         low_weights = sampling_weights < (1/Nf);
%         sampling_weights(low_weights) = 1/Nf;
%         sampling_weights(~low_weights) = sampling_weights(~low_weights)*((Nf-sum(low_weights))/Nf)/sum(sampling_weights(~low_weights));
        sampling_weights = log(sampling_weights);
        pf(kk).ancestor = sample_weights(algo, sampling_weights, Nf);
    else
        pf(kk).ancestor = 1:Nf;
    end
    
    % Initialise arrays
    pf(kk).weight = zeros(1, Nf);
    pf(kk).rb_mn = zeros(model.dw, Nf);
    pf(kk).rb_vr = zeros(model.dw, model.dw, Nf);
    pf(kk).cp_time = zeros(1, Nf);
    pf(kk).cp_param = zeros(model.dp, Nf);
    pf(kk).cp_rb_mn = zeros(model.dw, Nf);
    pf(kk).cp_rb_vr = zeros(model.dw, model.dw, Nf);
    pf(kk).last_cp_time = zeros(1, Nf);
    pf(kk).last_cp_param = zeros(model.dp, Nf);
    pf(kk).clut_indic = zeros(1, Nf);
    
    % Loop through particles
    for ii = 1:Nf
        
        % Get ancestor index
        a_idx = pf(kk).ancestor(ii);
        
        % Get most recent changepoint
        latest_cp_time = pf(kk-1).cp_time(a_idx);
        latest_cp_param = pf(kk-1).cp_param(:,a_idx);
        latest_cp_rb_mn = pf(kk-1).cp_rb_mn(:,a_idx);
        latest_cp_rb_vr = pf(kk-1).cp_rb_vr(:,:,a_idx);
        
        % Get last rb estimate
        last_rb_mn = pf(kk-1).rb_mn(:,pf(kk).ancestor(ii));
        last_rb_vr = pf(kk-1).rb_vr(:,:,pf(kk).ancestor(ii));
        
        % Get clutter indicator
        last_clut_indic = pf(kk-1).clut_indic(ii);
        
        % Sample changepoint transition density
        [ cp_time, cp_param, ~ ] = heartbeat_cptransition(model, latest_cp_time, latest_cp_param, time(kk-1), time(kk));
        
        if ~isempty(cp_time)
            % Changepoint has occured!
            
            % Changepoint shift
            last_cp_time = latest_cp_time;
            last_cp_param = latest_cp_param;
            Ncp = pf(kk-1).Ncp(a_idx) + 1;
            
            % Do a Kalman filter prediction
            [rb_mn, rb_vr] = kf_predict(last_rb_mn, last_rb_vr, eye(model.dw), model.w_trans_vr);
            
            cp_rb_mn = rb_mn;
            cp_rb_vr = rb_vr;
            
        else
            % Changepoint has not occured - keep the previous estimates
            rb_mn = last_rb_mn;
            rb_vr = last_rb_vr;
            cp_time = latest_cp_time;
            cp_param = latest_cp_param;
            cp_rb_mn = latest_cp_rb_mn;
            cp_rb_vr = latest_cp_rb_vr;
            last_cp_time = pf(kk-1).last_cp_time(a_idx);
            last_cp_param = pf(kk-1).last_cp_param(:,a_idx);
            Ncp = pf(kk-1).Ncp(a_idx);
            
        end
        
        % Store the latest changepoints
        pf(kk).cp_time(ii) = cp_time;
        pf(kk).cp_param(:,ii) = cp_param;
        pf(kk).cp_rb_mn(:,ii) = cp_rb_mn;
        pf(kk).cp_rb_vr(:,:,ii) = cp_rb_vr;
        pf(kk).last_cp_time(ii) = last_cp_time;
        pf(kk).last_cp_param(:,ii) = last_cp_param;
        pf(kk).Ncp(ii) = Ncp;
        
        % Linear observation and update bit
        
        % Interpolation vector
        H = cp_param(2)*heartbeat_interpolation(algo, model, time(kk), cp_time);
        
        % Calculate clutter proposal probabilities
        [rb_mn_noclut, rb_vr_noclut, ~, ~, ~, lh_noclut] = kf_update(rb_mn, rb_vr, observ(:,kk), H, model.y_obs_vr);
        [rb_mn_clut,   rb_vr_clut,   ~, ~, ~, lh_clut  ] = kf_update(rb_mn, rb_vr, observ(:,kk), H, model.y_clut_vr);
        clut_prob = [model.clut_trans(2, last_clut_indic+1) * lh_clut; model.clut_trans(1, last_clut_indic+1)*lh_noclut];
        lh_prob = log(sum(clut_prob));
        clut_prob = clut_prob/sum(clut_prob);
        
        % Propose a value for the clutter indicator
        clut_indic = rand<clut_prob(1);
        
        % Store updated values
        if clut_indic == 0
            pf(kk).rb_mn(:,ii) = rb_mn_noclut;
            pf(kk).rb_vr(:,:,ii) = rb_vr_noclut;
            pf(kk).clut_indic(ii) = 0;
        else
            pf(kk).rb_mn(:,ii) = rb_mn_clut;
            pf(kk).rb_vr(:,:,ii) = rb_vr_clut;
            pf(kk).clut_indic(ii) = 1;
        end
        
        % Weight
%         pf(kk).weight(ii) = lh_prob;
        if flag_resampling
            pf(kk).weight(ii) = pf(kk-1).weight(a_idx) - sampling_weights(a_idx) + lh_prob;
%             pf(kk).weight(ii) = lh_prob;
        else
            pf(kk).weight(ii) = pf(kk-1).weight(a_idx) + lh_prob;
        end
        
        if flag_resampling
            % Resample-move
            [cp_time, cp_param, rb_mn, rb_vr] = heartbeat_rm(algo, model, cp_time, cp_param, cp_rb_mn, cp_rb_vr, last_cp_time, last_cp_param, time(kk), time, observ, pf(kk).rb_mn(:,ii), pf(kk).rb_vr(:,:,ii));
            pf(kk).cp_time(ii) = cp_time;
            pf(kk).cp_param(:,ii) = cp_param;
            pf(kk).rb_mn(:,ii) = rb_mn;
            pf(kk).rb_vr(:,:,ii) = rb_vr;
        end
        
        % Diagnostics
        diagnostic_lastest_cp_time(ii) = cp_time;
        diagnostic_lastest_cp_param(:,ii) = cp_param;
        
    end
    
    assert(~any(isnan(pf(kk).weight)));
    assert(~all(isinf(pf(kk).weight)));
    
    % Diagnostics
    if display.plot_during% && (kk>60)
        figure(display.h_pf(1)); clf; hold on; hist(diagnostic_lastest_cp_time, 100);
        figure(display.h_pf(2)); clf; hold on; hist(diagnostic_lastest_cp_param(1,:), 100);
        figure(display.h_pf(3)); clf; hold on; hist(diagnostic_lastest_cp_param(2,:), 100);
        figure(display.h_pf(4)); clf; hold on; plot(pf(kk).rb_mn);
        pause
    end
    
    % Clear last particles (to save memory)
    pf(kk-1).rb_vr = [];
    pf(kk-1).cp_rb_vr = [];
    
end

end

