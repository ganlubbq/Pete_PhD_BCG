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

% Initialise arrays
pf = struct('cp_time', cell(K,1), ...
            'cp_param', cell(K,1), ...
            'rb_mn', cell(K,1), ...
            'rb_vr', cell(K,1), ...
            'ancestor', cell(K,1), ...
            'weight', cell(K,1));
pf(1).cp_time = repmat({0}, 1, Nf);
pf(1).cp_param = repmat({zeros(model.dp,1)}, 1, Nf);
pf(1).rb_mn = zeros(model.dw, Nf);
pf(1).rb_vr = zeros(model.dw, model.dw, Nf);
pf(1).ancestor = zeros(1, Nf);
pf(1).weight = zeros(1, Nf);

if display.text
    fprintf(1, 'Particle filter time step %u.\n', 1);
end

% Sample initial changepoint from prior
for ii = 1:Nf
    [pf(1).cp_param{ii}, ~] = heartbeat_cpparamprior(model);
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
    
    % Sample ancestors
    pf(kk).ancestor = sample_weights(algo, pf(kk-1).weight, Nf);
    
    % Initialise arrays
    pf(kk).weight = zeros(1, Nf);
    pf(kk).rb_mn = zeros(model.dw, Nf);
    pf(kk).rb_vr = zeros(model.dw, model.dw, Nf);
    pf(kk).cp_time = cell(1, Nf);
    pf(kk).cp_param = cell(1, Nf);
    
    % Loop through particles
    for ii = 1:Nf
        
        % Get last changepoint
        last_cp_time = [];
        kt = kk;
        anc = ii;
        while isempty(last_cp_time)
            anc = pf(kt).ancestor(anc);
            kt = kt - 1;
            last_cp_time = pf(kt).cp_time{anc};
        end
        last_cp_time = last_cp_time(end);
        last_cp_param = pf(kt).cp_param{anc}(:,end);
        
        % Get last rb estimate
        last_rb_mn = pf(kk-1).rb_mn(:,pf(kk).ancestor(ii));
        last_rb_vr = pf(kk-1).rb_vr(:,:,pf(kk).ancestor(ii));
        
        % Sample changepoint transition density
        [ cp_time, cp_param, ~ ] = heartbeat_cptransition(model, last_cp_time, last_cp_param, time(kk-1), time(kk));
        pf(kk).cp_time{ii} = cp_time;
        pf(kk).cp_param{ii} = cp_param;
        
        % Linear update bit
        
        if ~isempty(cp_time)
            % Changepoint has occured - do a Kalman filter prediction
            [rb_mn, rb_vr] = kf_predict(last_rb_mn, last_rb_vr, eye(model.dw), model.w_trans_vr);
        else
            % Changepoint has not occured - keep the previous estimate
            rb_mn = last_rb_mn;
            rb_vr = last_rb_vr;
            cp_time = last_cp_time;
            cp_param = last_cp_param;
        end
        
        % Kalman filter update
        H = cp_param(2)*heartbeat_interpolation(algo, model, time(kk), cp_time);
        [rb_mn, rb_vr, ~, ~, ~, lh] = kf_update(rb_mn, rb_vr, observ(:,kk), H, model.y_obs_vr);
        lh_prob = log(lh);
        pf(kk).rb_mn(:,ii) = rb_mn;
        pf(kk).rb_vr(:,:,ii) = rb_vr;
        
        % Weight
        pf(kk).weight(ii) = lh_prob;
        
        % Diagnostics
        diagnostic_lastest_cp_time(ii) = cp_time;
        diagnostic_lastest_cp_param(:,ii) = cp_param;
        
    end
    
    % Diagnostics
    if display.plot_during
        figure(display.h_pf(1)); clf; hold on; hist(diagnostic_lastest_cp_time);
        figure(display.h_pf(2)); clf; hold on; hist(diagnostic_lastest_cp_param(1,:));
        figure(display.h_pf(3)); clf; hold on; hist(diagnostic_lastest_cp_param(2,:));
        figure(display.h_pf(4)); clf; hold on; plot(pf(kk).rb_mn);
        pause
    end
    
end

end

