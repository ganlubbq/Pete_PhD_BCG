function [ beat ] = heartbeat_initialise( algo, model, time, observ, start_time, end_time )
%HEARTBEAT_INITIALISE Draw a sample from the posterior for the first chunk
%of data, using a simplified model

M = 100;

% Fixed things
wf_mn = model.w_prior_mn;
wf_vr = model.w_prior_vr;

% Start by sampling the prior
prior_prob = zeros(model.np,3);
trans_prob = [0;0];
for pp = 1:model.np
    [pre_time, pre_param, period, prior_prob(pp,:)] = sample_start_time(model);
    [beat(pp), trans_prob(pp)] = heartbeat_beattrans(model, pre_time, pre_param, [], start_time, end_time, []);
end

% Evaluate likelihood
[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);
R = H*wf_vr*H'+model.y_obs_vr*eye(length(Y));
lhood = loggausspdf(Y, H*wf_mn, R);

% M-H loop
for mm = 1:M
    
    ppsl_beat = beat;
    ppsl_prior_prob = prior_prob;
    ppsl_trans_prob = trans_prob;
    
    % Which person?
    pp = rem(mm,2)+1;
    
    % Propose a change
    [ ppsl_pre_param, ppsl_prior_prob(1) ] = heartbeat_paramtrans(model, [], []);
    [ ppsl_period, ppsl_prior_prob(2) ] = heartbeat_periodtrans(model, ppsl_pre_param, []);
    if beat(pp).pre_time < -ppsl_period
        ppsl_prior_prob(3) = -Inf;
    else
        ppsl_prior_prob(3) = log(1/period);
    end
    [ ppsl_beat(pp), ppsl_trans_prob(pp) ] = heartbeat_beattrans(model, beat(pp).pre_time, ppsl_pre_param, [], start_time, end_time, []);
    
    % Evaluate likelihood
    [ppsl_H, ~] = heartbeat_obsmat(algo, model, time, observ, ppsl_beat);
    ppsl_R = ppsl_H*wf_vr*ppsl_H'+model.y_obs_vr*eye(length(Y));
    ppsl_lhood = loggausspdf(Y, ppsl_H*wf_mn, ppsl_R);
    
    if log(rand) < ( (ppsl_lhood+ppsl_trans_prob(pp)+sum(prior_prob))-(lhood+trans_prob(pp)+sum(ppsl_prior_prob)) )
        beat = ppsl_beat;
        lhood = ppsl_lhood;
        trans_prob = ppsl_trans_prob;
        prior_prob = ppsl_prior_prob;
        H = ppsl_H;
        R = ppsl_R;
    end
    
end

end



function [pre_time, pre_param, period, prob] = sample_start_time(model)

% Sample parameter from prior
[ pre_param, param_prob ] = heartbeat_paramtrans(model, [], []);

% Sample a beat period from the prior
[ period, period_prob ] = heartbeat_periodtrans(model, pre_param, []);

% Sample beat time
pre_time = unifrnd(0-period, 0);
time_prob = log(1/period);

% Combine densities
prob = [param_prob, period_prob, time_prob];

end