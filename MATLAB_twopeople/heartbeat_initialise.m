function [ beat ] = heartbeat_initialise( algo, model, time, observ, start_time, end_time )
%HEARTBEAT_INITIALISE Draw a sample from the posterior for the first chunk
%of data, using a simplified model

M = 10;

% Fixed things
wf_mn = model.w_prior_mn;
wf_vr = model.w_prior_vr;

% Sampling start times from the prior
for pp = 1:model.np
    pre_time = sample_start_time(model);
    pre_param = heartbeat_paramtrans(model, [], []);
    beat(pp) = heartbeat_beattrans(model, pre_time, pre_param, [], start_time, end_time, []);
end

% Evaluate likelihood
[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);
R = H*wf_vr*H'+model.y_obs_vr*eye(length(Y));
lhood = loggausspdf(Y, H*wf_mn, R);

% M-H loop
for mm = 1:M
    
    ppsl_beat = beat;
    
    % Which person?
    if model.np == 2
        pp = rem(mm,2)+1;
    else
        pp = 1;
    end
    
    % Propose a change
    ppsl_pre_param = heartbeat_paramtrans(model, [], []);
    ppsl_beat(pp) = heartbeat_beattrans(model, beat(pp).pre_time, ppsl_pre_param, [], start_time, end_time, []);
    
    % Evaluate likelihood
    [ppsl_H, ~] = heartbeat_obsmat(algo, model, time, observ, ppsl_beat);
    ppsl_R = ppsl_H*wf_vr*ppsl_H'+model.y_obs_vr*eye(length(Y));
    ppsl_lhood = loggausspdf(Y, ppsl_H*wf_mn, ppsl_R);
    
    if log(rand) < ppsl_lhood-lhood
        beat = ppsl_beat;
        lhood = ppsl_lhood;
        H = ppsl_H;
        R = ppsl_R;
    end
    
end

end



function [pre_time] = sample_start_time(model)

% Sample parameter from prior
[ param, ~ ] = heartbeat_paramtrans(model, [], []);

% Sample a beat period from the prior
[ period, ~ ] = heartbeat_periodtrans(model, param, []);

% Sample beat time
pre_time = unifrnd(0-period, 0);

end