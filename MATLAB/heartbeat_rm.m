function [ cp_time, cp_param ] = heartbeat_rm( algo, model, cp_time, cp_param, last_cp_time, last_cp_param, current_time )
%HEARTBEAT_RM Do some resample move on the latest changepoint

% Propose a new value for the beat period from the transition density
% (Gibbs, so no acceptance probability needed.

if any(isnan(last_cp_param))
    ppsl_p = gamrnd(model.p_prior_shape, model.p_prior_scale);
else
    ppsl_p = gamrnd(last_cp_param(1)/model.p_trans_scale, model.p_trans_scale);
end

new_surv_prob = log(1-gamcdf(current_time-cp_time, ppsl_p/model.tau_trans_scale, model.tau_trans_scale));
old_surv_prob = log(1-gamcdf(current_time-cp_time, cp_param(1)/model.tau_trans_scale, model.tau_trans_scale));

if log(rand) < new_surv_prob-old_surv_prob
    cp_param(1) = ppsl_p;
end

end

