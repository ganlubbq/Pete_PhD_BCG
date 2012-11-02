function [ cp_time, cp_param ] = heartbeat_rm( algo, model, cp_time, cp_param, last_cp_time, last_cp_param )
%HEARTBEAT_RM Do some resample move on the latest changepoint

% Propose a new value for the beat period from the transition density
% (Gibbs, so no acceptance probability needed.

if any(isnan(last_cp_param))
    cp_param(1) = 0;
    while cp_param(1) <= model.p_min
        cp_param(1) = gamrnd(model.p_prior_shape, model.p_prior_scale);
    end
else
    cp_param(1) = 0;
    while cp_param(1) <= model.p_min
        cp_param(1) = mvnrnd(last_cp_param(1), model.p_trans_vr);
    end
end

end

