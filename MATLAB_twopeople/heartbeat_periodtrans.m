function [ period, prob ] = heartbeat_periodtrans( model, param, period )
%HEARTBEAT_PERIODTRANS Sample and evaluate transition density for heartbeat
%period.

% Set hyperparameters
tau_shape = model.tau_trans_shape;
tau_scale = model.tau_trans_scale;
tau_shift = param(1);

if isempty(period)
    period = tau_shift + gamrnd(tau_shape, tau_scale);
end

if nargout > 1
    prob = log(gampdf(period-tau_shift, tau_shape, tau_scale));
end

end

