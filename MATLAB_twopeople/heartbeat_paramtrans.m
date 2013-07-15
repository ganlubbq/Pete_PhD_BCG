function [ param, prob ] = heartbeat_paramtrans( model, last_param, param )
%HEARTBEAT_PARAMTRANS Sample and evaluate transition density for heartbeat
%parameters.

% Set hyperparameters depending on whether this is the first frame or not
if isempty(last_param)
    p_shape = model.p_prior_shape;
    p_scale = model.p_prior_scale;
else
    p_scale = model.p_trans_scale;
    p_shape = last_param/p_scale;
end

if isempty(param)
    param = gamrnd(p_shape, p_scale);
end

if nargout > 1
    prob = log(gampdf(param, p_shape, p_scale));
end

end

