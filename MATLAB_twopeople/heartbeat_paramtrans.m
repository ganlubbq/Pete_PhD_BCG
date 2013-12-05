function [ param, prob ] = heartbeat_paramtrans( model, last_param, param )
%HEARTBEAT_PARAMTRANS Sample and evaluate transition density for heartbeat
%parameters.

% Set hyperparameters depending on whether this is the first frame or not
if isempty(last_param)
    p_shift = model.p_prior_shift;
    p_shape = model.p_prior_shape;
    p_scale = model.p_prior_scale;
else
    p_shift = 0;
    p_scale = model.p_trans_scale;
    p_shape = last_param/p_scale;
end

if isempty(param)
    param = p_shift + gamrnd(p_shape, p_scale);
end

if nargout > 1
    prob = log(gampdf(param-p_shift, p_shape, p_scale));
end

end

