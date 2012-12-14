function [ cp_param, prob ] = heartbeat_cpparamprior( model, cp_param )
%HEARTBEAT_CPPARAMPRIOR Sample and/or calculate the probability for the 
% changepoint parameter prior.

% prob is a log-probability.

%%% ADD DISTURBANCE SAMPLING

% Sample state if not provided
if (nargin<2)||isempty(cp_param)
    cp_param = zeros(model.dp, 1);
    cp_param(1) = gamrnd(model.p_prior_shape, model.p_prior_scale);
    cp_param(2) = 0;
end

% Calculate probability if required
if nargout>1
    prob = loggampdf(cp_param(1), model.p_prior_shape, model.p_prior_scale);
else
    prob = [];
end

end

