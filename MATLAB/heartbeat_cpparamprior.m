function [ cp_param, prob ] = heartbeat_cpparamprior( model, cp_param )
%HEARTBEAT_CPPARAMPRIOR Sample and/or calculate the probability for the 
% changepoint parameter prior.

% prob is a log-probability.

% Sample state if not provided
if (nargin<2)||isempty(cp_param)
    cp_param = zeros(model.dp, 1);
    cp_param(1) = gamrnd(model.p_prior_shape, model.p_prior_scale);
    cp_param(2) = mvnrnd(model.a_prior_mn, model.a_prior_vr);
    cp_param(3) = exprnd(model.b_prior_mn);
end

% Calculate probability if required
if nargout>1
    prob = log(gampdf(cp_param(1), model.p_prior_shape, model.p_prior_scale)) ...
         + log(mvnpdf(cp_param(2), model.a_prior_mn, model.a_prior_vr)) ...
         + log(exppdf(cp_param(3), model.b_prior_mn));
else
    prob = [];
end

end

