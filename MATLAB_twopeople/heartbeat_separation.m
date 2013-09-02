function [mu, Sigma, H, Y] = heartbeat_separation( display, algo, model, time, observ, beat )
%HEARTBEAT_SEPARATION Separate heartbeats with a Gibbs Sampler

% Set up interpolation matrixes
[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);

% Other matrixes
Pw = model.w_prior_vr;
mw = model.w_prior_mn;

% Posterior
Sigma = inv( inv(Pw) + H'*H/model.y_obs_vr );
mu = Sigma*(H'*Y/model.y_obs_vr + Pw\mw);

% w_av{1} = reshape(mu(1:4*model.dw), model.dw, 4)';
% w_av{2} = reshape(mu(4*model.dw+1:end), model.dw, 4)';

end

