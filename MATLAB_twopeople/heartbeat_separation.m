function [mu, Sigma, H] = heartbeat_separation( display, algo, model, time, observ, beat )
%HEARTBEAT_SEPARATION Separate heartbeats with a Gibbs Sampler

% Set up interpolation matrixes
[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);

% Other matrixes
Pw = model.w_prior_vr;
if model.np == 1
    Pw = blkdiag(Pw, Pw, Pw, Pw);
elseif model.np == 2
    Pw = blkdiag(Pw, Pw, Pw, Pw, Pw, Pw, Pw, Pw);
end

% Posterior
Sigma = inv( inv(Pw) + H'*H/model.y_obs_vr );
mu = Sigma*H'*Y/model.y_obs_vr;

% w_av{1} = reshape(mu(1:4*model.dw), model.dw, 4)';
% w_av{2} = reshape(mu(4*model.dw+1:end), model.dw, 4)';

end

