function [ beat, period, prob ] = heartbeat_beatprior( model )
%HEARTBEAT_BEATPRIOR Sample a preceeding beat time and parameter.

% Sample parameter from prior
[ param, param_prob ] = heartbeat_paramtrans(model, [], []);

% Sample a beat period from the prior
[ period, period_prob ] = heartbeat_periodtrans(model, param, []);

% Sample beat time
tau0 = unifrnd(0-period, 0);
time_prob = log(1/period);

% Make a beat
beat = beat_init(model, tau0, param, [], [], []);

% Combine densities
prob = param_prob + period_prob + time_prob;

end

