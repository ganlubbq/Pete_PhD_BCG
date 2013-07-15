function [ beat, prob ] = heartbeat_beatprior( model, end_time )
%HEARTBEAT_BEATPRIOR Sample a preceeding beat time and parameter.

% Sample parameter from prior
[ param, param_prob ] = heartbeat_paramtrans(model, [], []);

% Sample a beat period from the prior
[ period, period_prob ] = heartbeat_periodtrans(model, param, []);

% Sample beat time
tau0 = unifrnd(end_time-period, end_time);
time_prob = log(1/period);

% Build a beat
beat = beat_init(model, tau0, param, [], [], []);

% Combine densities
prob = param_prob + period_prob + time_prob;

end

