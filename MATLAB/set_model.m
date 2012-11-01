% Set model parameters

% Basics
model.K = 800;              % Number of observations
model.fs = 30;              % Sampling frequency of observations (after load_and_calibrate, which downsamples)
model.dp = 2;               % Number of changepoint parameter dimensions (beat amplitude and period)
model.dw = 30;              % Number of 

% Priors
load('template_beat.mat');
model.w_prior_mn = template;
model.w_prior_vr = 0.1*eye(model.dw, model.dw);
model.p_prior_shape = 200;
model.p_prior_scale = 0.006;
model.a_prior_scale = 1;

% Transition models
model.tau_trans_scale = 1/1000;         % Changepoint time transition density (gamma) scale (shape is given by p/scale)
model.p_trans_vr = 0.01^2;              % Beat period transition density (normal) variance (mean is the previous value)
model.a_trans_vr = 0.05^2;              % Beat amplitude transition density (normal) variance (mean is the previous value)
model.w_trans_vr = 0.001*eye(model.dw); % Waveform transition density (normal) covariance matrix (mean is the previous value)

% Observation model
model.y_obs_vr = 0.1^2;                 % Observation variance
