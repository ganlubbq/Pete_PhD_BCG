% Set model parameters

% Basics
model.K = 6500;              % Number of observations
model.fs = 30;              % Sampling frequency of observations (after load_and_calibrate, which downsamples)
model.dp = 3;               % Number of changepoint parameter dimensions (beat amplitude and period)
model.dw = 30;              % Number of 

% Priors
% load('template_beat.mat');
% model.w_prior_mn = template;
% model.w_prior_vr = 0.1*eye(model.dw, model.dw);
model.w_prior_mn = zeros(model.dw,1);
model.w_prior_vr = 10*eye(model.dw, model.dw);
model.p_prior_shape = 40;
model.p_prior_scale = 0.03;
model.a_prior_mn = 1;
model.a_prior_vr = 1E-8;
model.b_prior_mn = 0.2;

% Transition models
model.tau_trans_scale = 1E-4;           % Changepoint time transition density (gamma) scale (shape is given by p/scale)
model.p_trans_scale = 1E-3;             % Beat period transition density (gamms) scale (shape is the previous value/scale)
% model.p_trans_vr = 0.01^2;              % Beat period transition density (normal) variance (mean is the previous value)
% model.p_min = 0.4;                      % Minimum beat period
model.a_trans_vr = 0.05^2;              % Beat amplitude transition density (normal) variance (mean is the previous value)
model.w_trans_vr = 0.01*eye(model.dw);  % Waveform transition density (normal) covariance matrix (mean is the previous value)
model.b_trans_mn = 0.2;                 % Breath delay density (exponential) mean

% Clutter
% model.pc = 1E-3;                           % Clutter probability
model.clut_trans = [0.999 0.01; 0.001 0.99]; % Clutter indicator transition matrix
model.y_clut_vr = 50^2;                    % Clutter observation variance

% Observation model
model.y_obs_vr = 0.5^2;                 % Observation variance
