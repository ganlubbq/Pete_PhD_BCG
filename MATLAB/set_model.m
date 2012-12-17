% Set model parameters

% Basics
model.K = 2000;                                  % Number of observations
model.fs = 30;                                  % Sampling frequency of observations (after load_and_calibrate, which downsamples)
model.dp = 2;                                   % Number of changepoint parameter dimensions (beat period)
model.dw = 40;                                  % Number of samples in the beat waveform

% Priors
% load('template_beat.mat');
% template = [template zeros(1,model.dw-length(template))];
% model.w_prior_mn = template';
% model.w_prior_vr = 0.5^2*eye(model.dw);
model.w_prior_mn = zeros(model.dw,1); model.w_prior_mn(10) = 1; model.w_prior_mn(15) = -1;
model.w_prior_vr = 0.5^2*eye(model.dw); %model.w_prior_vr(10,10) = 0.01;
model.p_prior_shape = 50;
model.p_prior_scale = 0.02;

% Transition models
model.tau_trans_shape = 3;                      % Changepoint time transition density (shifted inverse-gamma) shape paramter
model.tau_trans_scale = 0.4;                    % Changepoint time transition density (shifted inverse-gamma) scale paramter
model.p_trans_scale = 1E-4;                     % Beat period transition density (gamma) scale (shape is the previous value/scale)
model.w_trans_vr = 0.05^2*eye(model.dw);        % Waveform transition density (normal) covariance matrix (mean is the previous value)

% Clutter
model.clut_trans = [0.001 0.9; 0.999 0.1];      % Clutter indicator transition matrix
model.y_clut_vr = 2^2;                         % Clutter observation variance

% Disturbance
% model.dstb_trans = [0 0; 1 1];                  % Disturbance indicator transition matrix
model.dstb_trans = [0.01 0.95; 0.99 0.05];          % Disturbance indicator transition matrix
% model.dstb_prior = 0.2;                         % Prior probability of a disturbance beat
model.y_dstb_vr = 1^2;                         % Disturbance observation variance

% Observation model
model.y_obs_vr = 0.2^2;                         % Observation variance
