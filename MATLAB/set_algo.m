% Set algorithm options

algo.Nf = 200;           % Number of particle filter particles
algo.resam_type = 1;      % 0 = multinomial. 1 = systematic.
algo.no_resamp_period = 1.5*model.dw;    % Number of frames before resampling starts