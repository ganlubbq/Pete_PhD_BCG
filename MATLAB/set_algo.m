% Set algorithm options

algo.Nf = 200;                              % Number of particle filter particles
algo.resam_type = 2;                        % 1 = multinomial. 2 = systematic.
algo.L = 50;                                % Window length
algo.S = 20;                                % Step size

algo.cp_ppsl_lag = algo.L-1;
algo.kd_vr = 0.5/model.fs;

algo.min_noclut = 0;