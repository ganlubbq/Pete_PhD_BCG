% Set algorithm options

algo.Nf = 200;                              % Number of particle filter particles
algo.resam_type = 2;                        % 1 = multinomial. 2 = systematic.
algo.Nps = 100;                             % Number of pilot samples
algo.L = 50;                                % Pilot sampling horizon
algo.S = 25;                                % Observation batch size