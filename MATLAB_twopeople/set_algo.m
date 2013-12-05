% Set algorithm options

algo.Nf = 200;                              % Number of particle filter particles
algo.resam_type = 2;                        % 1 = multinomial. 2 = systematic.
algo.L = 50;                                % Window length
algo.S = 25;                                % Step size

algo.drop_fact = 5;                % Reduce number of particles by this factor after first 2 frames