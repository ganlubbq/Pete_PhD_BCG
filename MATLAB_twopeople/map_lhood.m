%% Two heartbeats

test_beat = beat;
tau1_offset = beat(1).pre_time + beat(1).pre_param;
tau2_offset = beat(2).pre_time + beat(2).pre_param;

tau_rng = 0:0.02:2;
lhood_map = zeros(length(tau_rng));

for xx = 1:length(tau_rng)
    tau1 = tau1_offset + tau_rng(xx);
    for yy = 1:length(tau_rng);
        tau2 = tau2_offset + tau_rng(yy);
        
        test_beat(1).time = tau1;
        test_beat(2).time = tau2;
        
        [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), test_beat);
        lhood_map(xx,yy) = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        
    end
end

figure, surf(tau1_offset+tau_rng, tau2_offset+tau_rng, lhood_map); shading interp;

%% One heartbeat
test_beat = beat;
% test_beat.pre_time = -Inf;
tau1_offset = 0;beat.pre_time + beat(1).pre_param;

tau_rng = -1.5:0.01:3;
lhood_map = zeros(1, length(tau_rng));

for xx = 1:length(tau_rng)
    tau1 = tau1_offset + tau_rng(xx);
    
    test_beat.time = tau1;
    
    [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), test_beat);
    lhood_map(xx) = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
    prior_map(xx) = log(gampdf(tau1-tau1_offset, model.tau_trans_shape, model.tau_trans_scale));
    
    R = model.y_obs_vr*eye(length(Y)) + H*wf_vr*H';
    Sigma = inv( inv(wf_vr)+H'*(R\H) );
    mu = Sigma*( H'*(R\Y) + wf_vr\wf_mn );
%     lhood_deriv_map(xx) = -( mu'*dH'*(R\(Y-H*mu)) - trace( H'*(R\dH)*Sigma ) );
    
end

figure, hold on, plot(tau_rng, lhood_map); plot(tau_rng, prior_map, 'r');
% figure, plot(tau_rng, lhood_deriv_map);

%% Variance

vr_map = zeros(1,length(tau_rng));
for xx = 1:length(tau_rng)
    tau1 = tau1_offset + tau_rng(xx);
    
    test_beat.time = tau1;
    
    [H, Y] = heartbeat_obsmat(algo, model, time(1), observ(:,1), test_beat);
    vr_map(xx) = H(1,:)*wf_vr*H(1,:)'+model.y_obs_vr;
    
end

figure, plot(tau_rng, sqrt(vr_map))

