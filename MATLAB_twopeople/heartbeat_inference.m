function [ ps ] = heartbeat_inference( display, algo, model, time, observ )
%HEARTBEAT_INFERENCE Separate BCG heartbeat signals. Joint estimation of
%beat timings and waveforms, using a variable rate particle filter on a
%rolling window.

% Make local copies of useful numbers
K = model.K;        % Number of time steps
Nf = algo.Nf;       % Number of filter particles
S = algo.S;         % Time increment (number of time steps to move forward each processing frame)
L = algo.L;         % Window length

% Number of particle filter steps we need to run
M = ceil(K/S);

% Initialise particle smoother structure
ps = ps_init(model, Nf);

% Initialise a blank particle filter structure
last_pf = pf_init(model, Nf);

% fprintf(1, 'Initialisation Loop: ');
% Particle loop
for ii = 1:Nf
    
%     fprintf(1, '%u ', ii);
%     
%     % Sample initial values for the parameters and changepoints
%     start_time = 0; end_time = algo.Lstart/model.fs;
%     last_pf(ii).beat = heartbeat_initialise( algo, model, time(1:algo.Lstart), observ(:,1:algo.Lstart), start_time, end_time );
    
    % People loop
    for pp = 1:model.np
        
        % Sample initial values for the parameters and changepoints from the prior
        start_time = 0; end_time = algo.L/model.fs;
        beat(pp) = heartbeat_beatprior(model);
        last_pf(ii).beat(pp) = heartbeat_beattrans(model, beat(pp).pre_time, beat(pp).pre_param, [], start_time, end_time, []);
        ps(ii).beat(pp).pre_time = last_pf(ii).beat(pp).pre_time;
        ps(ii).beat(pp).pre_param = last_pf(ii).beat(pp).pre_param;

    end
    
    % Probabilities
    [H, Y] = heartbeat_obsmat(algo, model, time(1:L), observ(:,1:L), ps(ii).beat);
    lhood = loggausspdf(Y, H*model.w_prior_mn, H*model.w_prior_vr*H'+model.y_obs_vr*eye(length(Y)));
    
    % Update weight
    last_pf(ii).weight = lhood;
    
end
% fprintf(1, ' COMPLETE.\n');

% Particle filter loop
for mm = 1:M-1
    
%     if mm == 3
%         Nf = Nf/10;
%     end
    
    % Time, step size and window length
    if mm == 1, kk = algo.S; else kk = kk + S; end
    L = min(algo.L, K-kk);
    S = min(algo.S, K-kk);
    start_time = time(kk+1);
    end_time = time(kk+L);
    
%     if mm < 4
%         L = 2*L;
%     end
        
    if display.text
        % Text output
        fprintf(1, 'Particle filter iteration %u of %u:\n', mm, M);
        fprintf(1, '     - Window from %f to %f (indexes %u to %u)\n',  start_time, end_time, kk+1, kk+L);
        tic;
    end
    
    % Initialise particle filter structure
    pf = pf_init(model, Nf);
    
    % Particle selection
    if mm > 0
        auxiliary_weights = [last_pf.weight];
    else
        auxiliary_weights = zeros(1,Nf);
    end
    [ancestor, selected_weights] = sample_weights(auxiliary_weights, Nf, algo.resam_type);
    ps = ps(ancestor);
    
    % Particle loop
    for ii = 1:Nf
        
        % Which ancestor?
        anc = ancestor(ii);
        
        % Initialise particle
        pt = pf_forwardparticle(model, anc, last_pf(anc), start_time, end_time);
        
        % Propose new parameters for the preceeding heartbeats from a MH kernel
        for pp = 1:model.np
            
            if pt.beat(pp).pre_time < 0
                
                beat = pt.beat(pp);
                ante_param = beat.ante_param;
                
                % Sample a new value for pre_param
                new_pre_param = heartbeat_paramtrans(model, ante_param, []);
                new_beat = beat;
                new_beat.pre_param = new_pre_param;
                
                % Probabilities
                new_pre_param_prob = heartbeat_preparamchange(model, new_beat, start_time);
                old_pre_param_prob = heartbeat_preparamchange(model, beat,     start_time);
                
                % Accept?
                if log(rand) < new_pre_param_prob-old_pre_param_prob
                    pt.beat(pp) = new_beat;
                    
                    % Update smoothed particle too
                    idx = find( ps(ii).beat(pp).time==pt.beat(pp).pre_time );
                    if ~isempty(idx)
                        ps(ii).beat(pp).param(:,idx) = new_pre_param;
                    end
                end
                
            end
            
        end
        
        % Find parameters of Gaussian posterior waveform distribution
        [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time(1:kk), observ(:,1:kk), ps(ii).beat );
        
        % Old probabilities
        if mm > 0
            [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L-S), observ(:,kk+1:kk+L-S), ps(ii).beat);
            old_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        else
            old_lhood = 0;
        end
        
        % Propose new heartbeats for the window
        beat = pt.beat;
        for pp = 1:model.np
            
            pre_time = pt.beat(pp).pre_time;
            pre_param = pt.beat(pp).pre_param;
            ante_param = pt.beat(pp).ante_param;
            [beat(pp), trans_prob(pp)] = heartbeat_beattrans(model, pre_time, pre_param, ante_param, start_time, end_time, []);
            
        end
        
        if 0;%(mm>2)&&(~isempty(beat.time))%
            % Optimise and re-propose
            [beat, ppsl_prob] = heartbeat_optimalgaussianproposal(algo, model, beat, time(kk+1:kk+L), observ(:,kk+1:kk+L), wf_mn, wf_vr);
            [~, trans_prob] = heartbeat_beattrans(model, beat.pre_time, beat.pre_param, beat.ante_param, start_time, end_time, beat);
        else
            ppsl_prob = trans_prob;
        end
        
%         % Composite proposal update
%         if model.np == 1
%             if ~isempty(beat.time)
%                 [beat, ppsl_prob_inc] = heartbeat_compositeproposal(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), beat, wf_mn, wf_vr);
%             end
%         else
%             ppsl_prob_inc = 0;
%         end
        
        % Store it
        pt.beat = beat;
        for pp = 1:model.np
            % Update smoothed particles too
            idx = find( ps(ii).beat(pp).time>start_time );
            ps(ii).beat(pp).time(idx) = [];
            ps(ii).beat(pp).param(:,idx) = [];
            ps(ii).beat(pp).time = [ps(ii).beat(pp).time, beat(pp).time];
            ps(ii).beat(pp).param = [ps(ii).beat(pp).param, beat(pp).param];
        end

        % Probabilities
        [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), ps(ii).beat);
        new_lhood = loggausspdf(Y, H*wf_mn, H*wf_vr*H'+model.y_obs_vr*eye(length(Y)));
        
        % Update weight
        pt.weight = pt.weight - selected_weights(ii) + new_lhood - old_lhood + sum(trans_prob)-sum(ppsl_prob);
        if isinf(pt.weight)||isnan(pt.weight)
            pt.weight = -inf;
        end
        
        % Store particle
        pf(ii) = pt;
        
        % Plot
        if display.plot_during
            x_vec = H*wf_mn;
            x = reshape(x_vec,4,L);
            figure, plot(time(kk+1:kk+L), x');
        end
        
    end

    if display.text
        % Text output
        fprintf(1, '     - Effective sample size: %f.\n', calc_ESS([pf.weight]));
        fprintf(1, '     - That took %f seconds.\n', toc);
    end
    
    % Keep things that we'll need in the next interation
    last_pf = pf;
    
end

% Final smoother particle selection step

end

% function [beat, ppsl_prob] = heartbeat_compositeproposal(algo, model, time, observ, beat, wf_mn, wf_vr)
% 
% %%% THIS WILL ONLY WORK IF THERE'S ONLY 1 BEAT IN THE WINDOW %%%
% 
% dl_start = 1E-5;
% dl_min = 1E-8;
% dl_max = 0.5;
% err_thresh = 0.1;
% dl_sf = 0.8;
% dl_pow = 0.7;
% 
% tau_shift = beat.pre_time + beat.pre_param;
% 
% % Initialise loop variables
% ppsl_prob = 0;
% dl = dl_start;
% lam = 0;
% 
% lam_evo = 0;
% tau_evo = beat.time;
% 
% % Loop
% ll_count = 0;
% while lam < 1
%     
%     if ll_count > 10
%         ppsl_prob = 1E10;
%         break;
%     end
%     
%     % Pseudo-time and step-size
%     lam0 = lam;
%     lam1 = lam + dl;
%     if lam1 > 1
%         lam1 = 1;
%     end
%     
%     % Starting point
%     tau0 = beat.time;
%     
%     % Interpolation
%     [intvec, Y, d_intvec] = heartbeat_obsmat(algo, model, time, observ, beat);
% 
%     % Linearise observation model
%     R = model.y_obs_vr*eye(length(Y)) + intvec*wf_vr*intvec';
%     Yhat = Y - intvec*wf_mn + d_intvec*wf_mn*tau0;
%     H = d_intvec*wf_mn;
%     
%     % Linearise prior
%     grad_prior = (model.tau_trans_shape-1)/(tau0-tau_shift) - 1/model.tau_trans_scale;
%     hess_prior = -(model.tau_trans_shape-1)/((tau0-tau_shift)^2);
%     P = -inv(hess_prior);
%     m = tau0 + P*grad_prior;
%     
%     % Analytical flow
%     [ tau, prob_ratio, drift] = linear_flow_move( lam1, lam0, tau0, m, P, Yhat, H, R);
%     
%     % Error estimate
%     beat_new = beat; beat_new.time = tau;
%     [interp_new, Y_new, d_interp_new] = heartbeat_obsmat(algo, model, time, observ, beat_new);
%     R_new = model.y_obs_vr*eye(length(Y)) + interp_new*wf_vr*interp_new';
%     Yhat_new = Y_new - interp_new*wf_mn + d_interp_new*wf_mn*tau;
%     H_new = d_interp_new*wf_mn;
%     grad_prior = (model.tau_trans_shape-1)/(tau-tau_shift) - 1/model.tau_trans_scale;
%     hess_prior = -(model.tau_trans_shape-1)/((tau-tau_shift)^2);
%     P_new = -inv(hess_prior);
%     m_new = tau + P*grad_prior;
% 
%     [drift_new] = linear_drift( lam1, tau, m_new, P_new, Yhat_new, H_new, R_new );
%     
%     err_est = 0.5*(lam1-lam0)*(drift_new-drift);
%     err_crit = err_est'*err_est;
%     
%     % Step size adjustment
%     if (err_crit > err_thresh) || (lam1 < 1)
%         dl = min(dl_max, min(10*dl, dl_sf * (err_thresh/err_crit)^dl_pow * dl));
%         if dl < dl_min
%             warning('nlng_smoothupdatebyparticle:ErrorTolerance', 'Minimum step size reached. Local error tolerance exceeded.');
%             ppsl_prob = 1E10;
%             break;
%         end
%     end
%     
%     % Accept/reject step
%     if err_crit < err_thresh
%         
%         ll_count = ll_count + 1;
%         
%         % Update time
%         lam = lam1;
%         
%         % Update state
%         beat.time = tau;
%         
%         lam_evo = [lam_evo lam];
%         tau_evo = [tau_evo tau];
%         
%         % Update probability
%         ppsl_prob = ppsl_prob - log(prob_ratio);
%         
%     else
%         
% %         disp('Error too large. Reducing step size');
%         
%     end
%     
% end
% 
% end

function [ beat, ppsl_prob ] = heartbeat_optimalgaussianproposal(algo, model, beat, time, observ, wf_mn, wf_vr)

%%% JUST FOR ONE HEARTBEAT AT THE MOMENT %%%

tau0 = beat.time;

% Maximise OID
h_of = @(tau) log_oid_with_grad(algo, model, time, observ, beat, wf_mn, wf_vr, tau);
% options = optimset('GradObj','on', 'Display','notify-detailed', 'TolX',0.001);
options = optimset('Display','off', 'TolX',0.01);
warning('off','optim:fminunc:SwitchingMethod');
tau_opt = fminunc(h_of, tau0, options);

% % Which way shall we look?
% [func0, grad0] = log_oid_with_grad(algo, model, time, observ, beat, wf_mn, wf_vr, tau0);
% dir = sign(grad0);
% 
% % A couple of points
% Tpp = 0.25;
% tau1 = tau0+0.25*dir*Tpp;
% tau2 = tau0+0.5*dir*Tpp;
% [func1] = log_oid_with_grad(algo, model, time, observ, beat, wf_mn, wf_vr, tau1);
% [func2] = log_oid_with_grad(algo, model, time, observ, beat, wf_mn, wf_vr, tau2);
% 
% % Fit a parabola
% cub_fit = [tau0^2 tau0^2 tau0 1; tau1^2 tau1^2 tau1 1; tau2^2 tau2^2 tau2 1; 3*tau0^2 2*tau0 1 0]\[func0; func1; func2; grad0];
% tau_opt(1) = (-2*cub_fit(2)+sqrt(4*cub_fit(2)^2-12*cub_fit(1)*cub_fit(3)))/(6*cub_fit(1));
% tau_opt(2) = (-2*cub_fit(2)-sqrt(4*cub_fit(2)^2-12*cub_fit(1)*cub_fit(3)))/(6*cub_fit(1));
% 
% % Find the best peak
% tau_set = [tau0, tau1, tau2];
% func_set = [func0, func1, func2];
% tau_opt( isnan(tau_opt)|~isreal(tau_opt)|(tau_opt<min(tau_set))|(tau_opt>max(tau_set)) ) = [];
% if isempty(tau_opt)
%     tau_opt = tau_set( func_set==max(func_set) );
% end
% if length(tau_opt)>1
%     tau_opt = tau_opt( abs(tau_opt-tau0)==min(abs(tau_opt-tau0)) );
% end

% Propose a new time
ppsl_vr = 0.002^2;
tau = mvnrnd(tau_opt, ppsl_vr);
time_ppsl_prob = loggausspdf(tau, tau_opt, ppsl_vr);

% Propose a new parameter
[param, param_ppsl_prob] = heartbeat_paramtrans(model, beat.pre_param, []);

% Store
beat.time = tau;
beat.param = param;

ppsl_prob = time_ppsl_prob + param_ppsl_prob;

end

function [func, grad, hess] = log_oid_with_grad(algo, model, time, observ, beat, wf_mn, wf_vr, tau)

beat.time = tau;
[H, Y, dH] = heartbeat_obsmat(algo, model, time, observ, beat);
R = model.y_obs_vr*eye(length(Y)) + H*wf_vr*H';

func = -loggausspdf(Y, H*wf_mn, R);
if nargout > 1
    Sigma = inv( inv(wf_vr)+H'*(R\H) );
    mu = Sigma*( H'*(R\Y) + wf_vr\wf_mn );
    grad = -( mu'*dH'*(R\(Y-H*mu)) - trace( H'*(R\dH)*Sigma ) );
end
if nargout > 2
    hess = 0;
end

end



