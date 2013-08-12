function [ ps ] = heartbeat_inference( display, algo, model, time, observ )
%HEARTBEAT_INFERENCE Separate BCG heartbeat signals. Joint estimation of
%beat timings and waveforms, using a variable rate particle filter on a
%rolling window.

Nf = algo.Nf;

% Initialise particle smoother structure
ps = ps_init(model, Nf);

% Initialise a blank particle filter structure
last_pf = pf_init(model, Nf);

fprintf(1, 'Initialisation Loop: ');
% Particle loop
for ii = 1:Nf
    
    fprintf(1, '%u ', ii);
    
    % Sample initial values for the parameters and changepoints
    start_time = 0; end_time = algo.Lstart/model.fs;
    last_pf(ii).beat = heartbeat_initialise( algo, model, time(1:algo.Lstart), observ(:,1:algo.Lstart), start_time, end_time );
    
    % Store
    ps(ii).beat = last_pf(ii).beat;
    
%     % Probabilities
%     [H, Y] = heartbeat_obsmat(algo, model, time(1:algo.Lstart), observ(:,1:algo.Lstart), ps(ii).beat);
%     lhood = loggausspdf(Y, H*model.w_prior_mn, H*model.w_prior_vr*H'+model.y_obs_vr*eye(length(Y)));
    
    % Update weight
    last_pf(ii).weight = 0;%lhood;%
    
end
fprintf(1, ' COMPLETE.\n');

% Particle filter loop
kk = 0;
mm = 0;
while 1
    
    % Time, step size and window length
    mm = mm + 1;
    kk = kk + algo.S;
    if kk >= model.K
        break;
    end
    
    L = min(algo.L, model.K-kk);
    S = min(algo.S, model.K-kk);
    start_time = time(kk+1);
    end_time = time(kk+L);
    if mm == 1
        last_L = algo.Lstart;
    else
        last_L = algo.L;
    end
    
%     if mm == 3
%         Nf = Nf/10;
%     end
        
    if display.text
        % Text output
        fprintf(1, 'Particle filter iteration %u:\n', mm);
        fprintf(1, '     - Window from %f to %f (indexes %u to %u of %u)\n',  start_time, end_time, kk+1, kk+L, model.K);
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
                
        % Find parameters of Gaussian posterior waveform distribution
        [wf_mn, wf_vr] = heartbeat_separation( display, algo, model, time(1:kk), observ(:,1:kk), ps(ii).beat );
        
        % Old probabilities
        if mm > 1
            [H, Y] = heartbeat_obsmat(algo, model, time(kk+1:kk+last_L-S), observ(:,kk+1:kk+last_L-S), ps(ii).beat);
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
            
            % MH step to refresh pre_param
            beat_period_offset = start_time - pre_time - pre_param;
            survive_prob = log(1 - gamcdf( beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
            ppsl_pre_param = heartbeat_paramtrans(model, ante_param, []);
            beat_period_offset = start_time - pre_time - ppsl_pre_param;
            ppsl_survive_prob = log(1 - gamcdf( beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
            if log(rand) < ppsl_survive_prob-survive_prob
                pre_param = ppsl_pre_param;
                
                % Update smoothed particle too
                idx = find( ps(ii).beat(pp).time==pt.beat(pp).pre_time );
                if ~isempty(idx)
                    ps(ii).beat(pp).param(:,idx) = ppsl_pre_param;
                end
            end
            
            [beat(pp), init_trans_prob(pp)] = heartbeat_beattrans(model, pre_time, pre_param, ante_param, start_time, end_time, []);
            
        end
        
        if 1%model.np == 1
            % Optimise and re-propose
            [beat, ppsl_prob] = heartbeat_optimalgaussianproposal(algo, model, beat, time(kk+1:kk+L), observ(:,kk+1:kk+L), wf_mn, wf_vr);
            for pp = 1:model.np
                [~, trans_prob(pp)] = heartbeat_beattrans(model, beat(pp).pre_time, beat(pp).pre_param, beat(pp).ante_param, start_time, end_time, beat(pp));
            end
        else
            trans_prob = 0;
            ppsl_prob = 0;
        end

%         % Composite proposal update
%         if model.np == 1
%             if ~isempty(beat.time)&&(mm>1)
%                 [beat, ppsl_prob_inc] = heartbeat_compositeproposal(algo, model, time(kk+1:kk+L), observ(:,kk+1:kk+L), beat, wf_mn, wf_vr);
%                 [~, trans_prob(pp)] = heartbeat_beattrans(model, beat(pp).pre_time, beat(pp).pre_param, beat(pp).ante_param, start_time, end_time, beat);
%             end
%         else
%             trans_prob = init_trans_prob;
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
        pt.weight = pt.weight - selected_weights(ii) ...
                    + new_lhood - old_lhood ...
                    + sum(trans_prob) - sum(init_trans_prob) ...
                    - ppsl_prob;
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






function [ beat, ppsl_prob ] = heartbeat_optimalgaussianproposal(algo, model, beat, time, observ, wf_mn, wf_vr)

% Create a beat structure with just the first beat of each person
original_beat = beat;
for pp = 1:model.np
    beat(pp).time(2:end) = [];
    beat(pp).param(2:end) = [];
end

if (model.np==1) && ~isempty(beat.time)
    p_idx = 1;
    b_idx = 1;
elseif(model.np==1) && isempty(beat.time)
    p_idx = [];
    b_idx = [];
elseif isempty(beat(1).time) && ~isempty(beat(2).time)
    p_idx = 2;
    b_idx = 1;
elseif ~isempty(beat(1).time) && isempty(beat(2).time)
    p_idx = 1;
    b_idx = 1;
else
    p_idx = [1 2];
    b_idx = [1 1];
end

% [tau_opt, hess] = optimise_beat_times(algo, model, beat, p_idx, b_idx, time, observ, wf_mn, wf_vr);
[tau_opt] = optimise_beat_times(algo, model, beat, p_idx, b_idx, time, observ, wf_mn, wf_vr);

if ~isempty(p_idx)
%     if isposdef(-hess)
%         ppsl_vr = -inv(hess);
%     else
        ppsl_vr = 0.001^2 * eye(length(p_idx));         % Should be using the Hessian, but a fixed value seems to do fine
%     end
    tau = mvnrnd(tau_opt', ppsl_vr)';
    ppsl_prob = loggausspdf(tau, tau_opt, ppsl_vr);
else
    tau = [];
    ppsl_prob = 0;
end

for pp = 1:length(tau)
    beat(p_idx(pp)).time(b_idx(pp)) = tau(pp);
end

% Put the remaining beats back in, optimise marginally and proprose a time
for pp = 1:model.np
    beat(pp).time = [beat(pp).time original_beat(pp).time(2:end)];
    beat(pp).param = [beat(pp).param original_beat(pp).param(2:end)];
end

end

function [tau, hess] = optimise_beat_times(algo, model, beat, p_idx, b_idx, time, observ, wf_mn, wf_vr)

% p_idx is the index(es) of the person for whom a heartbeat needs optimising
% b_idx is the index(es) of the heartbeat
% They should have the same length

if isempty(p_idx)
    % No search required
    
    tau = [];
    hess = [];

elseif length(p_idx) == 1
    % Line search
    
    tau0 = beat(p_idx).time(b_idx);
    
    % Minimise negative log-likelihood
    h_of = @(tau) log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau);
    LB = tau0-0.15;
    UB = tau0+0.15;
%     tau = fminbnd(h_of, LB, UB);
    options = optimset('GradObj','on', 'Display','notify-detailed');%, 'TolX',0.001);
    tau = fmincon(h_of, tau0, [], [], [], [], LB, UB, [], options);
    
    % Approximate Hessian?
    if nargout > 1
        dt = 0.0001;
        [~, grad0] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau);
        [~, grad1] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau+dt);
        hess = (grad1-grad0)/dt;
    else
        hess = [];
    end
    
else
    % 2D search

    for bb = 1:length(p_idx)
        tau0(bb) = beat(p_idx(bb)).time(b_idx(bb));
    end
    tau0 = tau0';
    
    % Minimise negative log-likelihood
    h_of = @(tau) log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau);
    options = optimset('GradObj','on', 'Display','notify-detailed', 'TolX',0.001);
    LB = tau0-0.15;
    UB = tau0+0.15;
    tau = fmincon(h_of, tau0, [], [], [], [], LB, UB, [], options);
    
    % Hessian?
    hess = [];
    
end

end



function [func, grad, hess] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau)

for bb = 1:length(p_idx)
    beat(p_idx(bb)).time(b_idx(bb)) = tau(bb);
end

[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);
R = model.y_obs_vr*eye(length(Y));
S = R + H*wf_vr*H';

func = -loggausspdf(Y, H*wf_mn, S);

if nargout > 1
    grad = zeros(length(p_idx),1);
    for bb = 1:length(p_idx)
        dH = heartbeat_obsmatderiv(algo, model, time, beat, p_idx(bb), b_idx(bb));
        invSigma = inv(wf_vr)+H'*(R\H);
        mu = invSigma\( H'*(R\Y) + wf_vr\wf_mn );
        grad(bb) = -( mu'*dH'*(R\(Y-H*mu)) - trace( H'*(R\dH)/invSigma ) );
    end
end
if nargout > 2
    hess = 0;
end

end

% function [beat, ppsl_prob] = heartbeat_compositeproposal(algo, model, time, observ, beat, wf_mn, wf_vr)
% 
% %%% THIS WILL ONLY WORK IF THERE'S ONLY 1 BEAT IN THE WINDOW %%%
% 
% dl_start = 1E-5;
% dl_min = 1E-8;
% dl_max = 0.5;
% err_thresh = 0.01;
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
%     if ll_count > 50
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
% %     grad_prior = (model.tau_trans_shape-1)/(tau0-tau_shift) - 1/model.tau_trans_scale;
% %     hess_prior = -(model.tau_trans_shape-1)/((tau0-tau_shift)^2);
% %     P = -inv(hess_prior);
% %     m = tau0 + P*grad_prior;
%     P = model.tau_trans_shape*model.tau_trans_scale^2;
%     m = model.tau_trans_shape*model.tau_trans_scale+tau_shift;
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
% %     grad_prior = (model.tau_trans_shape-1)/(tau-tau_shift) - 1/model.tau_trans_scale;
% %     hess_prior = -(model.tau_trans_shape-1)/((tau-tau_shift)^2);
% %     P_new = -inv(hess_prior);
% %     m_new = tau + P*grad_prior;
%     P_new = model.tau_trans_shape*model.tau_trans_scale^2;
%     m_new = model.tau_trans_shape*model.tau_trans_scale;
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

