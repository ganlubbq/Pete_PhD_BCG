function [ cp_time, cp_param, rb_mn, rb_vr ] = heartbeat_rm( algo, model, kk, cp_time, cp_param, cp_rb_mn, cp_rb_vr, last_cp_time, last_cp_param, current_time, time, observ, rb_mn, rb_vr, cp_idx_start, clut_hist )
%HEARTBEAT_RM Do some resample move on the latest changepoint

% Propose a new value for the last changepoint time, restricting it to lie
% between the same pair of observations

if cp_time > 0
    
    clut_indic = zeros(size(time));
    clut_indic(cp_idx_start:kk) = flipud( clut_hist(1:kk+1-cp_idx_start) );
    
    a = cp_param(2);
    
    % Make a vector of observations since the last changepoint and interpolate
    % the expected signal at these points
    to_use = (clut_indic==0)&(time<=current_time)&(time>cp_time);
    t_vec =   time( to_use );
    y_vec = observ( to_use )';
    
    % Propose a change to the time
    t_min = max(time(time<cp_time));
    t_max = min(time(time>cp_time));
    ppsl_cp_time = unifrnd(t_min, t_max);
    
    % Calculate likelihoods
    old_H = heartbeat_interpolation(algo, model, t_vec, cp_time);
    new_H = heartbeat_interpolation(algo, model, t_vec, ppsl_cp_time);
    new_lhood = loggausspdf(y_vec, a*new_H*cp_rb_mn, a^2*new_H*cp_rb_vr*new_H'+model.y_obs_vr*eye(length(t_vec)));
    old_lhood = loggausspdf(y_vec, a*old_H*cp_rb_mn, a^2*old_H*cp_rb_vr*old_H'+model.y_obs_vr*eye(length(t_vec)));
    
    % Calculate priors
    last_p = last_cp_param(1);
    last_b = last_cp_param(3);
    new_trans_prob = loggampdf(ppsl_cp_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale);
    old_trans_prob = loggampdf(cp_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale);
    
    % Posterior
    new_post_prob = new_lhood + new_trans_prob;
    old_post_prob = old_lhood + old_trans_prob;
    
    % Acceptance?
    if log(rand) < new_post_prob - old_post_prob        % Uniform proposal cancels out
        cp_time = ppsl_cp_time;
    end
    
end



% Propose a new value for the beat period from the transition density
% (Gibbs, so no acceptance probability needed.
if any(isnan(last_cp_param))
    ppsl_p = gamrnd(model.p_prior_shape, model.p_prior_scale);
else
    ppsl_p = gamrnd(last_cp_param(1)/model.p_trans_scale, model.p_trans_scale);
end

b = cp_param(3);
% b = 0;

new_surv_prob = log(1-gamcdf(current_time-cp_time, (ppsl_p+b)/model.tau_trans_scale, model.tau_trans_scale));
old_surv_prob = log(1-gamcdf(current_time-cp_time, (cp_param(1)+b)/model.tau_trans_scale, model.tau_trans_scale));

if log(rand) < new_surv_prob-old_surv_prob
    cp_param(1) = ppsl_p;
end

% Propose a new value for the beat period adjustment from the transition density
% (Gibbs, so no acceptance probability needed.
if any(isnan(last_cp_param))
    ppsl_b = raylrnd(model.b_prior_mn);
else
    ppsl_b = invgamrnd(model.b_trans_shape, model.b_trans_scale);
end

p = cp_param(1);

new_surv_prob = log(1-gamcdf(current_time-cp_time, (p+ppsl_b)/model.tau_trans_scale, model.tau_trans_scale));
old_surv_prob = log(1-gamcdf(current_time-cp_time, (p+cp_param(3))/model.tau_trans_scale, model.tau_trans_scale));

if log(rand) < new_surv_prob-old_surv_prob
    cp_param(3) = ppsl_b;
end

% % Propose a new value for the beat amplitude using MH
% 
% if ~isnan(last_cp_time)
%     
%     % Make a vector of observations since the last changepoint and interpolate
%     % the expected signal at these points
%     t_vec = time( (time<current_time)&(time>cp_time) );
%     y_vec = observ( (time<current_time)&(time>cp_time) )';
%     H = heartbeat_interpolation(algo, model, t_vec, cp_time);
%     
%     % Proposal
%     [ppsl_mn, ppsl_vr] = kf_update(last_cp_param(2), model.a_trans_vr, y_vec, H*cp_rb_mn, model.y_obs_vr*eye(length(t_vec)));
%     ppsl_a = mvnrnd(ppsl_mn, ppsl_vr);
%     frwd_ppsl_prob = loggausspdf(ppsl_a, ppsl_mn, ppsl_vr);
%     bkwd_ppsl_prob = loggausspdf(cp_param(2), ppsl_mn, ppsl_vr);
%     
%     % Posterior
%     frwd_post_prob = loggausspdf(y_vec, ppsl_a*H*cp_rb_mn, model.y_obs_vr*eye(length(t_vec))+ppsl_a^2*H*cp_rb_vr*H') ...
%         + loggausspdf(ppsl_a, last_cp_param(2), model.a_trans_vr);
%     bkwd_post_prob = loggausspdf(y_vec, cp_param(2)*H*cp_rb_mn, model.y_obs_vr*eye(length(t_vec))+cp_param(2)^2*H*cp_rb_vr*H') ...
%         + loggausspdf(cp_param(2), last_cp_param(2), model.a_trans_vr);
%     
%     % Acceptance
%     ap = (frwd_post_prob-frwd_ppsl_prob)-(bkwd_post_prob-bkwd_ppsl_prob);
%     if log(rand) < ap
%         cp_param(2) = ppsl_a;
%         
%         % Update latest rb estimate - THIS OUGHT TO TAKE ACCOUNT OF
%         CLUTTER
%         [rb_mn, rb_vr] = kf_update(cp_rb_mn, cp_rb_vr, y_vec, ppsl_a*H, model.y_obs_vr*eye(length(t_vec)));
%         
%     end
%     
% end

% % First sample a template
% w = mvnrnd(cp_rb_mn', cp_rb_vr)';
%
% % Make a vector of observations since the last changepoint and interpolate
% % the expected signal at these points
% t_vec = time( (time<current_time)&(time>cp_time) );
% y_vec = observ( (time<current_time)&(time>cp_time) )';
% H = zeros(length(t_vec), model.dw);
% for ii = 1:length(t_vec)
%     H(ii,:) = heartbeat_interpolation(algo, model, t_vec(ii), cp_time);
% end
%
% % Use this to do a KF update
% [a_mn, a_vr] = kf_update(last_cp_param(2), model.a_trans_vr, y_vec, H*w, model.y_obs_vr*eye(length(t_vec)));
%
% % Sample this distribution
% cp_param(2) = mvnrnd(a_mn, a_vr);

end

