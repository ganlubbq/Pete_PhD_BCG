function [ cp_time, cp_param, rb_mn, rb_vr ] = heartbeat_rm( algo, model, cp_time, cp_param, cp_rb_mn, cp_rb_vr, last_cp_time, last_cp_param, current_time, time, observ, rb_mn, rb_vr )
%HEARTBEAT_RM Do some resample move on the latest changepoint

% Propose a new value for the beat period from the transition density
% (Gibbs, so no acceptance probability needed.

if any(isnan(last_cp_param))
    ppsl_p = gamrnd(model.p_prior_shape, model.p_prior_scale);
else
    ppsl_p = gamrnd(last_cp_param(1)/model.p_trans_scale, model.p_trans_scale);
end

new_surv_prob = log(1-gamcdf(current_time-cp_time, ppsl_p/model.tau_trans_scale, model.tau_trans_scale));
old_surv_prob = log(1-gamcdf(current_time-cp_time, cp_param(1)/model.tau_trans_scale, model.tau_trans_scale));

if log(rand) < new_surv_prob-old_surv_prob
    cp_param(1) = ppsl_p;
end

% Propose a new value for the beat amplitude using MH

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

