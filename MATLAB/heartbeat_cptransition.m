function [ cp_time, cp_param, prob ] = heartbeat_cptransition( model, last_cp_time, last_cp_param, known_time, current_time, cp_time, cp_param )
%HEARTBEAT_CPTRANSITION Sample and/or calculate the probability for the 
% changepoint transition density.

% prob is a log-probability.
% We know that no changepoint occurs between last_cp_time and known_time.

% Sample state if not provided
if (nargin<7)||isempty(cp_time)||isempty(cp_param)
    
    last_p = last_cp_param(1);
    last_b = last_cp_param(3);
%     last_b = 0;
    
    % Rejection sample new changepoint time
    lower_lim = gamcdf(known_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale);
    if lower_lim < 1
        u = unifrnd(lower_lim, 1);
        cp_time = last_cp_time + gaminv(u, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale);
    else
        % Distribution is pathologically peaky
        cp_time = last_cp_time + (last_p+last_b);
    end
    
    if cp_time < current_time
        
        % Sample a new parameter for the new changepoint
%         cp_param = zeros(model.dp,1);
%         while cp_param <= model.p_min
%             cp_param(1) = mvnrnd(last_cp_param(1), model.p_trans_vr);
%         end
        cp_param(1) = gamrnd(last_cp_param(1)/model.p_trans_scale, model.p_trans_scale);
        cp_param(2) = mvnrnd(last_cp_param(2), model.a_trans_vr);
        cp_param(3) = raylrnd(model.b_trans_mn);
%         cp_param(3) = 0;
        
    else
        
        % No changepoint occured between known_time and current_time
        cp_time = zeros(1,0);
        cp_param = zeros(model.dp,0);
        
    end
end

% Calculate probability if required
if nargout>2
    if isempty(cp_time)
        prob = log(1-gamcdf(known_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale));
    else
        prob = loggampdf(cp_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale) ...
              -log(1-gamcdf(known_time-last_cp_time, (last_p+last_b)/model.tau_trans_scale, model.tau_trans_scale)) ...
              +loggampdf(cp_param(1), last_cp_param(1)/model.p_trans_scale, model.p_trans_scale) ...
              +loggausspdf(cp_param(2), last_cp_param(2), model.a_trans_vr) ...
              +log(exppdf(cp_param(3), model.b_trans_mn));
%         prob = log(gampdf(cp_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale)) ...
%               -log(1-gamcdf(known_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale)) ...
%               +loggausspdf(cp_param, last_cp_param, diag([model.p_trans_vr,  model.a_trans_vr])) ...
%               -log(1-normcdf(model.p_min, last_cp_param(1), model.p_trans_vr));
    end
else
    prob = [];
end


end

