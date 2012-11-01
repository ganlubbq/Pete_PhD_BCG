function [ cp_time, cp_param, prob ] = heartbeat_cptransition( model, last_cp_time, last_cp_param, known_time, current_time, cp_time, cp_param )
%HEARTBEAT_CPTRANSITION Sample and/or calculate the probability for the 
% changepoint transition density.

% prob is a log-probability.
% We know that no changepoint occurs between last_cp_time and known_time.

% Sample state if not provided
if (nargin<7)||isempty(cp_time)||isempty(cp_param)
    
    last_p = last_cp_param(1);
    
    % Rejection sample new changepoint time
    lower_lim = gamcdf(known_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale);
    if lower_lim < 1
        u = unifrnd(lower_lim, 1);
        cp_time = last_cp_time + gaminv(u, last_p/model.tau_trans_scale, model.tau_trans_scale);
    else
        % Distribution is pathologically peaky
        cp_time = last_cp_time + last_p;
    end
    
    if cp_time < current_time
        
        % Sample a new parameter for the new changepoint
        cp_param = mvnrnd(last_cp_param, diag([model.p_trans_vr,  model.a_trans_vr]))';
        
    else
        
        % No changepoint occured between known_time and current_time
        cp_time = zeros(1,0);
        cp_param = zeros(model.dp,0);
        
    end
end

% Calculate probability if required
if nargout>2
    if isempty(cp_time)
        prob = log(1-gamcdf(known_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale));
    else
        prob = log(gampdf(cp_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale)) ...
              -log(1-gamcdf(known_time-last_cp_time, last_p/model.tau_trans_scale, model.tau_trans_scale)) ...
              +loggausspdf(cp_param, last_cp_param, diag([model.p_trans_vr,  model.a_trans_vr]));
    end
else
    prob = [];
end


end

