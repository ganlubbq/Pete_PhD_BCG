function [ cp_time, cp_param, prob ] = heartbeat_cptransition( model, last_cp_time, last_cp_param, known_time, current_time, cp_time, cp_param )
%HEARTBEAT_CPTRANSITION Sample and/or calculate the probability for the 
% changepoint transition density.

% prob is a log-probability.
% We know that no changepoint occurs between last_cp_time and known_time.

% Sample state if not provided
if (nargin<7)||isempty(cp_time)||isempty(cp_param)
    
    last_p = last_cp_param(1);
    last_d = last_cp_param(2);
    
    % Rejection sample new changepoint time
    lower_lim = invgamcdf(known_time-last_cp_time-last_p, model.tau_trans_shape, model.tau_trans_scale);
    u = unifrnd(lower_lim, 1);
    cp_time = last_cp_time + last_p + invgaminv(u, model.tau_trans_shape, model.tau_trans_scale);
    
    if cp_time < current_time
        
        % Sample a new parameter for the new changepoint
        cp_param = zeros(model.dp,1);
        cp_param(1) = gamrnd(last_cp_param(1)/model.p_trans_scale, model.p_trans_scale);
        cp_param(2) = rand < model.dstb_trans(1, last_d+1);
        
    else
        
        % No changepoint occured between known_time and current_time
        cp_time = zeros(1,0);
        cp_param = zeros(model.dp,0);
        
    end
end

% Calculate probability if required
if nargout>2
    if isempty(cp_time)
        prob = log(1-invgamcdf(known_time-last_cp_time-last_p, model.tau_trans_shape, model.tau_trans_scale));
    else
        prob = log(invgampdf(cp_time-last_cp_time-last_p, model.tau_trans_shape, model.tau_trans_scale)) ...
              -log(1-invgamcdf(known_time-last_cp_time-last_p, model.tau_trans_shape, model.tau_trans_scale)) ...
              +loggampdf(cp_param(1), last_p/model.p_trans_scale, model.p_trans_scale);
        if cp_param(2)==1
            prob = prob + log(model.dstb_trans(1, last_d+1));
        else
            prob = prob + log(model.dstb_trans(2, last_d+1));
        end
    end
else
    prob = [];
end


end

