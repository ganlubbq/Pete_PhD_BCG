function [ beat, prob ] = heartbeat_beattrans(model, pre_time, pre_param, ante_param, start_time, end_time, beat)
%HEARTBEAT_BEATTRANS Sample a new heartbeat sequence and calculate the
%probability.

if isempty(beat)
    
    % Sample heartbeats through the window
    latest_time = pre_time; latest_param = pre_param;
    period = 0;
    time = []; param = zeros(model.dp,0);
    
    % Loop until we're through the window
    while true
        
        % Rejection sample a beat time
        lower_lim = gamcdf(start_time-latest_time-latest_param, model.tau_trans_shape, model.tau_trans_scale);
        u = unifrnd(lower_lim, 1);
        period = latest_param + gaminv(u, model.tau_trans_shape, model.tau_trans_scale);
        
        latest_time = latest_time + period;
        latest_param = heartbeat_paramtrans(model, latest_param, []);
        
        % Store it
        if latest_time < end_time
            time = [time, latest_time];
            param = [param, latest_param];

%             %%% FUDGE TO STOP MORE THAN 1 BEAT PER WINDOW!!! %%%
%             break;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            break;
        end
        
    end
    
    % Make a beat
    beat = beat_init(model, pre_time, pre_param, ante_param, time, param);
    
end

if nargout > 1
    
    % Calculate probability
    prob = 0;
    latest_time = pre_time; latest_param = pre_param;
    
    for ll = 1:length(beat.time)
        period = beat.time(ll) - latest_time;
        param = beat.param(:,ll);
        [~, period_prob] = heartbeat_periodtrans(model, latest_param, period);
        [~, param_prob] = heartbeat_paramtrans(model, latest_param, param);
        prob = prob + period_prob + param_prob;
        
        latest_time = beat.time(ll);
        latest_param = param;
    end
    
    % End effects
    prob = prob + log(1 - gamcdf( end_time-latest_time, model.tau_trans_shape, model.tau_trans_scale ));
    prob = prob - log(1 - gamcdf( start_time-pre_time, model.tau_trans_shape, model.tau_trans_scale ));
    
end

end

