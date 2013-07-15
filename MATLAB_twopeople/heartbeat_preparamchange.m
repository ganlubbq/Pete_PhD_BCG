function [ prob ] = heartbeat_preparamchange( model, beat, end_time )
%HEARTBEAT_PREPARAMPROB Calculate the probability associated with a change
%in the beat pre_param value.

% Dig out parameters
% ante_param = beat.ante_param;
pre_param = beat.pre_param;
pre_time = beat.pre_time;

if isempty(beat.time)
    
    % No following beats
    beat_period_offset = end_time - pre_time - pre_param(1);
    prob = log(1 - gamcdf( beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
    
else
    
    % A beat does follow
    [~, prob] = heartbeat_paramtrans(model, pre_param, beat.param(:,1));
    beat_period_offset = beat.time(1) - pre_time - pre_param(1);
    prob = prob + log(gampdf(beat_period_offset, model.tau_trans_shape, model.tau_trans_scale ));
    
end

end

