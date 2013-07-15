function [ beat ] = beat_init( model, pre_time, pre_param, ante_param, time, param )
%BEAT_INIT Initialise a heartbeat structure

if isempty(pre_time)
    pre_time = -inf;
end
if isempty(pre_param)
    pre_param = zeros(model.dp,1);
end
if isempty(ante_param)
    ante_param = zeros(model.dp,0);
end
if isempty(time)
    time = zeros(1,0);
end
if isempty(param)
    param = zeros(model.dp,0);
end

beat = struct(...
    'pre_time', pre_time, ...
    'pre_param', pre_param, ...
    'ante_param', ante_param, ...
    'time', time, ...
    'param', param);

end

