function [ cp_time, cp_param, prob ] = heartbeat_cpprior( model, known_time, current_time, cp_time, cp_param )
%HEARTBEAT_CPPRIOR Sample and/or calculate the probability for the 
% changepoint prior density.

% prob is a log-probability.
% We know that the previous changepoint occured before known_time

% Sample changepoint if not provided
if (nargin<7)||isempty(cp_time)||isempty(cp_param)
    
    % Sample previous parameter from prior
    [last_cp_param, ~] = heartbeat_cpparamprior(model);
    
    % Sample interval from transition density
    [tau01, cp_param, ~]  = heartbeat_cptransition(model, 0, last_cp_param, 0, inf);
    
    % Sample a time for the changepoint before known_time
    last_cp_time = unifrnd(known_time-tau01, known_time);
    
    % Find the changepoint time after time 0
    cp_time = last_cp_time + tau01;
    if cp_time > current_time
        cp_time = [];
        cp_param = [];
    end
    
end

prob = NaN; %%%%% YOU HAVEN'T WRITTEN THIS YET!!! %%%%%

end

