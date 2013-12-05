function [ prob ] = heartbeat_survivorfunc( model, cp_time, cp_param, time )
%HEARTBEAT_SURVIVORFUNC Summary of this function goes here
%   Detailed explanation goes here

prob = log(1-gamcdf(time-(cp_time+cp_param), model.tau_trans_shape, model.tau_trans_scale ));

end

