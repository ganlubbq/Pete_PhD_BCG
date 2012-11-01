function [ interp_vector ] = heartbeat_interpolation( algo, model, t, cp_time )
%HEARTBEAT_INTERPOLATION Construct a vector of interpolation factors which
%allow a heartbeat waveform to be interpolated to any time.

n = 1:model.dw;
interp_vector = sinc((t - cp_time - n/model.fs)*model.fs);

end

