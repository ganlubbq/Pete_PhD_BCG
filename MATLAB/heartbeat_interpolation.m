function [ interp_vector ] = heartbeat_interpolation( algo, model, t, cp_time )
%HEARTBEAT_INTERPOLATION Construct a vector of interpolation factors which
%allow a heartbeat waveform to be interpolated to any time.

t = t(:);
n = 0:model.dw-1;
% interp_vector = sinc((t - cp_time - n/model.fs)*model.fs);
interp_vector = sinc( (bsxfun(@minus, t, n/model.fs)-cp_time) *model.fs);

end

