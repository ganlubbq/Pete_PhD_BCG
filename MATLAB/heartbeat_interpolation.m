function [ interp_vector ] = heartbeat_interpolation( algo, model, t, cp_time )
%HEARTBEAT_INTERPOLATION Construct a vector of interpolation factors which
%allow a heartbeat waveform to be interpolated to any time.

t = t(:);
n = 0:model.dw-1;
% interp_vector = sinc((t - cp_time - n/model.fs)*model.fs);
% grid = bsxfun(@minus, t*model.fs, n) - cp_time*model.fs;
grid = (t*model.fs*ones(1,model.dw)-ones(length(t),1)*n) - cp_time*model.fs;
interp_vector = sinc(grid);
% neg = abs(grid)>5;
% interp_vector = zeros(size(grid));
% interp_vector(~neg) = sinc(grid(~neg));

end

