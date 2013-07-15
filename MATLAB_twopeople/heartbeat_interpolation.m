function [ interp_vector ] = heartbeat_interpolation( algo, model, t, cp_time )
%HEARTBEAT_INTERPOLATION Construct a matrix of interpolation factors which
%allow a heartbeat waveform to be interpolated to any time.

dt = t - cp_time;

dt = dt(:);
n = 0:model.dw-1;
grid = dt(:,ones(1,model.dw))*model.fs - n(ones(length(t),1),:);

interp_vector = sinc(grid);
interp_vector(isinf(grid))=0;

% neg = abs(grid)>5;
% interp_vector = zeros(size(grid));
% interp_vector(~neg) = sinc(grid(~neg));
% interp_vector(isinf(grid))=0;

end

