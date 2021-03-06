function [ interp_vector ] = heartbeat_interpolation( algo, model, t, cp_time, ante_cp_time )
%HEARTBEAT_INTERPOLATION Construct a matrix of interpolation factors which
%allow a heartbeat waveform to be interpolated to any time.

if numel(cp_time)==1
    cp_time = cp_time*ones(size(t));
end

dt = t - cp_time;

dt = dt(:);
n = 1:model.dw;
grid = dt(:,ones(1,model.dw))*model.fs - n(ones(length(t),1),:);

interp_vector = sinc(grid);
interp_vector(isinf(grid))=0;

window = ones(size(dt));
zto_idx = (t-cp_time)<(1/model.fs);
window(zto_idx) = 0.5*( 1-cos(pi*model.fs*(t(zto_idx)-cp_time(zto_idx))) );
window((t-cp_time)<0)=0;

dt2 = t - ante_cp_time;
dt2 = dt2(:);
grid2 = dt2(:,ones(1,model.dw))*model.fs - n(ones(length(t),1),:);
blend_term = zeros(size(interp_vector));
blend_term(zto_idx,:) = sinc(grid2(zto_idx,:));
blend_term(isinf(grid2))=0;

interp_vector = window(:,ones(1,model.dw)).*interp_vector + (1-window(:,ones(1,model.dw))).*blend_term;

% neg = abs(grid)>5;
% interp_vector = zeros(size(grid));
% interp_vector(~neg) = sinc(grid(~neg));
% interp_vector(isinf(grid))=0;

end

