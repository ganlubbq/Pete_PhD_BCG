function [ deriv_vector ] = heartbeat_interpolationderiv( algo, model, t, cp_time )
%HEARTBEAT_INTERPOLATIONDERIV Construct a matrix of the derivatives of 
%the interpolation factors which allow a heartbeat waveform to be
%interpolated to any time.

if numel(cp_time)==1
    cp_time = cp_time*ones(size(t));
end

dt = t - cp_time;

dt = dt(:);
n = 1:model.dw;
grid = dt(:,ones(1,model.dw))*model.fs - n(ones(length(t),1),:);

deriv_vector = -model.fs*sincd(grid);
deriv_vector(isinf(grid))=0;

% window = ones(size(dt));
% zto_idx = (t-cp_time)<(1/model.fs);
% window(zto_idx) = 0.5*( 1-cos(pi*model.fs*dt(zto_idx)) );
% window((t-cp_time)<0)=0;
% 
% deriv_vector = window(:,ones(1,model.dw)).*deriv_vector;
% 
% other_term = zeros(size(grid));
% other_term(zto_idx,:) = sinc(grid(zto_idx,:));
% 
% window = ones(size(dt));
% window(zto_idx) = 0.5*pi*model.fs*sin(pi*model.fs*dt(zto_idx));
% window((t-cp_time)<0)=0;
% 
% deriv_vector = deriv_vector - window(:,ones(1,model.dw)).*other_term;

end

