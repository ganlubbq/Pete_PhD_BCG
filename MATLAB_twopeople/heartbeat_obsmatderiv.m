function [ dH ] = heartbeat_obsmatderiv( algo, model, time, beat, p_idx, b_idx )
%HEARTBEAT_OBSMATDERIV

last_beat = -Inf(size(time));
ante_beat = -Inf(size(time));

tau_list = [beat(p_idx).pre_time beat(p_idx).time]';
tdiff = time(ones(length(tau_list),1),:)-tau_list(:,ones(length(time),1));
tdiff(tdiff<0)=0;
if length(time) > 0
    preceding_beat_flag = tdiff & not([zeros(length(tau_list),1) tdiff(:,1:end-1)]);
else
    preceding_beat_flag = false(size(tdiff));
end
[tau_idx, t_idx] = find(preceding_beat_flag);
t_idx = [t_idx; length(time)+1];

tau_idx = tau_idx - 1;
K = find(tau_idx == b_idx);
if ~isempty(K)
    last_beat(t_idx(K):t_idx(K+1)-1) = beat(p_idx).time(tau_idx(K));
    if tau_idx(K) > 1
        ante_beat(t_idx(K):t_idx(K+1)-1) = beat(p_idx).time(tau_idx(K)-1);
    else
        ante_beat(t_idx(K):t_idx(K+1)-1) = beat(p_idx).pre_time;
    end
end

% for nn = 1:length(time)
%     [tauK, K] = max(beat(p_idx).time( beat(p_idx).time<=time(nn) ));
%     if (~isempty(K))&&(K==b_idx)
%         last_beat(nn) = tauK;
%         if K > 1
%             ante_beat(nn) = beat(p_idx).time(K-1);
%         else
%             ante_beat(nn) = beat(p_idx).pre_time;
%         end
%     else
%         last_beat(nn) = -inf;
%         ante_beat(nn) = -inf;
%     end
%     
% end

interp_mat = heartbeat_interpolationderiv(algo, model, time, last_beat, ante_beat);

tmp = repmat({interp_mat},model.num_sens,1);
dH = blkdiag(tmp{:});

if model.np == 2
    if p_idx == 1
        dH = [dH, zeros(size(dH))];
    elseif p_idx == 2
        dH = [zeros(size(dH)), dH];
    end
end

end

