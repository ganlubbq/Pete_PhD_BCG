function [ dH ] = heartbeat_obsmatderiv( algo, model, time, beat, p_idx, b_idx )
%HEARTBEAT_OBSMATDERIV

last_beat = zeros(size(time));
ante_beat = zeros(size(time));
for nn = 1:length(time)
    [tauK, K] = max(beat(p_idx).time( beat(p_idx).time<=time(nn) ));
    if (~isempty(K))&&(K==b_idx)
        last_beat(nn) = tauK;
        if K > 1
            ante_beat(nn) = beat(p_idx).time(K-1);
        else
            ante_beat(nn) = beat(p_idx).pre_time;
        end
    else
        last_beat(nn) = -inf;
        ante_beat(nn) = -inf;
    end
    
end

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

