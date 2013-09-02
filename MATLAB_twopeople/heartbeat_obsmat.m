function [ H, Y ] = heartbeat_obsmat( algo, model, time, observ, beat )
%HEARTBEAT_OBSMAT Build the linear interpolation matrix and vectorise the
%corresponding observations

last_beat = cell(model.np,1);
ante_beat = cell(model.np,1);
interp_mat = cell(model.np,1);
vec_interp_mat = cell(model.np,1);
for pp = 1:model.np
    last_beat{pp} = zeros(1,length(time));
    ante_beat{pp} = zeros(1,length(time));
    
    tau_list = [beat(pp).pre_time beat(pp).time]';
    tdiff = time(ones(length(tau_list),1),:)-tau_list(:,ones(length(time),1));
    tdiff(tdiff<0)=0;
    if length(time) > 0
        preceding_beat_flag = tdiff & not([zeros(length(tau_list),1) tdiff(:,1:end-1)]);
    else
        preceding_beat_flag = false(size(tdiff));
    end
    [tau_idx, t_idx] = find(preceding_beat_flag);
    t_idx = [t_idx; length(time)+1];
    
    for kk = 1:length(tau_idx)
        if tau_idx(kk) == 1
            last_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).pre_time;
            ante_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).pre_time;
        elseif tau_idx(kk) == 2
            last_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).time(1);
            ante_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).pre_time;
        else
            last_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).time(tau_idx(kk)-1);
            ante_beat{pp}(t_idx(kk):t_idx(kk+1)-1) = beat(pp).time(tau_idx(kk)-2);
        end
    end
    
%     for nn = 1:length(time)
%         [tauK, K] = max(beat(pp).time( beat(pp).time<=time(nn) ));
%         if ~isempty(K)
%             last_beat{pp}(nn) = tauK;
%             if K > 1
%                 ante_beat{pp}(nn) = beat(pp).time(K-1);
%             else
%                 ante_beat{pp}(nn) = beat(pp).pre_time;
%             end
%         else
%             last_beat{pp}(nn) = beat(pp).pre_time;
%             ante_beat{pp}(nn) = beat(pp).pre_time;
%         end
%     end
    interp_mat{pp} = heartbeat_interpolation(algo, model, time, last_beat{pp}, ante_beat{pp});
    tmp = repmat(interp_mat(pp),model.num_sens,1);
    vec_interp_mat{pp} = blkdiag(tmp{:});

end
H = [vec_interp_mat{:}];

Y = observ';
Y = Y(:);

end

