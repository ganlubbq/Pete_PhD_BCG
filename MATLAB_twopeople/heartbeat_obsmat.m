function [ H, Y ] = heartbeat_obsmat( algo, model, time, observ, beat )
%HEARTBEAT_OBSMAT Build the linear interpolation matrix and vectorise the
%corresponding observations

last_beat = cell(model.np,1);
interp_mat = cell(model.np,1);
vec_interp_mat = cell(model.np,1);
for pp = 1:model.np
    last_beat{pp} = zeros(1,0);
    for nn = 1:length(time)
        [tauK, K] = max(beat(pp).time( beat(pp).time<=time(nn) ));
        if ~isempty(K)
            last_beat{pp}(nn) = tauK;
        else
            last_beat{pp}(nn) = beat(pp).pre_time;
        end
    end
    interp_mat{pp} = heartbeat_interpolation(algo, model, time, last_beat{pp});
    tmp = repmat(interp_mat(pp),model.num_sens,1);
    vec_interp_mat{pp} = blkdiag(tmp{:});
end
H = [vec_interp_mat{:}];

Y = observ';
Y = Y(:);

end

