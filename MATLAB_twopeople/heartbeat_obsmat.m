function [ H, Y, dH ] = heartbeat_obsmat( algo, model, time, observ, beat )
%HEARTBEAT_OBSMAT Build the linear interpolation matrix and vectorise the
%corresponding observations

last_beat = cell(model.np,1);
interp_mat = cell(model.np,1);
vec_interp_mat = cell(model.np,1);
deriv_interp_mat = cell(model.np,1);
deriv_vec_interp_mat = cell(model.np,1);
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

    if nargout > 2
        deriv_interp_mat{pp} = heartbeat_interpolationderiv(algo, model, time, last_beat{pp});
        tmp = repmat(deriv_interp_mat(pp),model.num_sens,1);
        deriv_vec_interp_mat{pp} = blkdiag(tmp{:});
    end

end
H = [vec_interp_mat{:}];
if nargout > 2
    dH = [deriv_vec_interp_mat{:}];
else
    dH = [];
end

Y = observ';
Y = Y(:);

end

