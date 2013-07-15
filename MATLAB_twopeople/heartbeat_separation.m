function [mu, Sigma] = heartbeat_separation( display, algo, model, time, observ, beat )
%HEARTBEAT_SEPARATION Separate heartbeats with a Gibbs Sampler

% Set up interpolation matrixes
last_beat = cell(model.np,1);
interp_mat = cell(model.np,1);
vec_interp_mat = cell(model.np,1);
for pp = 1:model.np
    for nn = 1:length(time)
        [tauK, K] = max(beat(pp).time( beat(pp).time<=time(nn) ));
        if ~isempty(K)
            last_beat{pp}(nn) = tauK;
        else
            last_beat{pp}(nn) = -inf;
        end
    end
    interp_mat{pp} = heartbeat_interpolation(algo, model, time, last_beat{pp});
    tmp = repmat(interp_mat(pp),model.num_sens,1);
    vec_interp_mat{pp} = blkdiag(tmp{:});
end

% Vectorise model
H = [vec_interp_mat{:}];
Y = observ';
Y = Y(:);
Pw = model.w_prior_vr;
Pw = blkdiag(Pw, Pw, Pw, Pw, Pw, Pw, Pw, Pw);

% Posterior
Sigma = inv( inv(Pw) + H'*H/model.y_obs_vr );
mu = Sigma*H'*Y/model.y_obs_vr;

w_av{1} = reshape(mu(1:4*model.dw), model.dw, 4)';
w_av{2} = reshape(mu(4*model.dw+1:end), model.dw, 4)';

if display.plot_during
    figure, plot(w_av{1}');
    figure, plot(w_av{2}');
    
    x = zeros(model.num_sens,length(time));
    for pp = 1:model.np
        x = x + w_av{pp}*interp_mat{pp}';
    end
    
    figure, plot(time, x'); ylim([-1 1]);
    figure, plot(time, observ); ylim([-1 1]);
    
    err = observ-x;
    figure, plot(time, err'); ylim([-1 1]);
    figure,plot(sqrt(diag(err'*err)))
end

end

