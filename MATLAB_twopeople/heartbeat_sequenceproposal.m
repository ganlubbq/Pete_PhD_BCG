function [beat, sf_ratio] = heartbeat_sequenceproposal(algo, model, beat, time, observ, wf_mn, wf_vr)

% We're going to sequentially propose a new beat and then optimise it for
% each person until we overun the window.

% Throw away old heartbeats
for pp = 1:model.np
    beat(pp).time = [];
    beat(pp).param = [];
    beat(pp).pqratio = [];
end

% Sample a new pre_param, to refresh the diversity
for pp = 1:model.np
%     if ~isempty(beat(pp).ante_param)
        beat(pp).pre_param = heartbeat_paramtrans(model, beat(pp).ante_param, []);
%     end
end

% Start and end times of window
start_time = time(1);
end_time = time(end);

% Which person and beat indices
if (model.np==1)
    p_idx = 1;
    b_idx = 0;
else
    p_idx = [1 2];
    b_idx = [0 0];
end

% Get the latest beat time and parameters
latest_time = [beat.pre_time];
latest_param = [beat.pre_param];

% Initialise arrays
previous_time = zeros(size(latest_time));
previous_param = zeros(size(latest_param));
interm_time = zeros(size(latest_time));
interm_param = zeros(size(latest_param));

interm_sf = zeros(1,model.np);
new_interm_sf = zeros(1,model.np);
for pp = 1:model.np
    interm_sf(pp) = heartbeat_survivorfunc(model, latest_time(pp), latest_param(pp), end_time);
end
final_sf = interm_sf;

% Loop
while true
    
    % Sample intermediate times
    for pp = 1:model.np
        if ~isnan(b_idx(pp))
            b_idx(pp) = b_idx(pp) + 1;
            lower_lim = gamcdf(start_time-(latest_time(pp)+latest_param(pp)), model.tau_trans_shape, model.tau_trans_scale);
            u = unifrnd(lower_lim, 1);
            interm_time(pp) = latest_time(pp) + latest_param(pp) + gaminv(u, model.tau_trans_shape, model.tau_trans_scale);
            interm_param(pp) = heartbeat_paramtrans(model, latest_param(pp), []);
            new_interm_sf(pp) = heartbeat_survivorfunc(model, interm_time(pp), interm_param(pp), end_time);
        end
    end
    
    % Make an intermediate beat structure
    interm_beat = beat;
    
    % See whose left in
    for pp = 1:model.np
        if ~isnan(b_idx(pp))
            % Have we overshot the window?
            if interm_time(pp) > end_time
                % Overshot - remove from the sampling list
                b_idx(pp) = NaN;
                p_idx(p_idx==pp) = [];
            else
                interm_beat(pp).time = [interm_beat(pp).time interm_time(pp)];
                interm_beat(pp).param = [interm_beat(pp).param interm_param(pp)];
                interm_sf(pp) = new_interm_sf(pp);
            end
        end
    end
    
    % Break if everyone's done
    if isempty(p_idx)
        break;
    end
    
    % Optimise times
    [optimal_time, hess, LB, UB] = optimise_beat_times(algo, model, interm_beat, p_idx, b_idx(~isnan(b_idx)), time, observ, wf_mn, wf_vr);
    
    % Propose times
    previous_time(p_idx) = latest_time(p_idx);
    previous_param(p_idx) = latest_param(p_idx);
    tmp = diag(inv(full(hess)));
    ppsl_vr = NaN(model.np,1);
    ppsl_vr(p_idx) = tmp;
    ppsl_vr(ppsl_vr<0) = 0.001^2;
    for pp = 1:model.np
        if ~isnan(b_idx(pp))
            
            % Propose a new time
            latest_time(pp) = Inf;
            while (latest_time(pp)<LB(pp))||(latest_time(pp)>UB(pp))
                latest_time(pp) = mvnrnd(optimal_time(pp), ppsl_vr(pp));
            end
            latest_param(pp) = interm_param(pp);
            final_sf(pp) = heartbeat_survivorfunc(model, latest_time(pp), latest_param(pp), end_time);
            
            % Calculate probabilities
            ppsl_dens = loggausspdf(latest_time(pp), optimal_time(pp), ppsl_vr(pp));
            period_offset = latest_time(pp)-(previous_time(pp)+previous_param(pp));
            trans_dens = log(gampdf(period_offset, model.tau_trans_shape, model.tau_trans_scale));
            
            % Store
            beat(pp).time = [beat(pp).time, latest_time(pp)];
            beat(pp).param = [beat(pp).param, latest_param(pp)];
            beat(pp).pqratio = [beat(pp).pqratio, trans_dens - ppsl_dens];
            
        end
    end
    
end

sf_ratio = final_sf - interm_sf;

end



function [tau, hess, LB, UB] = optimise_beat_times(algo, model, beat, p_idx, b_idx, time, observ, wf_mn, wf_vr)

% p_idx is the index(es) of the person for whom a heartbeat needs optimising
% b_idx is the index(es) of the heartbeat
% They should have the same length

tol = 0.001;
tolf = 0.1;
opt_alg = 'trust-region-reflective' ;

if any([beat.pre_time]<0)
    half_width = 1;
else
    half_width = 0.15;
end

LB = NaN(model.np,1);
UB = NaN(model.np,1);
tau = NaN(model.np,1);

if isempty(p_idx)
    % No search required
    
    hess = [];

elseif length(p_idx) == 1
    % Line search
    
    tau0 = beat(p_idx).time(b_idx);
    
    % Minimise negative log-likelihood
    if b_idx > 1
        tau_offset(p_idx) = beat(p_idx).time(b_idx-1) + beat(p_idx).param(b_idx-1);
    else
        tau_offset(p_idx) = beat(p_idx).pre_time + beat(p_idx).pre_param;
    end
    h_of = @(tau) log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau, tau_offset);
    LB(p_idx) = max(tau0-half_width, time(1));
    LB(p_idx) = max(LB(p_idx), tau_offset(p_idx));
    UB(p_idx) = min(tau0+half_width, time(end));
%     options = optimset('Display','notify-detailed', 'TolX',tol);
%     tau = fminbnd(h_of, LB, UB, options);
    options = optimset('GradObj','on', 'Display','notify-detailed', 'TolFun',tolf, 'TolX',tol,'algorithm',opt_alg);
    [tau(p_idx), ~, ~, ~, ~, ~, hess] = fmincon(h_of, tau0, [], [], [], [], LB(p_idx), UB(p_idx), [], options);
    
    % Approximate Hessian?
%     if nargout > 1
%         dt = 0.0001;
%         [~, grad0] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau);
%         [~, grad1] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau+dt);
%         hess = (grad1-grad0)/dt;
%     else
%         hess = [];
%     end
    
else
    % 2D search

    for bb = 1:length(p_idx)
        tau0(bb) = beat(p_idx(bb)).time(b_idx(bb));
    end
    tau0 = tau0';
    
    % Minimise negative log-likelihood
    tau_offset = zeros(1,length(p_idx));
    for pp = 1:2
        if b_idx(pp) > 1
            tau_offset(pp) = beat(p_idx(pp)).time(b_idx(pp)-1) + beat(p_idx(pp)).param(b_idx(pp)-1);
        else
            tau_offset(pp) = beat(p_idx(pp)).pre_time + beat(p_idx(pp)).pre_param;
        end
    end
    h_of = @(tau) log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau, tau_offset);
    LB = max(tau0-half_width, time(1));
    for pp = 1:2
        LB(pp) = max(LB(pp), tau_offset(pp));
    end
    UB = min(tau0+half_width, time(end));
    options = optimset('GradObj','on', 'Display','notify-detailed', 'TolX', tol, 'algorithm',opt_alg);
    [tau, ~, ~, ~, ~, ~, hess] = fmincon(h_of, tau0, [], [], [], [], LB, UB, [], options);
    
%     % Hessian?
%     hess = [];
    
end

end



function [func, grad, hess] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau, tau_min)

for bb = 1:length(p_idx)
    beat(p_idx(bb)).time(b_idx(bb)) = tau(bb);
end

[H, Y] = heartbeat_obsmat(algo, model, time, observ, beat);
R = model.y_obs_vr*eye(length(Y));
% S = R + H*wf_vr*H';

% y_minus_Hm = Y-H*wf_mn;
% G = H*wf_vr;
% S = R + G*H';
% nu = S\y_minus_Hm;

y_minus_Hm = Y-mtimesx(H,'n',wf_mn,'n','SPEED');

if tau < 3

    G = mtimesx(H,'n',wf_vr,'n','SPEED');
    S = R + mtimesx(G,'n',H,'t','SPEED');
    % U = S\G;
    nu = S\y_minus_Hm;
    
    % func = -loggausspdf(Y, H*wf_mn, S);
    % func = -loggausspdf(y_minus_Hm, 0, S);
    
    % log_2(pi) is 1.83787706640935;
    func = 0.5*( y_minus_Hm'*nu + 1.83787706640935*length(Y) +  sum(log(diag(chol(S)))) );
    
    for bb = 1:length(p_idx)
        func = func + (1-model.tau_trans_shape)*log(tau(bb)-tau_min(bb)) + (tau(bb)-tau_min(bb))/model.tau_trans_scale;
    end
    
    if nargout > 1
        
        grad = zeros(length(p_idx),1);
        hess = zeros(length(p_idx));
        
        %     %%%%%
        %     M = S\G;
        %     zeta = G'*nu;
        %     for bb = 1:length(p_idx)
        %         dH = heartbeat_obsmatderiv(algo, model, time, beat, p_idx(bb), b_idx(bb));
        %         grad(bb) = -(  wf_mn'*dH'*nu + zeta'*dH'*nu - trace(dH'*M)  );
        %     end
        %     %%%%%
        
        for bb = 1:length(p_idx)
            dH = heartbeat_obsmatderiv(algo, model, time, beat, p_idx(bb), b_idx(bb));
            
            %         T = G*dH';
            %         M = T/S;
            T = mtimesx(G,'n',dH,'t','SPEED');
            M = zeros(size(T));
            L = min(algo.L, length(time));
            for ii = 1:model.num_sens
                idx_rng = (L*(ii-1)+1):(L*ii);
                M(idx_rng,idx_rng) = T(idx_rng,idx_rng)/S(idx_rng,idx_rng);
            end
            grad(bb) = -(  wf_mn'*dH'*nu + nu'*M*y_minus_Hm - trace(M)  );
            
            grad(bb) = grad(bb) + (1-model.tau_trans_shape)/(tau(bb)-tau_min(bb)) + 1/model.tau_trans_scale;
            
            %         V = U*dH';
            %         grad(bb) = -(  wf_mn'*dH'*nu + y_minus_Hm'*V*nu - trace(V)  );
            
            %         invSigma = inv(wf_vr)+H'*(R\H);
            %         mu = invSigma\( H'*(R\Y) + wf_vr\wf_mn );
            %         grad(bb) = -( mu'*dH'*(R\(Y-H*mu)) - trace( H'*(R\dH)/invSigma ) );
            
            %         if nargout > 2
            %
            %             d2H =
            %
            %             h_Omega = T+T';
            %             h_zeta = h_Omega*nu;
            %             h_Theta = G*d2H;
            %             h_xi = dH'*nu;
            %             h_mtilde = dH'*m;
            %             h_Gamma = h_Omega*V;
            %             h_Lambda = dH*wf_vr*dH';
            %             h_Xi = (h_Lambda + h_Theta - h_Gamma)/S;
            %
            %             hess(bb) = - h_zeta'*(S\h_zeta) ...
            %                        - h_xi'*wf_vr*h_xi ...
            %                        - nu'*h_Theta*nu ...
            %                        + wf_mn'*d2H*nu ...
            %                        + h_mtilde'*(S\h_mtilde) ...
            %                        + trace(h_Xi);
            %
            %         end
            
        end
        
    end

else
    
    
    nu = R\y_minus_Hm;
    % log_2(pi) is 1.83787706640935;
    func = 0.5*( y_minus_Hm'*nu + 1.83787706640935*length(Y) +  2*sum(log(diag(chol(R)))) );
    
    if nargout > 1
        grad = zeros(length(p_idx),1);
        for bb = 1:length(p_idx)
            dH = heartbeat_obsmatderiv(algo, model, time, beat, p_idx(bb), b_idx(bb));
            grad(bb) = -wf_mn'*dH'*nu;
        end
    end

end

end
