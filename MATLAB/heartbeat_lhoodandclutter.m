function [ pf_win ] = heartbeat_lhoodandclutter( algo, model, L, time, observ, pre_cp_time, pre_cp_param, pre_rb_mn, pre_rb_vr, pre_clut_indic, win_cp_time, win_cp_param, flag_no_clut )
%HEARTBEAT_LHOODANDCLUTTER Kalman filter heartbeat template through the
%window, sampling clutter indicators and evaluating likelihood as we go.

last_cp_time = pre_cp_time;
last_cp_param = pre_cp_param;
if isempty(last_cp_time)
    last_cp_time = -inf;
end
rb_mn = pre_rb_mn;
rb_vr = pre_rb_vr;
clut_indic = pre_clut_indic;

for ll = 1:L
    
    % Update changepoint if we've past one
    if ~isempty(win_cp_time) && (time(ll) > win_cp_time) && (last_cp_time < win_cp_time)
        if (last_cp_param(2) == 1)&&(win_cp_param(2) == 0)
            rb_vr = rb_vr + model.w_trans_vr;
        else
            rb_vr = rb_vr + model.w_trans_vr;
        end
        last_cp_time = win_cp_time;
        last_cp_param = win_cp_param;
    end
            
    % Interpolation and Kalman filtering
    H = heartbeat_interpolation(algo, model, time(ll), last_cp_time);
            
    if isempty(last_cp_param) || (last_cp_param(2) == 0)
                
        [noclut_rb_mn, noclut_rb_vr, ~, s_mn, s_vr] = kf_update(rb_mn, rb_vr, observ(ll), H, model.y_obs_vr);
        noclut_lhood = loggausspdf(observ(ll), s_mn, s_vr);
                
        % Clutter sampling
        clut_rb_mn = rb_mn; clut_rb_vr = rb_vr;
        clut_lhood = loggausspdf(observ(ll), 0, model.y_clut_vr);
%         clut_lhood = log(tpdf(observ(ll)/sqrt(model.y_clut_vr), 1)/sqrt(model.y_clut_vr));
        if flag_no_clut
            clut_prior = log([0; 1]);
        else
            clut_prior = log(model.clut_trans(:,clut_indic+1));
        end
        clut_prob = clut_prior + [clut_lhood; noclut_lhood];
        lhood = logsumexp(clut_prob);
        clut_prob = clut_prob - lhood;
        clut_indic = log(rand)<clut_prob(1);
        if clut_indic
            rb_mn = clut_rb_mn;
            rb_vr = clut_rb_vr;
        else
            rb_mn = noclut_rb_mn;
            rb_vr = noclut_rb_vr;
        end
                
    elseif last_cp_param(2) == 1
        
        lhood = log(tpdf(observ(ll)/sqrt(model.y_dstb_vr), 1)/sqrt(model.y_dstb_vr));
        s_mn = H * rb_mn;
        s_vr = H * rb_vr * H' + model.y_obs_vr;
        clut_indic = 0;
        
    end

    % Store everything
    pf_win.rb_mn(:,ll) = rb_mn;
    pf_win.rb_vr(:,:,ll) = rb_vr;
    pf_win.obslhood(ll) = lhood;
    pf_win.signal_mn(ll) = s_mn;
    pf_win.signal_vr(ll) = s_vr;
    pf_win.clut(ll) = clut_indic;
    
end



end

