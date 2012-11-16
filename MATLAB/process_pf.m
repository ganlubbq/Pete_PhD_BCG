function [ cp_list, pf_cp, pf_p, pf_a, pf_clut, rb_est, clut_indic, reconstructed ] = process_pf( algo, model, time, pf )
%PROCESS_PF Take a particle filter output structure and collate useful
%arrays of things

% List of all changepoints
% cp_list = cell2mat(arrayfun(@(x) {cell2mat(x.cp_time)'}, pf));
cp_list = unique(cat(2,pf.cp_time));

% Particle changepoint sets
pf_cp = cell(algo.Nf,1);
pf_p = cell(algo.Nf,1);
pf_a = cell(algo.Nf,1);
pf_clut = zeros(algo.Nf, model.K);
for ii = 1:algo.Nf
    anc = ii;
    for kk = model.K:-1:1
        pf_clut(ii,kk) = pf(kk).clut_indic(anc);
        if ~ismember( pf(kk).cp_time(anc), pf_cp{ii} )
            pf_cp{ii} = [pf(kk).cp_time(anc) pf_cp{ii}];
            pf_p{ii} = [pf(kk).cp_param(1,anc) pf_p{ii}];
            pf_a{ii} = [pf(kk).cp_param(2,anc) pf_a{ii}];
        end
        anc = pf(kk).ancestor(anc);
    end
end

% Reconstruct estimated signal
reconstructed = zeros(algo.Nf,model.K);
for kk = 1:model.K
    weight = pf(kk).weight;
    weight = weight - max(weight);
    weight = exp(weight);
    weight = weight/sum(weight);
    for ii = 1:algo.Nf
        cp_time = pf(kk).cp_time(ii);
        cp_a = pf(kk).cp_param(2,ii);
        rb_mn = pf(kk).rb_mn(:,ii);
        interp_vector = heartbeat_interpolation(algo, model, time(kk), cp_time);
        reconstructed(ii,kk) = weight(ii)*cp_a*interp_vector*rb_mn;
        
    end
end
reconstructed = sum(reconstructed);

%% Fixed rate estimates
rb_est = zeros(model.dw, model.K);
clut_indic = zeros(1, model.K);
for kk = 1:model.K
    weight = pf(kk).weight - max(pf(kk).weight);
    lin_weight = exp(weight);
    lin_weight = lin_weight/sum(lin_weight);
    rb_est(:,kk) = (lin_weight * pf(kk).rb_mn')';
    clut_indic(1,kk) = lin_weight * pf(kk).clut_indic';
end

end

