function [ cp_list, pf_cp ] = process_pf( algo, model, pf )
%PROCESS_PF Take a particle filter output structure and collate useful
%arrays of things

% List of all changepoints
cp_list = cell2mat(arrayfun(@(x) {cell2mat(x.cp_time)'}, pf));

% Particle changepoint sets
pf_cp = cell(algo.Nf,1);
for ii = 1:algo.Nf
    cp_set = [];
    anc = ii;
    for kk = model.K:-1:1
        if ~isempty(pf(kk).cp_time{anc})
            cp_set = [pf(kk).cp_time{anc} cp_set];
        end
        anc = pf(kk).ancestor(anc);
    end
    pf_cp{ii} = cp_set;
end


end

