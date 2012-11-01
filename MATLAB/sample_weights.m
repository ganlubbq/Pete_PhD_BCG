function [ ancestor ] = sample_weights( algo, weight, N )
%SAMPLE_WEIGHTS Samples an array of N ancestor indexes from a set of
%weights

% Row vectors only here, please
weight = weight(:)';

% Convert weights to linear domain and normalise
weight = weight - max(weight);
weight = exp(weight);
weight = weight/sum(weight);

% Catch NaNs
assert(~any(isnan(weight)));

% Create bin boundaries
edges = min([0 cumsum(weight)],1);
edges(end) = 1;

% Sample indexes
if algo.resam_type == 1
    % Multinomial
    idx = rand(N,1);
elseif algo.resam_type == 2
    % Systematic
    idx = (1/N)*ones(N, 1);
    idx(1) = rand/N;
    idx = cumsum(idx);
end

% Draw particles
[~, ancestor] = histc(idx,edges);
ancestor = ancestor';

end

