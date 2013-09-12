function [ pf ] = pf_init( model, N )
%PF_INIT Initialise particle filter structure for heartbeat inference

beat = beat_init(model, [], [], [], [], [], []);
beats = repmat(beat, 1, model.np);

pt = part_init(model, 0, 0, beats);
pf = repmat(pt, N, 1);

end

