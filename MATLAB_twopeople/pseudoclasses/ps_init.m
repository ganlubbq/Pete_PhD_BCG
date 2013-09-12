function [ ps ] = ps_init( model, N )
%PS_INIT Initialise particle smoother structure for heartbeat inference

beat = beat_init(model, [], [], [], [], [], []);
beats = repmat(beat, 1, model.np);

ps = repmat( struct(...
    'beat', beats ),...
    N, 1);

end

