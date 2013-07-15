function [ pt ] = part_init( model, ancestor, weight, beats )
%PART_INIT Initialise a particle

pt = struct(...
    'ancestor', ancestor, ...
    'weight', weight, ...
    'beat', beats );

end

