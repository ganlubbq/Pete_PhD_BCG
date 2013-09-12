function [ new_pt ] = pf_forwardparticle( model, anc, old_pt, start_time, end_time )
%PF_FORWARDPARTICLE Create a new particle in a particle filter structure
%derived from its ancestor particle from the previous processing frame.

% People loop
for pp = 1:model.np
    
    % Get the old set of heartbeats
    old_beat = old_pt.beat(pp);
    
    % Find the latest beat which has passed out of the window
    fixed_time = old_beat.time(old_beat.time<start_time);
    if ~isempty(fixed_time)
        % If there is one, then use it as the new preceeding beat
        latest = find( old_beat.time==max(fixed_time) );
        pre_time = old_beat.time(latest);
        pre_param = old_beat.param(:,latest);
        time = old_beat.time(latest+1:end);
        param = old_beat.param(:,latest+1:end);
        pqratio = old_beat.pqratio(latest+1:end);
        if latest > 1
            ante_param = old_beat.param(:,latest-1);
        else
            ante_param = old_beat.pre_param(:);
        end 
    else
        % Otherwise, nothing's changed
        ante_param = old_beat.ante_param;
        pre_time = old_beat.pre_time;
        pre_param = old_beat.pre_param;
        time = old_beat.time;
        param = old_beat.param;
        pqratio = old_beat.pqratio;
    end
    
    new_beats(pp) = beat_init(model, pre_time, pre_param, ante_param, time, param, pqratio);
    
end

% Create particle
new_pt = part_init(model, anc, old_pt.weight, new_beats);
    
end

