function ti = ai2ti(ai)
%AI2TI given thorsync recordings, will find frametimes and other
% trial-related events
% 
%   Detailed explanation goes here


if isfield(ai, 'FrameOut')
    ai.Frame_Out = ai.FrameOut;
    ai = rmfield(ai, 'FrameOut');
end
    
[frame_times, indexes] = thorSync.generate_frame_time(ai.Frame_Out, ai.time);
ti = struct();
ti.frame_times = frame_times;


[rLT, fLT, fT] = thorSync.get_frameout_times(ai.Frame_Out', ai.time');
ti.frameouts.rT = rLT;   % rise time
ti.frameouts.fT = fLT;   % fall time of frame outs
ti.frameouts.T = (rLT + fLT)/2;


% get timing information about, and number of, blocks and stim
% acquisition block info
[W,ict,fct,~] = pulsewidth(ai.scopePin, ai.time);
ti.num_blocks = length(W);                 % # of blocks
ti.block_len = W;
ti.block_ict = ict;
ti.block_fct = fct;

[W,ict,fct,~] = pulsewidth(ai.olfDispPin, ai.time);  
ti.num_stim = length(W);                  % # of odor pulses
ti.stim_len = W;
ti.stim_ict = ict;
ti.stim_fct = fct;

ti.timepoints = length(ti.frameouts.rT);
%% make sure all fields are rows

ti = make_fields_rows(ti);
ti = ti_frame_corrector(ti);
ti.spt = ti.num_stim / ti.num_blocks;               % stimuli per trial

end

function scopeIdx = get_scope_idx(frame_times, block_ict, block_fct)
    
end