function frame_times = new_frame_times(ai, averageNum)
%NEW_FRAME_TIMES Computes frame times, using the field 'FrameCounter'
% inputs: 
%   ai
%    - time
%    - FrameCounter
%   averageNum
%
% outputs:
%   frame_times
arguments
    ai struct
    averageNum (1,1) {mustBeNumeric} = 1
end
    
new_frame = diff([0; ai.FrameCounter]);
frame_times = ai.time(find(new_frame));
frame_times = frame_times(1:averageNum:end);

frame_times = frame_times(:)';
end

