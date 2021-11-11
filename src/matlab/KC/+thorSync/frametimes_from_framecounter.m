function frametimes = frametimes_from_framecounter(frame_counter, time)
%FRAMETIMES_FROM_FRAMECOUNTER Summary of this function goes here
%   Detailed explanation goes here

[C,ia,~] = unique(frame_counter,'first');

idx = ia(C>0);
frametimes = time(idx);
end

