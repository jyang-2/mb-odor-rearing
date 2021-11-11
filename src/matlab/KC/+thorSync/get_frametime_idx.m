function idx = get_frametime_idx(frame_times, ict, fct)
%GET_SCOPEP Summary of this function goes here
%   Detailed explanation goes here


nBlocks = numel(ict);

idx = zeros(size(frame_times));
for i = 1:nBlocks
    rT = ict(i);
    fT = fct(i);
   idx(frame_times>=rT & frame_times<fT) = i;
end

end

