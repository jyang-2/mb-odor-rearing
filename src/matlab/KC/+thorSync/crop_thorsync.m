function ai_crop = crop_thorsync(ai, start_idx, end_idx)
%CROP_THORSYNC crop thorsync signal to desired window of time
% inputs:
%   -ai : loaded thorSync traces (struct)
%   - start_idx: starting index of time window
%   - end_idx: stopping index of time window
% outputs
%   -ai_cropped : thorsync traces cropped to desired time/index window

f = fields(ai);
ai_crop = struct();

siz = size(ai.time);

for i = 1:numel(f)
    x = ai.(f{i});
    if size(x)==siz
        ai_crop.(f{i}) = x(start_idx:end_idx);    
    else
        ai_crop.(f{i}) = x;
    end
    
end


end

