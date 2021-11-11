function ai_ds = downsample_ai(ai, k)
%DOWNSAMPLE_AI downsample all lines in ai struct by a factor of k
% inputs:
%   -ai : loaded thorSync traces
%   - k : downsampling factor
% outputs
%   -ai_ds: downsampled ai struct

f = fields(ai);
ai_ds = struct();

for i = 1:numel(f)
    x = downsample(ai.(f{i}), k);
    ai_ds.(f{i}) = x;
end

end

