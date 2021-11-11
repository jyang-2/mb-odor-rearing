function AI = splice_thorsync_files(filepaths)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nFiles = numel(filepaths);


for i = 1:nFiles    
    fname = filepaths(i);
    ai = load(fname);
    
    if i == 1
        AI = ai;
        flds = fields(ai);
        nFields = numel(flds);
    else
        ai.time = ai.time + AI.time(end) + .001;
        ai.FrameCounter = ai.FrameCounter + AI.FrameCounter(end);
        
        for j = 1:nFields
           f = flds{j};
           AI.(f) = [AI.(f); ai.(f)];
        end
    end
end    

end

