%% open this folder and navigate to the processed_data folder for that fly
%
% flyinfo.xlsx should contain info about which thorimage and thorsync files
% go together
%
% 
stimtype = load_stimtype_from_xls('movie/metadata.xlsx');

tbl_flyinfo = load_sessioninfo_from_xls('flyinfo.xlsx');
analysis_steps = {'metadata', 'make_ti', 'make_frametimeinfo', 'bulkoverwrite'};

% make new empty directories for thorsync files
for i = 1:height(tbl_flyinfo)
    
    if strcmp(astep, 'bulkoverwrite')
    tiname = tbl_flyinfo.thorimagename(i);
    fname = fullfile(thorimagename, 'metadata.xlsx');
    writetable(stimtype, fname, 'Sheet', 'stimtype');
    
    %[sessioninfo_tbl, stimulus_tbl, olfactometer_tbl] =  helper_fill_xls_meta(thorimagename);
end
%%

function [sessioninfo_tbl, stimulus_tbl, olfactometer_tbl]=  helper_fill_xls_meta(thorimagename)
    fname = fullfile(thorimagename, 'metadata.xlsx');
    fprintf('Updating: \n\t %s\n', fname)
    [sessioninfo_tbl, stimulus_tbl, olfactometer_tbl] = fill_xls_meta(fname);
end

