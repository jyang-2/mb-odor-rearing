%% Run to copy Experiment.xml files from raw data folder to processed_data folders
% need to have already created flyinfo.xlsx (run
% mb_odor_rearing/src/make_excel_metadata.py)
% 

close('all');
clear, clc;

flydate = '2021-08-24';
flynum = 1;

PRJ_DIR =  "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
RAW_DIR = fullfile(PRJ_DIR, 'data', '2p', flydate, num2str(flynum));
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data', flydate, num2str(flynum));
tbl_flyinfo = load_sessioninfo_from_xls(fullfile(PROC_DIR, 'flyinfo.xlsx'));

for i = 1:height(tbl_flyinfo)
    tiname = tbl_flyinfo.thorimagename(i);
    tsname = tbl_flyinfo.thorsyncname(i);
    copyfile(fullfile(RAW_DIR, tiname, 'Experiment.xml'), fullfile(PROC_DIR, tiname, 'Experiment.xml'));
end
%%
