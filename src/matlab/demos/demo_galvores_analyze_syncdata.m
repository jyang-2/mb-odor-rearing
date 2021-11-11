close('all'), clear, clc;

flydate = '2021-07-27';
flynum = 1;

%PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";



% set directory paths
if ispc
    listing = dir(fullfile("G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data", flydate, flynum, "**\flyinfo.xlsx"));
    PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing";
    DATA_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\2p";
    PROC_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data";
elseif isunix
    listing = dir(fullfile("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-06-26",...
                            "**/flyinfo.xlsx"));
    DATA_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/2p";
    PROC_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data";
end

RAW_DIR = fullfile(PRJ_DIR, 'data', '2p', flydate, num2str(flynum));
FLY_DIR = fullfile(PRJ_DIR, 'data', 'processed_data', flydate, num2str(flynum));

tbl_flyinfo = load_sessioninfo_from_xls(fullfile(FLY_DIR, 'flyinfo.xlsx'));

if ismember('preprocess', tbl_flyinfo.Properties.VariableNames)
   tbl_flyinfo = tbl_flyinfo((tbl_flyinfo.preprocess==1),:);
end

disp(tbl_flyinfo)
%% Initial frame time calculation
% saves ti0.mat and frametimeinfo0.csv
load_syncdata_type = 'h5'; % {'h5', 'mat'}

% options for thorSync.make_frametimeinfo()
opts = struct('baseline_win', 15, 'peakresp_win', 15, 'trial_win', 20);

for i = 1:height(tbl_flyinfo)    
    tiname = tbl_flyinfo.thorimagename(i);
    tsname = tbl_flyinfo.thorsyncname(i);
    
    if ~ismissing(tsname)
        PROC_DIR = fullfile(FLY_DIR, tiname);
        disp(PROC_DIR);
        if strcmp(load_syncdata_type, 'mat') % load downsampled .mat file
            ai = load(fullfile(PROC_DIR, 'syncdata.mat'));
        elseif strcmp(load_syncdata_type, 'h5') % load raw .h5 file
            h5filepath = fullfile(RAW_DIR, tsname, 'Episode001.h5');
            ai = thorSync.load_sync_episode(char(h5filepath));
        end
       
        ti = thorSync.ai2ti(ai);
        T = thorSync.make_frametimeinfo(ti.frame_times, ti, opts);
       
        fprintf('\n%s: %s\n', tiname, tsname)
        fprintf('\tSaving ti0.mat\n')
        save(fullfile(PROC_DIR, 'ti0.mat'), '-struct', 'ti');
        
        fprintf('\tSaving frametimeinfo0.csv\n');
        writetable(T, fullfile(PROC_DIR, 'frametimeinfo0.csv'));        
    end
end

%% correct frametimes using frame averaging info from Experiment.xml
% saves ti.mat, frametimeinfo.csv
for i = 1:height(tbl_flyinfo)   
    tiname = tbl_flyinfo.thorimagename(i);
    tsname = tbl_flyinfo.thorsyncname(i);
    
    if ~ismissing(tsname)
        PROC_DIR = fullfile(FLY_DIR, tiname);
        disp(PROC_DIR);
        
        ti0 = load(fullfile(PROC_DIR, 'ti0.mat'));
        T0 = readtable(fullfile(PROC_DIR, 'frametimeinfo0.csv'));        

        % parse Experiment.xml info
        xmlpath = fullfile(PROC_DIR, 'Experiment.xml');
        exml = readstruct(xmlpath);
        s = parse_experiment_xml(exml); % s : experiment xml info
        
        % # TIMEPOINTS WRONG
        if height(T0) ~= s.num_timepoints
            % single plane
            if s.fast_z == 0 
                % downsample frametimeinfo (table T) w/ frame averaging #
                if s.average_num>1
                    tbls = cell(1, ti0.num_blocks);
                    for iblock = 1:ti0.num_blocks
                       tbl0 = T0(T0.scope_pulse_idx==iblock,:);
                       tbls{iblock} = tbl0(1:s.average_num:end,:);
                    end
                    T = vertcat(tbls{:});
                    T.frame_num = (1:height(T))';
                elseif s.average_num == 1
                    if (height(T0) ~= s.num_timepoints)
                        % trim # of frames, using rising/falling time to check against ai.scopePin
                        good_idx = frames_within_scopepulse(ti0.frameouts.rT, ti0.frameouts.fT, ti0.block_ict, ti0.block_fct);
                        T = T0(good_idx, :);
                        T.frame_num = (1:height(T))';
                    end
                end
                
            % multiplane, volumetric, fast-z mode
            elseif s.fast_z > 0 
                tbls = cell(1, ti0.num_blocks);
                for iblock = 1:ti0.num_blocks
                    tbl0 = T0(T0.scope_pulse_idx == iblock, :);
                    num_block_frames = floor(height(tbl0)/s.total_planes) * s.total_planes;
                    
                    tbl = tbl0(1:num_block_frames, :);
                    tbl.plane = mod((1:height(tbl))' - 1, s.total_planes);
                    tbls{iblock} = tbl;
                end
                T = vertcat(tbls{:});
                T.frame_num = (1:height(T))';
            end       
        % TIMEPOINTS CORRECT - frametimeinfo0, ti0 frame times are correct
        else    
            T = T0; 
        end        
        
        % if # of timepoints in T is correct, make ti
        if height(T) == s.num_frames
            ti = ti0;
            ti.frame_times = T.frame_times';
            ti.timepoints = length(ti.frame_times);
            
            fprintf('\n%s: %s\n', tiname, tsname)
            fprintf('\tSaving ti.mat\n')
            save(fullfile(PROC_DIR, 'ti.mat'), '-struct', 'ti');

            fprintf('\tSaving frametimeinfo.csv\n');
            writetable(T, fullfile(PROC_DIR, 'frametimeinfo.csv'));        
        else
            fprintf('\n%s: %s\n', tiname, tsname)
            fprintf('\tUnable to resolve timepoint/frame # conflict.\n')
        end
    end
end


%% check that # of timepoints is correct

tbl_timepoints = tbl_flyinfo(~ismissing(tbl_flyinfo.thorsyncname), {'date', 'fly', 'thorimagename', 'thorsyncname'});

for i = 1:height(tbl_timepoints)    
    tiname = tbl_timepoints.thorimagename(i);
    tsname = tbl_timepoints.thorsyncname(i);
    
    
    PROC_DIR = fullfile(FLY_DIR, tiname);
    disp(PROC_DIR);

    ti = load(fullfile(PROC_DIR, 'ti.mat'));
    T = readtable(fullfile(PROC_DIR, 'frametimeinfo.csv'));        

    xmlpath = fullfile(PROC_DIR, 'Experiment.xml');
    exml = readstruct(xmlpath);
    s = parse_experiment_xml(exml); % s : experiment xml info
    
    tbl_timepoints.fast_z(i) = s.fast_z;
    tbl_timepoints.xml_timepoints(i) = s.num_timepoints;
    tbl_timepoints.ti_timepoints(i) = numel(ti.frame_times);
    tbl_timepoints.T_timepoints(i) = height(T);
    tbl_timepoints.xml_frames(i) = s.num_frames;
    
    if s.fast_z == 0
        tbl_timepoints.correct(i) =  (tbl_timepoints.xml_timepoints(i) == tbl_timepoints.ti_timepoints(i)) & ...
                                        (tbl_timepoints.xml_timepoints(i) == tbl_timepoints.T_timepoints(i));
    elseif s.fast_z == 1
        tbl_timepoints.correct(i) =  (tbl_timepoints.xml_frames(i) == tbl_timepoints.ti_timepoints(i)) & ...
                                        (tbl_timepoints.xml_frames(i) == tbl_timepoints.T_timepoints(i));
    end
    
end

%tbl_timepoints.correct = (tbl_timepoints.xml_timepoints == tbl_timepoints.ti_timepoints) & (tbl_timepoints.xml_timepoints == tbl_timepoints.T_timepoints);
disp(tbl_timepoints);

