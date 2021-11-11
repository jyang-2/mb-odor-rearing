close all, clear, clc

flydate = '2021-06-14';
flynum = 1;

PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
RAW_DIR = fullfile(PRJ_DIR, 'data', '2p', flydate, num2str(flynum));
FLY_DIR = fullfile(PRJ_DIR, 'data', 'processed_data', flydate, num2str(flynum));

%RAW_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/2p/2021-06-26/2";
%FLY_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-06-26/2";

tbl_flyinfo = load_sessioninfo_from_xls(fullfile(FLY_DIR, 'flyinfo.xlsx'));
disp(tbl_flyinfo);
%%
%disp(tbl_flyinfo)
%tbl_flyinfo = tbl_flyinfo(6:end,:);


% make new empty directories for thorsync files
% for i = 1:5
%     tiname = tbl_flyinfo.thorimagename(i);
%     
%    if ~isfolder(tiname)
%        disp(i)
%        mkdir(tiname)
%    end
% end
%%
%close('all'), clc
baseline_win = 10;
peakresp_win = 5;
trial_win = 20;

for i = 4:height(tbl_flyinfo)
    tiname = tbl_flyinfo.thorimagename(i);
    tsname = tbl_flyinfo.thorsyncname(i);
    
%     tsnames = strtrim(strsplit(tbl_flyinfo.thorsyncname(i), ','));
%     for j = 1:length(tsnames)
%         tsname = tsnames(j);
% 
%         h5filepath = sprintf('%s.mat', tsname);
%         if tbl_flyinfo.insubdir(i)==0
%             h5filepath = fullfile('thorSync', h5filepath);
%         elseif tbl_flyinfo.insubdir(i)==1
%             h5filepath = fullfile(tiname, h5filepath);
%         end
%         disp(h5filepath)t
        
    if ~ismissing(tsname)
        h5filepath = fullfile(RAW_DIR, tsname, 'Episode001.h5');
        ai = thorSync.load_sync_episode(h5filepath);
        % load and downsample file    
        k = 1; % downsample by a factor of 30  
        ai = thorSync.downsample_ai(ai, k);
        
        fig = plot_thorsync_lines(ai);

        savename = sprintf('fig_thorsync_lines.png');
        exportgraphics(fig, fullfile(FLY_DIR, tiname, savename));
        
        savename = sprintf('syncdata.mat');
        save(fullfile(FLY_DIR, tiname, savename), '-struct',  'ai');        
    end
        

    
    
    
%     
%     
%     
    
%     savename = fullfile(tiname, 'fig_thorsync_lines'
%     fig = plot_thorsync_lines(ai);
% 
%     exportgraphics(fig, fullfile(tiname, 'fig_thorsync_lines.png'));
%     save(fullfile(tiname, 'syncdata.mat'), '-struct',  'ai');
%    
end
%%

function fig = plot_thorsync_lines(ai)
    if isfield(ai, 'FrameOut')
        ai.Frame_Out = ai.FrameOut;
        ai = rmfield(ai, 'FrameOut');
    end
    fig = figure;
    plot(ai.time, ai.scopePin); hold on; 
    plot(ai.time, ai.olfDispPin);
    plot(ai.time, ai.Frame_Out);
    xlabel('time (sec)');
    ylabel('V');
    ylim([-1 6])
    title('thorsync lines', 'Interpreter', 'none');
    legend({'scopePin', 'olfDispPin', 'Frame_Out'}, 'Location', 'best', 'Interpreter', 'none');
    legend('boxoff');
end
