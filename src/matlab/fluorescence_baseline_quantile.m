if isunix
    PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
elseif ispc
    PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing";
end
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data');

file_list = fullfile(PROC_DIR, "**", "suite2p", "combined", "Fall.mat");
listing = dir(file_list);
%%
for i = 1:height(listing)
    try
        % load suite2p .mat file
        Fall_filepath = fullfile(listing(i).folder, listing(i).name);
        disp(Fall_filepath);
        Fall = load(Fall_filepath);

        % create a vector of length Fall.stat
        iplane = nan(size(Fall.stat));

        % assign every stat element (1 per cell) to it's plane #
        for icell = 1:length(Fall.stat)
            iplane(icell)=Fall.stat{icell}.iplane;
        end

        iplane = int32(iplane);
        writeNPY(iplane, fullfile(listing(i).folder, 'iplane.npy'));
        clearvars('iplane') 
    catch
        warning("Could not create iplane.npy from Fall.mat");
        disp(fullfile(listing(i).folder, listing(i).name))
    end
end
%%

statfolder = fullfile(PROC_DIR, '2021-08-30/1/movie_002/0/suite2p/combined');

Fall = load(fullfile(statfolder, 'Fall.mat'));
Fc = Fall.F - Fall.Fneu*.7;
cellprob=Fall.iscell(:,2);
iplane = readNPY(fullfile(statfolder, 'iplane.npy'));
%%
% 
% iplane = nan(size(Fall.stat));
% for i = 1:numel(iplane)
%     iplane(i)=Fall.stat{i}.iplane;
% end
% 
% iplane = int32(iplane);
% writeNPY(iplane, fullfile(statfolder, 'iplane.npy'));

%% lowpass filter in one piece
Fc_lowpass = lowpass(Fc',0.5, 3)';



% Calculate the 0.3 quantile for each column of X (dim = 1).
qt = quantile(Fc_lowpass, 0.6, 2); % .6 quantile, along 2nd dimension
Fc_lowpass_qt = Fc_lowpass./qt; % F corrected, lowpass filted, quantile normed

%writeNPY(Fc_lowpass, fullfile(folder, 'Fc_lowpass.npy'));
%writeNPY(Fc_lowpass_qt, fullfile(folder, 'Fc_lowpass_qt.npy'));
%% lowpass filter acquisitions separately

[n_cell, T] = size(Fc);
Fc_split = mat2cell(Fc,n_cell, (T/3)*ones(1,3));

Fc_lowpass_block = cell(1, 3);
qt_block = cell(1, 3);
Fc_lowpass_qt_block = cell(1, 3);

for i = 1:3
    Fc_lowpass_block{i} = lowpass(Fc_split{i}', 0.2, 3)';
    qt_block{i} = quantile(Fc_lowpass_block{i}, 0.6, 2) - quantile(Fc_lowpass_block{i}, 0.05, 2);
    Fc_lowpass_qt_block{i} = Fc_lowpass_block{i}./qt_block{i};
    %qt_block{i} = quantile(Fc_split{i}, 0.7, 2);
    %Fc_qt_block{i} = Fc_split{i}./qt_block{i};
end

Fc_lowpass = cell2mat(Fc_lowpass_block);
qt = cell2mat(qt_block);
% Fc_lowpass_qt = cell2mat(Fc_lowpass_qt_block);
%%
Fc_lowpass_qt_cell = cell(1,3);
for i = 1:3
    fc = Fc_split{i};
    Fc_lowpass_qt_cell{i} = qt_filt_traces(fc, 3.0);
end

%% 


planecids = (iplane'==5);
iscell = (cellprob > 0.5);

fc_lowpass_qt = Fc_lowpass_qt(iscell, :);
fc_lowpass_qtqt = fc_lowpass_qt./quantile(fc_lowpass_qt, 0.6, 2);
cid0 = 1;
cidlist = (1:10)+cid0-1;
figure, s=stackedplot(fc_lowpass_qt(cidlist+cid0,:)', 'DisplayLabels', string(cidlist));

for i = 1:10
    s.AxesProperties(i).YLimits = [-0.5, 5];
end
%%
% plot first 10 rows
iblock = 1;

figure
%s = stackedplot(Fc_lowpass_qt_cell{1}(1:10,:)');
s = stackedplot(Fc_lowpass_qt_block{iblock}(1:10, :)');

for i = 1:10
    s.AxesProperties(i).YLimits = [-0.5, 10];
end
%%
TENSOR_TRIALS = load(fullfile(statfolder, 'TENSOR_TRIALS.mat'));

zdfof_tensor = TENSOR_TRILAS.zdfof_tensor;
zdfof_tensor = permute(zdfof_tensor, [1 3 2]);
data = zdfof_tensor;

model = cp_als(tensor(data), 10); % fit CP model with 10 components
viz_ktensor(model) % produces a nice plot for you
%%
 [~,index_in_b,index_in_Y] = intersect(brushedData,Y,'rows');
for row = 1:size(brushedData,1)
    
end
%%
% fit the cp decomposition from random initial guesses

R = 15;
n_fits = 10;
err = zeros(n_fits,1);
for n = 1:n_fits
    % fit model
    est_factors = cp_als(tensor(data),R, 'tol', 1e-6);
    
    if n==1
        first_factors = est_factors;
    end
    
    % store error
    err(n) = norm(full(est_factors) - data)/norm(data);
    
    
    % visualize fit for first several fits
    if n > 1
        % score aligns the cp decompositions
        [sc, est_factors] = score(est_factors, first_factors);
        
        % plot the estimated factors
        viz_ktensor(est_factors, ... 
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'})
        set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
   end
end
%%
function fc_lowpass_qt = qt_filt_traces(fc, fs)
    fc_lowpass = lowpass(fc, 0.5, fs);
    qt = quantile(fc_lowpass, 0.5, 2, 'Method', 'approximate');
    fc_lowpass_qt = fc_lowpass_qt./qt;
    
end
