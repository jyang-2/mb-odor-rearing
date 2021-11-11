PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data"; 
cd(PRJ_DIR)
listing = dir('**/TENSOR_RAW.mat*');


%%

for i = 27
    TENSOR_RAW = load(fullfile(listing(i).folder, listing(i).name));
    Fall = load(fullfile(listing(i).folder, 'Fall.mat'));
    
    iscell = Fall.iscell(:,1);
    cellprob = Fall.iscell(:,2);
    
    % fix df_stimulus, make odors categorical
    ntmeta = load(fullfile(listing(i).folder, "nt_meta_data.mat"));
    ntmeta.meta_yaml = jsondecode(ntmeta.meta_yaml);

    df_stimulus = struct2table(ntmeta.df_stimulus);
    for field = ["odor_a", "odor_b", "stimname"]
        df_stimulus.(field) = string(df_stimulus.(field));
    end
    valueset = ["pfo", "1-5ol", "VA", "1-6ol", "EP", "1-6ol+EP"];
    df_stimulus.stimnames = categorical(df_stimulus.stimname,valueset,'Ordinal',true);
    [sorted_trials, sort_idx] = sort(df_stimulus.stimname);

    F = Fall.F;
    Fneu = Fall.F;
    
    [n_cells, T] = size(F);
    F_block = mat2cell(F, n_cells, [1, 1, 1] * T/3);
    F_block = cat(3, F_block{:});
    %%
    
    
    tensor_F = TENSOR_RAW.tensor_F;
    tensor_F = permute(tensor_F, [1 3 2]);  
    
    tensor_Fneu = TENSOR_RAW.tensor_Fneu;
    tensor_Fneu = permute(tensor_Fneu, [1 3 2]);  
%%
[~,index_A,index_B] = intersect(Y,brushedData,'rows');
end
R = 20;
model_blk = cp_als(tensor(F_block), R, maxiter=100);

info = viz_ktensor(model_blk, ...
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'});
        
%%
tensor_Fc = tensor_F - 0.7*tensor_Fneu;


R = 15;
model_Fc = cp_als(tensor(tensor_Fc(cid, :, :)), R, 'tol', 1e-6); % fit CP model with 10 components

info = viz_ktensor(model_Fc, ...
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'trials', 'trials'});
%%

        
model_Fneu = cp_als(tensor(tensor_Fneu), 5); % fit CP model with 10 components

viz_ktensor(model_Fneu, ...
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'});
        
viz_ktensor(model_Fneu) % produces a nice plot for you
%%
tfc = tensor_Fc(iscell, :);  

model_F = cp_als(tensor(tensor_F), 10); % fit CP model with 10 components
F_info = viz_ktensor(model_F, ...
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Plotcolors', {'k', 'k', 'k'},...
            'Modetitles', {'neurons', 'time', 'trials'});

err_F = norm(model_F - tensor_F)/norm(tensor_F);        

%%

hBrush = brush(figure(1));
set( hBrush,'ActionPostCallback', @(ohf, s) brushDataCallback)

%%


function brushDataCallback(~, ~)
    [~,cid, ~] = intersect(Y,brushedData,'rows');
    print(cid)
end


