PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb-odor-rearing";
tiff_folder = fullfile(PRJ_DIR, "1", "1");
tiff_files = dir(fullfile(tiff_folder, "**/stk*.tif"));
%%
for i = 1:length(tiff_files)
    Y = loadtiff(tiff_files(i).listing, tiff_files(i).name);
    
    Yr = reshape(Y, size(Y,1), size(Y, 2), 16, []);
    
    Ysub = zeros(size(Y,1), size(Y,2), 4, 
    
    for t = 1:size(Yr,4)
        yr = Yr(:,:,:,t)
        yr= reshape(size(yr,1), size(yr,2), 3, []);
        
    end
%     
end