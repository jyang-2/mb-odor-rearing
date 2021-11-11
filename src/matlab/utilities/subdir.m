function [subFolders] = subdir(varargin)
%SUBDIR - % Get a list of all folders in this folder
%   subsubdir() - returns list of subfolders in current working directory
%   subdir(folder) - returns list of subfolders in folder
%   subdir(folder, 'verbose', false) - returns list of subfolder names w/o
%   printing output
%   
% 
% inputs:
%   folder
% optional inputs:
%   subdir(..., 'verbose', false) is default

default_folder = pwd;
default_verbose = true;

p = inputParser;
validFolderName = @(x) ischar(x) || isStringScalar(x);
%addRequired(p,'folder',@isfolder);
addOptional(p,'folder', default_folder, validFolderName);
addParameter(p,'verbose',default_verbose, @isboolean);
parse(p,varargin{:});
folder = p.Results.folder;
verbose = p.Results.verbose;
   

files = dir(folder);

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

% Print folder names to command window.
if verbose
    fprintf('Directory: %s\n', folder);
    for k = 1 : length(subFolders)
        fprintf('\tSub folder #%d = %s\n', k, subFolders(k).name);
    end
end

end

