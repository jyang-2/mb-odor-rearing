function pathout = remyDropboxPath()
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

s = regexp(dropboxPath, filesep, 'split');
s(cellfun(@isempty, s)) = [];

pathout = fullfile('/', s{1:end-1}, 'Remy');


end

