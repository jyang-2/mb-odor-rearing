function dbpath = hongLabDropboxRoot()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ismac
    % Code to run on Mac platform
elseif isunix
    [~,name] = system('hostname');
    name = strtrim(name);

elseif ispc
    name = getenv('COMPUTERNAME');
else
    disp('Platform not supported')
end

switch name
    case 'DARWIN' % home desktop computer
        dbpath = 'G:\HongLab @ Caltech Dropbox\Remy';
    case 'gerty' % honglab analysis workstation
         dbpath = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy";
    case 'hex'   % remy's laptop (ubuntu)
         dbpath = "/home/remy/HongLab @ Caltech Dropbox/Remy";
    otherwise 
        
end

end

