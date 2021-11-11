function ti = load_ti(folder)
%LOAD_TI Loads ti from 'ti.mat' in folder
%
% Inputs
%   - folder
%
% Outputs
%   - ti

ti = load(fullfile(folder, 'ti.mat'));
end

