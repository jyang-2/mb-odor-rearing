function tbl = load_sessioninfo_from_xls(fname)
%LOAD_STIMULUS_FROM_XLS Loads the contents of sheet 'stimulus' from a .xlsx
%                        file
% inputs: 
%   - fname : filename, should be .xlsx file
% outputs:
%   -   tbl : returns table containing the contents of 'sessioninfo' sheet
%   Detailed explanation goes here
tbl = readtable(fname,'Sheet', 'sessioninfo','TextType', 'string');
end

