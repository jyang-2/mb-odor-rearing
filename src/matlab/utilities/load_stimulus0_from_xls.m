function tbl = load_stimulus0_from_xls(fname)
%LOAD_STIMULUS_FROM_XLS Loads the contents of sheet 'stimulus' from a .xlsx
%                        file
% inputs: 
%   - fname : filename, should be .xlsx file
% outputs:
%   -   tbl : returns table containing the contents of 'stimulus' sheet
%   Detailed explanation goes here
tbl = readtable(fname,'Sheet', 'stimulus0','TextType', 'string');
end

