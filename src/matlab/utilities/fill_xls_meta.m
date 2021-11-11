function [sessioninfo_tbl, stimulus_tbl, olfactometer_tbl] = fill_xls_meta(fname)
%FILL_XLS_META Reads in .xlsx file w/ stimulus, olfactometer information
%              and fills out the remaining data
%   Detailed explanation goes here
%
%   required functions:
%      - load_sessioninfo_from_xls.......src/utility/load_sessioninfo_from_xls.m
%      - load_stimulus_from_xls..........src/utility/load_stimulus_from_xls.m
%      - load_olfactometer_from_xls......src/utility/load_olfactometer_from_xls.m
%      - load_stimtype_from_xls..........src/utility/load_stimtype_from_xls.m
%      - make_fields_rows .,.............src/utility/make_fields_rows.m
%      - getChannelInfo..................src/arduino_stimulus/getChannelInfo.m

if isfolder(fname)
    fname = fullfile(fname, 'metadata.xlsx');
elseif isfile(fname)
    % do nothing
end
sessioninfo_tbl = load_sessioninfo_from_xls(fname);

stimulus_tbl0 = load_stimulus0_from_xls(fname);
olfactometer_tbl = load_olfactometer_from_xls(fname);
writetable(stimulus_tbl0, fname, 'Sheet', 'stimulus0');

stimulus = make_fields_rows(table2struct(stimulus_tbl0,'ToScalar',true));
olfactometer = make_fields_rows(table2struct(olfactometer_tbl,'ToScalar',true));

chanA = df_(stimulus.pinA, olfactometer);
chanB = getChannelInfo(stimulus.pinB, olfactometer);

stimulus_tbl = table(chanA.pin', chanA.hi', chanA.odor', chanA.conc', ...
                    chanB.pin', chanB.hi', chanB.odor', chanB.conc',...
    'VariableNames', {'pinA', 'hi1', 'odorA', 'concA', 'pinB', 'hi2', 'odorB', 'concB'});

writetable(stimulus_tbl, fname, 'Sheet', 'stimulus');

stimtype_tbl = load_stimtype_from_xls(fname);
if ~isempty(stimtype_tbl)
   stimulus_tbl = join(stimulus_tbl, stimtype_tbl); 
   writetable(stimulus_tbl, fname, 'Sheet', 'stimulus');
end


end




