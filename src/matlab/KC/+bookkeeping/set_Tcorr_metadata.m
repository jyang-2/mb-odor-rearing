function Tcorr = set_Tcorr_metadata(Tcorr, varargin)
%SET_TCORR_METADATA Adds metadata (like 'date','fly','nam') to simple Tcorr table
%   INPUTS
%         Tcorr: simple table (no metadata), w/ stim. info and correlation vals
%      varargin: parameter names followed by values
%   
%   OUTPUTS
%         Tcorr: Table w/ metadata added.


p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'Tcorr', @istable);
addParameter(p, 'date', @(x) isdatetime(x) || isstring(x) || ischar(x));
addParameter(p, 'fly', @isnumeric);
addParameter(p, 'nam', @(x) isstring(x) || ischar(x));
parse(p, Tcorr, varargin{:});

dat = repelem( string(p.Results.date) , height(Tcorr), 1);
nam = repelem( string(p.Results.nam)  , height(Tcorr), 1);
fly = repelem( p.Results.fly, height(Tcorr), 1);

Tcorr = addvars(Tcorr,...
    dat, fly, nam, ...
    'Before', 'hi1',...
    'NewVariableNames', {'date', 'fly', 'nam'});


end

