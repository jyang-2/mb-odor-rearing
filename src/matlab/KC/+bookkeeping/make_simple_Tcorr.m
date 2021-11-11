function Tcorr = make_simple_Tcorr(corrlist, response, ti, bad_trials)
%MAKE_SIMPLE_TCORR     builds table of correlations from struct response
%
%  INPUTS:
%       corrlist: cellstr, list of correlation/distance matrix types
%       response: response struct containing corr. matrices of different types
%             ti: struct holding stimulus, frame information for movie
%     bad_trials: trials that should be ignored, [1 x num_stim], logical
%
%  OUTPUTS
%          Tcorr: table, w/ row for each good odor pair 
%                 variables (columns) = hi1, hi2, si1, si2, various correlations
%                 does not include date, fly, or experiment name
%% check if bad_trials exists

if ~exist('bad_trials', 'var') || isempty(bad_trials)
   bad_trials = false(size(ti.si)); 
end

num_pairs = ti.num_stim*(ti.num_stim-1)/2;
%% build Rs (struct) and Tcorr (table) with stimulus info and corr. values
Rs = struct();
for i = 1:numel(corrlist)
   Rs.(corrlist{i}) = nan(num_pairs,1);
end

hi1 = [];
hi2 = [];
si1 = [];
si2 = [];
bad_pair = [];


for i = 1:ti.num_stim
   for j = (i+1):ti.num_stim
       bad_pair = [bad_pair; any([bad_trials(i), bad_trials(j)])];

       [hi_s, idx] = sort([ti.hi(i) ti.hi(j)], 'ascend');
       hi1 = [hi1; hi_s(1)];
       hi2 = [hi2; hi_s(2)];
       
       stim = [ti.si(i) ti.si(j)];
       stim = stim(idx);
       si1 = [si1; stim(1)];
       si2 = [si2; stim(2)];
       
       itemp = length(hi1);
       for ic = 1:numel(corrlist)
          Rs.(corrlist{ic})(itemp) = response.(corrlist{ic})(i,j);
       end
   end
end

Tcorr = table(hi1, hi2, si1, si2);

for i = 1:numel(corrlist)
   Tcorr.(corrlist{i}) = Rs.(corrlist{i});
end

bad_pair = logical(bad_pair);
Tcorr(bad_pair,:) = [];
end