function ccids = multi_movie_cells(all_matchings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% find the cell IDs that are constant across movies

ccids=struct();
n_expt = max(size(all_matchings));

for i = 1:n_expt
    K = size(all_matchings(i,i),1);
    common_cids = zeros(K, n_expt);
    for j = 1:n_expt
        matching = all_matchings{i,j};
        common_cids(matching(:,1),j) = 1;
    end
    ccids(i).mat = common_cids;
    ccids(i).cids = find(all(common_cids,2));
end

for i=1:n_expt
   ccids(i).mcids = [];
   for j = 1:n_expt
       matching = all_matchings{i,j};
       
       [umatching, ia] = unique(matching(:,1), 'first');
       matching = matching(ia,:);
       
       Lia = ismember(matching(:,1), ccids(i).cids);
       disp([i j sum(Lia)])
       mcids = matching(:,2);
       mcids = mcids(Lia,:);
       
%        if numel(mcids) ~= numel(ccids(i).cids)
%           keyboard; 
%        end
       ccids(i).mcids = [ccids(i).mcids mcids];
   end
end
%%
end

