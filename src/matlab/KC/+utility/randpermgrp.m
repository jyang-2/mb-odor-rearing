function perm_idx = randpermgrp(grp, niter)
%RANDPERMGRP Summary of this function goes here
%   Detailed explanation goes here
[ugrp, grpcounts] = uniquecount(grp);
ngrp = numel(ugrp); % # of groups

if nargin == 1
   niter =  1; 
end

perm_idx = zeros(size(grp,1), niter);

for ni = 1:niter
    for i = 1:ngrp
        igrp = ugrp(i);
        grpidx = find(grp == igrp);
        perm_idx(grpidx,ni) = grpidx(randperm(grpcounts(i)));
    end
end

end

