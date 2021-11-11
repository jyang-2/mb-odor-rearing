function FF2 = get_component_adjacency(A, C, A_thr)
%GET_COMPONENT_ADJACENCY Summary of this function goes here
%   Detailed explanation goes here

A = CNMc.A;
C = CNMc.C;
A_thr = 0;

%% % normalize spatial components to unit energy
nA = full(sqrt(sum(A.^2))');
[K,~] = size(C);
normA = A/spdiags(nA,0,K,K);    
normC = bsxfun(@times,C,nA(:));

%% find graph of overlapping spatial components
nr = size(normA,2);
A_corr = (normA'*normA);                
A_corr(1:nr+1:nr^2) = 0;
FF2 = (A_corr > A_thr);          

end