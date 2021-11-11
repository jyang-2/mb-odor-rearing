function A_cat = row_horzcat(mcids, A)
%UNTITLED Summary of this function goes here
%   A : cell array of matrices to concatenate

[K,n] = size(mcids);
assert(n==numel(A));

A_cat = [];

for i = 1:n
    ridx = mcids(:,i);
    dat = A{i};
    A_cat = [A_cat, dat(ridx,:)];
end

