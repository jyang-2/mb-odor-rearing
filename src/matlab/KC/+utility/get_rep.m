function rep = get_rep(hi)

[uhi, counts] = uniquecount(hi);

rep = zeros(size(hi));
for i = 1:numel(uhi)
    idx = hi==uhi(i);
    rep(idx) = 1:counts(i);
end

end

