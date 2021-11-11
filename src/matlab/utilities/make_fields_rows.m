function S2 = make_fields_rows(S)
%MAKE_FIELDS_ROWS Summary of this function goes here
%   Detailed explanation goes here
    fns = fields(S);
    S2 = S;
    for i = 1:length(fns)
        f = fns{i};
        S2.(f) = torow(S2.(f));
    end

end

