function pairs = pdistpairs(X)
%PDISTPAIRS returns order of pairwise distances returned by pdist
% like (2,1), (3,1), (3,2), etc. as a [k x 2] matrix
%
% Syntax:  pairs = pdistpairs(X)
%
% Inputs:
%    X - matrix input to pdist, or # of columns in a matrix input to pdist
%
% Outputs:
%    pairs - [k x 2] each row corresponds to a pair of columns
%
% Examples: 
%    pairs = pdistpairs(27);
%    
%    X = rand(2000,27);
%    pairs = pdistpairs(X);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
if isscalar(X)
    m = X;
elseif ismatrix(X)
    m = size(X,1);
end

pairs = [];
for i = 1:m
    for j = i:m
        if i~=j
            pairs = [pairs; j i];
        end
    end
end


end

