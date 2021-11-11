function newhi = renumber(hi, newmapping)
%RENUMBER  renumber(A, newmapping)
%
%   reassigns values in matrix according to mapping

uhi = unique(hi);

if nargin==1
    newmapping = 1:numel(uhi);
elseif nargin==2
    if ~iscolumn(newmapping) && ~isrow(newmapping)
        assert(numel(uhi)==size(newmapping,1),...
            'Size of newmapping does not match # of values required');

    end
    
end

newhi = hi;


for i = 1:numel(uhi)
    idx = hi==uhi(i);
    if size(newmapping,2)==2
        l = newmapping(newmapping(:,1)==uhi(i),2);
        newhi(idx) = l;
    elseif isvector(newmapping)
        newhi(idx) = newmapping(i);
    end
end



end

