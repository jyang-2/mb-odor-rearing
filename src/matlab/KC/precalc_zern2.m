function [mmnts, nm] = precalc_zern2(cm, img3d, order)
%PRECALC_ZERN2 Summary of this function goes here
%   Detailed explanation goes here

%simfun = @vcorr;

% Precalculate moment vectors for each point
% order = 25; % Maximum order for moment calculation
nm = string(); % List of subscripts
for n=0:order
    for m=0:n
        if mod(n-abs(m), 2)==0
            nm = [nm; [num2str(n),num2str(m)]];
        end
    end
end
nm = nm(2:end);


K = size(cm,1);  % number of centroids
icm = round(cm);    % rounded centroids
z = icm(:,3);       % nearest plane of each centroid

% Pre-calculate Zernike moments for all scene and model points
mmnts = zeros(K,numel(nm));
for i=1:K
    im = img3d(:,:,z(i));
    mmnts(i,:) = zern2(submask(im, cm(i, 1:2)),order).';
end


end

