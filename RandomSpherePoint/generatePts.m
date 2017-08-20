% generate points on a unit sphere from a uniform
% distribution: n is number of pts, d is dimension of n-sphere

function pts = generatePts(n, d)
    pts = [];
    for i=1:n
        pts = horzcat(pts, randomSpherePt(d));
    end
end