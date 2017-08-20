% generate a point on a unit n-sphere from a uniform distribution

function pt = randomSpherePt(dimension)
    rN = normrnd(0,1,[dimension,1]);
    pt = zeros(dimension,1);
    denominator = 0;
    for i=1:length(rN)
        denominator = denominator + rN(i)^2;
    end
    denominator = denominator^(0.5);
    for i=1:length(rN)
        pt(i) = rN(i)/denominator;
    end
end