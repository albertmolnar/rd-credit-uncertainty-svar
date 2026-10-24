function Q = pctl(X, q)
% Percentiles along the first dimension with linear interpolation
% (base MATLAB has no prctile without the Statistics Toolbox).
sz = size(X);
N = sz(1);
Xs = sort(reshape(X, N, []), 1);
Q = zeros(numel(q), size(Xs, 2));
for i = 1:numel(q)
    pos = 1 + (N - 1)*q(i)/100;
    lo = floor(pos);
    hi = ceil(pos);
    w = pos - lo;
    Q(i, :) = (1 - w)*Xs(lo, :) + w*Xs(hi, :);
end
Q = reshape(Q, [numel(q) sz(2:end)]);
end
