function fe = fevd_draws(R, M, H)
% Forecast-error variance shares for every accepted draw:
% fe(draw, h+1, variable, shock), shocks ordered as in R.A (identified first).
nd = size(R.A, 3);
n = M.n;
fe = zeros(nd, H+1, n, n);
for d = 1:nd
    irf = var_irf(R.B(:, :, d), n, M.p, M.ndet, R.A(:, :, d), H);
    c = cumsum(irf.^2, 1);
    fe(d, :, :, :) = reshape(c./sum(c, 3), [1, H+1, n, n]);
end
end
