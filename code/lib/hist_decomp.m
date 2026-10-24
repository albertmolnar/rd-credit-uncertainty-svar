function [contrib, baseline, e] = hist_decomp(M, B, A)
% Historical decomposition. contrib(t, i, s) is the part of variable i at date t
% due to structural shock s since the start of the sample; baseline is the rest
% (initial conditions and deterministic terms). A must be a full impact matrix
% with A*A' = Sigma; e holds the structural shocks (T x n).

n = M.n;
T = M.T;
U = M.Y - M.X*B;
e = (A\U')';
[~, Psi] = var_irf(B, n, M.p, M.ndet, [], T-1);
contrib = zeros(T, n, n);
for t = 1:T
    acc = zeros(n);
    for j = 0:t-1
        acc = acc + (Psi(:, :, j+1)*A).*e(t-j, :);
    end
    contrib(t, :, :) = reshape(acc, [1 n n]);
end
baseline = M.Y - sum(contrib, 3);
end
