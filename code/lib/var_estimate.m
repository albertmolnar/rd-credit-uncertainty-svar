function M = var_estimate(Y, p, Z)
% OLS estimation of a VAR(p) in levels, Y_t = C z_t + A_1 Y_{t-1} + ... + A_p Y_{t-p} + u_t.
% Rows of B: deterministic terms first, then lag 1 (n rows), lag 2, ...;
% A_l = B(ndet + (l-1)*n + (1:n), :)'.

[T, n] = size(Y);
if nargin < 3 || isempty(Z), Z = ones(T, 1); end
Ye = Y(p+1:end, :);
X = Z(p+1:end, :);
for l = 1:p
    X = [X, Y(p+1-l:end-l, :)]; %#ok<AGROW>
end
B = X\Ye;
U = Ye - X*B;
[Te, m] = size(X);

M.B = B;
M.Bhat = B;
M.U = U;
M.X = X;
M.Y = Ye;
M.XXi = inv(X'*X);
M.n = n;
M.p = p;
M.ndet = size(Z, 2);
M.m = m;
M.T = Te;
M.Sigma = (U'*U)/(Te - m);
M.SigmaML = (U'*U)/Te;
end
