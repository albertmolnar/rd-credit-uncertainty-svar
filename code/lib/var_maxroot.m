function r = var_maxroot(B, n, p, ndet)
% Largest modulus of the companion-matrix eigenvalues (< 1: stable VAR).
F = zeros(n*p);
F(1:n, :) = B(ndet+1:ndet+n*p, :)';
if p > 1
    F(n+1:end, 1:n*(p-1)) = eye(n*(p-1));
end
r = max(abs(eig(F)));
end
