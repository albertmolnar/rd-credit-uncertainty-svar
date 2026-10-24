function [B, Sigma] = var_posterior_draw(M)
% One draw from the normal-inverse-Wishart posterior of the VAR under the
% diffuse prior (Uhlig, 2005): Sigma ~ IW(U'U, T-m), vec(B) | Sigma ~ N(vec(Bhat), Sigma (x) (X'X)^-1).

n = M.n;
nu = M.T - M.m;
S = M.U'*M.U;
Ci = chol(inv(S), 'lower');
Zw = randn(nu, n)*Ci';
Sigma = inv(Zw'*Zw);
Sigma = (Sigma + Sigma')/2;
L = chol(M.XXi, 'lower');
B = M.Bhat + L*randn(M.m, n)*chol(Sigma);
end
