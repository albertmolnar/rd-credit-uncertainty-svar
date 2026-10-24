function R = chol_identify(M, cfg)
% Recursive (Cholesky) identification in the variable order of cfg.vars, used as a
% robustness check. Returns the credit and uncertainty shocks in the same format
% as sign_identify.

n = M.n;
ic = find(strcmp(M.vars, cfg.roleCredit));
iu = find(strcmp(M.vars, cfg.roleUnc));
rest = setdiff(1:n, [ic iu]);
R.irf = zeros(cfg.nDraws, cfg.H+1, n, 2);
R.B = zeros(M.m, n, cfg.nDraws);
R.A = zeros(n, n, cfg.nDraws);
kept = 0; tried = 0;
while kept < cfg.nDraws && tried < cfg.maxPostDraws
    tried = tried + 1;
    [B, Sigma] = var_posterior_draw(M);
    if var_maxroot(B, n, M.p, M.ndet) >= 1, continue; end
    A = chol(Sigma, 'lower');
    kept = kept + 1;
    R.irf(kept, :, :, :) = reshape(var_irf(B, n, M.p, M.ndet, A(:, [ic iu]), cfg.H), [1, cfg.H+1, n, 2]);
    R.B(:, :, kept) = B;
    R.A(:, :, kept) = A(:, [ic iu rest]);
end
R.irf = R.irf(1:kept, :, :, :);
R.B = R.B(:, :, 1:kept);
R.A = R.A(:, :, 1:kept);
R.shockNames = {'Credit supply', 'Uncertainty'};
R.nKept = kept;
R.nTried = tried;
R.accRate = kept/tried;
R.explosiveShare = 1 - kept/tried;
end
