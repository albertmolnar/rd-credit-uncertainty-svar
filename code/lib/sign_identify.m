function R = sign_identify(M, S, cfg)
% Posterior draws of a set-identified SVAR (Rubio-Ramirez, Waggoner and Zha, 2010).
% For each reduced-form draw, Haar-distributed rotations Q are tried until the
% impact matrix A = chol(Sigma)*Q satisfies every restriction in S over
% h = 0..cfg.Hrestr. Columns may match shocks in any order and with either sign.
% Optional narrative restrictions fix the sign of a structural shock at given
% dates (Antolin-Diaz and Rubio-Ramirez, 2018). Their importance weight is the
% inverse of the probability, under N(0,I) shocks, that the restrictions hold.
% For restrictions on the SIGN of shocks only, that probability is 2^-k for
% every draw (our derivation, not stated in the paper), so the weights are
% constant and reweighting has no effect. Restrictions on the historical
% decomposition would need reweighting.

n = M.n; p = M.p; ndet = M.ndet; H = cfg.H; Hr = cfg.Hrestr;
k = numel(S);
nKeep = cfg.nDraws;

narr = zeros(0, 3);                      % [shock, row in sample, sign]
if cfg.useNarrative
    for i = 1:size(cfg.narrative, 1)
        s = find(strcmp({S.name}, cfg.narrative{i, 1}));
        qq = cfg.narrative{i, 2}(1)*4 + cfg.narrative{i, 2}(2) - 1;
        row = find(M.qidx == qq);
        if ~isempty(s) && ~isempty(row)
            narr(end+1, :) = [s row cfg.narrative{i, 3}]; %#ok<AGROW>
        end
    end
end

R.irf = zeros(nKeep, H+1, n, k);
R.B = zeros(M.m, n, nKeep);
R.A = zeros(n, n, nKeep);
kept = 0; tried = 0; explosive = 0;
while kept < nKeep && tried < cfg.maxPostDraws
    tried = tried + 1;
    [B, Sigma] = var_posterior_draw(M);
    if var_maxroot(B, n, p, ndet) >= 1
        explosive = explosive + 1;
        continue
    end
    P = chol(Sigma, 'lower');
    [~, Psi] = var_irf(B, n, p, ndet, [], Hr);
    sd = sqrt(diag(Sigma));
    for it = 1:cfg.nRotMax
        [Q, RR] = qr(randn(n));
        Q = Q*diag(sign(diag(RR)));
        A = P*Q;
        [ok, cols, sgn] = check_restrictions(A, Psi, S, sd, Hr);
        if ~ok, continue; end
        Aid = A(:, cols).*sgn;
        if ~isempty(narr)
            e = A\(M.Y - M.X*B)';                   % all structural shocks, n x T
            es = e(cols, :).*sgn(:);                % identified shocks, k x T
            if any(sign(es(sub2ind(size(es), narr(:, 1), narr(:, 2)))) ~= narr(:, 3))
                continue
            end
        end
        kept = kept + 1;
        R.irf(kept, :, :, :) = reshape(var_irf(B, n, p, ndet, Aid, H), [1, H+1, n, k]);
        R.B(:, :, kept) = B;
        R.A(:, :, kept) = [Aid, A(:, setdiff(1:n, cols))];
        break
    end
end
if kept < nKeep
    warning('Only %d of %d draws accepted (%d reduced-form draws).', kept, nKeep, tried);
    R.irf = R.irf(1:kept, :, :, :);
    R.B = R.B(:, :, 1:kept);
    R.A = R.A(:, :, 1:kept);
end
R.shockNames = {S.name};
R.nKept = kept;
R.nTried = tried;
R.accRate = kept/max(tried - explosive, 1);
R.explosiveShare = explosive/tried;
fprintf('  accepted %d draws; acceptance %.3f; explosive draws discarded %.1f%%\n', ...
    kept, R.accRate, 100*R.explosiveShare);
end

function [ok, cols, sgn] = check_restrictions(A, Psi, S, sd, Hr)
n = size(A, 1);
k = numel(S);
Rh = zeros(n, n, Hr+1);
for h = 0:Hr
    Rh(:, :, h+1) = Psi(:, :, h+1)*A;
end
free = true(1, n);
cols = zeros(1, k);
sgn = zeros(1, k);
ok = false;
for s = 1:k
    found = false;
    for j = find(free)
        resp = reshape(Rh(:, j, :), n, Hr+1);
        for sg = [1 -1]
            if shock_matches(sg*resp, S(s), sd)
                cols(s) = j; sgn(s) = sg; free(j) = false; found = true;
                break
            end
        end
        if found, break; end
    end
    if ~found, return; end
end
ok = true;
end

function tf = shock_matches(resp, s, sd)
tf = false;
for i = 1:size(s.sign, 1)
    if any(s.sign(i, 2)*resp(s.sign(i, 1), :) <= 0), return; end
end
for i = 1:size(s.rel, 1)
    a = s.rel(i, 1); b = s.rel(i, 2);
    if resp(a, 1)/sd(a) <= resp(b, 1)/sd(b), return; end
end
tf = true;
end
