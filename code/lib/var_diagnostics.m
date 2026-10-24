function G = var_diagnostics(M, hmax)
% Stability and residual diagnostics for the OLS VAR:
% - largest companion root;
% - multivariate Breusch-Godfrey LM test for residual autocorrelation up to
%   order h = 1..hmax (Lutkepohl, 2005, Sec. 4.4.4), chi2(h n^2);
% - Jarque-Bera normality and ARCH(4) LM tests per equation.

n = M.n;
T = M.T;
U = M.U;
G.maxroot = var_maxroot(M.B, n, M.p, M.ndet);

G.lm = zeros(hmax, 3);                         % [stat df pval]
SigR = (U'*U)/T;
for h = 1:hmax
    Ul = zeros(T, n*h);
    for l = 1:h
        Ul(l+1:end, (l-1)*n + (1:n)) = U(1:end-l, :);
    end
    Xa = [M.X, Ul];
    E = U - Xa*(Xa\U);
    SigE = (E'*E)/T;
    stat = T*(n - trace(SigR\SigE));
    df = h*n^2;
    G.lm(h, :) = [stat df chi2_pval(stat, df)];
end

G.jb = zeros(n, 2);                            % [stat pval]
G.arch = zeros(n, 2);
for i = 1:n
    u = U(:, i);
    z = (u - mean(u))/std(u, 1);
    sk = mean(z.^3);
    ku = mean(z.^4);
    jb = T/6*(sk^2 + (ku - 3)^2/4);
    G.jb(i, :) = [jb chi2_pval(jb, 2)];

    q = 4;
    e2 = u.^2;
    Xa = ones(T - q, 1);
    for l = 1:q
        Xa = [Xa, e2(q+1-l:end-l)]; %#ok<AGROW>
    end
    yy = e2(q+1:end);
    res = yy - Xa*(Xa\yy);
    r2 = 1 - sum(res.^2)/sum((yy - mean(yy)).^2);
    G.arch(i, :) = [(T - q)*r2 chi2_pval((T - q)*r2, q)];
end
end

function pv = chi2_pval(x, k)
pv = gammainc(x/2, k/2, 'upper');
end
