function S = make_restrictions(cfg, vars)
% Sign and relative-magnitude restrictions. Rows of S(k).sign are [variable, sign];
% rows of S(k).rel are [a b], meaning the impact response of a exceeds that of b
% when both are measured in units of their own residual standard deviation.
% The R&D variable is deliberately left unrestricted for every shock.

ix = @(nm) find(strcmp(vars, nm));
c = ix(cfg.roleCredit);
u = ix(cfg.roleUnc);
y = ix('y');
p = ix('p');
r = ix('r');

S = struct('name', {}, 'sign', {}, 'rel', {});

S(1).name = 'Credit supply';
S(1).sign = [c 1; y -1; p -1];
if cfg.creditRateRestr, S(1).sign = [S(1).sign; r -1]; end
S(1).rel = [c u];

S(2).name = 'Uncertainty';
if cfg.uncUseMacro
    S(2).sign = [u 1; y -1; p -1];
else
    S(2).sign = [u 1];
end
S(2).rel = [u c];

S(3).name = 'Aggregate supply';
S(3).sign = [y -1; p 1];
S(3).rel = zeros(0, 2);

ird = ix(cfg.roleRD);
for k = 1:numel(S)
    assert(~any(S(k).sign(:, 1) == ird), 'R&D must not be restricted.');
end
end
