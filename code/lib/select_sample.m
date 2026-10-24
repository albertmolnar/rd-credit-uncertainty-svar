function [Y, Z, qidx] = select_sample(D, cfg)
% Selects the estimation sample (contiguous, no missing values) and builds the
% deterministic terms: a constant, plus COVID impulse dummies if requested.

q0 = cfg.sampleStart(1)*4 + cfg.sampleStart(2) - 1;
if isempty(cfg.sampleEnd)
    q1 = inf;
else
    q1 = cfg.sampleEnd(1)*4 + cfg.sampleEnd(2) - 1;
end
Ym = D{:, cfg.vars};
in = D.qidx >= q0 & D.qidx <= q1 & all(~isnan(Ym), 2);
first = find(in, 1);
last = find(in, 1, 'last');
if isempty(first), error('No observations for variables %s.', strjoin(cfg.vars, ', ')); end
if any(~in(first:last))
    error('Missing values inside the sample %s-%s.', D.quarter(first), D.quarter(last));
end
Y = Ym(first:last, :);
qidx = D.qidx(first:last);
Z = ones(size(Y, 1), 1);
if cfg.covidDummies
    for qq = (2020*4):(2021*4 + 3)
        if any(qidx == qq), Z = [Z, double(qidx == qq)]; end %#ok<AGROW>
    end
end
fprintf('Sample %s-%s (%d quarters incl. pre-sample), variables: %s\n', ...
    D.quarter(first), D.quarter(last), numel(qidx), strjoin(cfg.vars, ', '));
end
