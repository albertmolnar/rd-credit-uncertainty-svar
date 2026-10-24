function [q, v] = read_series_quarterly(file, valCol)
% Reads a CSV whose first column is a date and aggregates the value column
% (index or case-insensitive name) to quarterly averages.
% q is the quarter index year*4 + quarter - 1. Incomplete quarters at the edges
% (fewer than 60% of the usual number of observations) are dropped.

opts = detectImportOptions(file, 'TextType', 'string');
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(file, opts);
if isnumeric(valCol)
    vcol = T{:, valCol};
else
    k = find(strcmpi(T.Properties.VariableNames, valCol), 1);
    if isempty(k), error('Column "%s" not found in %s', valCol, file); end
    vcol = T{:, k};
end
d = parse_dates(T{:, 1});
x = str2double(vcol);              % FRED codes missing values as "."
ok = ~isnat(d) & ~isnan(x);
d = d(ok);
x = x(ok);
[d, ord] = sort(d);
x = x(ord);

qi = year(d)*4 + quarter(d) - 1;
[q, ~, g] = unique(qi);
v = accumarray(g, x, [], @mean);
cnt = accumarray(g, 1);
keep = cnt >= 0.6*median(cnt);
q = q(keep);
v = v(keep);
end

function d = parse_dates(s)
s = strtrim(string(s));
fmts = {'yyyy-MM-dd', 'MM/dd/yyyy', 'M/d/yyyy', 'yyyy-MM', 'yyyyMMdd'};
d = NaT(size(s));
for k = 1:numel(fmts)
    todo = isnat(d);
    if ~any(todo), break; end
    try
        d(todo) = datetime(s(todo), 'InputFormat', fmts{k});
    catch
    end
end
end
