function s = quarter_label(q)
% Quarter index (year*4 + quarter - 1) to labels such as "1985Q1".
s = string(floor(q/4)) + "Q" + string(mod(q, 4) + 1);
end
