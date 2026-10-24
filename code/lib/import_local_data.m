function imported = import_local_data(cfg)
% Converts the manually downloaded FRED files in cfg.localDir (CSV or XLSX as
% saved from the FRED website) into data_raw/<id>.csv, the format used by
% build_dataset. Daily and monthly series are kept at their native frequency;
% read_series_quarterly averages them to quarters.
% Optional files are picked up when present, e.g. an R&D price index or real
% R&D series saved under its FRED ID.

map = {   % local id      candidate file stems in cfg.localDir
    'GDP',      {'GDP'}
    'GDPC1',    {'GDPC1'}
    'GDPDEF',   {'GDPDEF'}
    'RD_NOM',   {'Y006RC1Q027SBEA'}
    'RD_REAL',  {'Y006RX1Q020SBEA'}
    'RD_PRICE', {'Y006RD3Q086SBEA', 'Y006RG3Q086SBEA'}
    'BAA10Y',   {'BAA10Y'}
    'GS1',      {'DGS1', 'GS1'}
    'EPU',      {'USEPUINDXD', 'USEPUINDXM'}
    'VIXCLS',   {'VIXCLS'}
    'EBP',      {'ebp_csv', 'EBP'}
    };
imported = strings(0, 2);
for i = 1:size(map, 1)
    id = map{i, 1};
    for stem = map{i, 2}
        f = '';
        for ext = {'.csv', '.xlsx', '.xls'}
            cand = fullfile(cfg.localDir, [stem{1} ext{1}]);
            if isfile(cand), f = cand; break; end
        end
        if isempty(f), continue; end
        [d, v] = read_any(f, id);
        out = fullfile(cfg.rawDir, [id '.csv']);
        T = table(string(d, 'yyyy-MM-dd'), v, 'VariableNames', {'observation_date', id});
        writetable(T, out);
        imported(end+1, :) = [string(id), string(f)]; %#ok<AGROW>
        break
    end
end
disp(array2table(imported, 'VariableNames', {'id', 'source_file'}));
end

function [d, v] = read_any(f, id)
[~, ~, ext] = fileparts(f);
if strcmpi(ext, '.csv')
    T = readtable(f, 'TextType', 'string');
else
    sh = sheetnames(f);
    sh = sh(~strcmpi(sh, 'README'));
    T = readtable(f, 'Sheet', sh(1));
end
d = T{:, 1};
if ~isdatetime(d), d = datetime(string(d), 'InputFormat', 'yyyy-MM-dd'); end
if strcmp(id, 'EBP')
    v = T.ebp;
else
    v = T{:, 2};
end
if ~isnumeric(v), v = str2double(string(v)); end
ok = ~isnat(d) & ~isnan(v);
d = d(ok);
v = v(ok);
end
