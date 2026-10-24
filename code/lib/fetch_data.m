function src = fetch_data(cfg)
% Downloads every configured source series once into cfg.rawDir. Files already
% present are kept, so the raw folder of a run is a fixed data vintage; delete a
% file to refresh it. Returns a table documenting source, URL and vintage date.

opts = weboptions('Timeout', 60);
rows = cell(0, 5);
for i = 1:size(cfg.series, 1)
    id = cfg.series{i, 1};
    url = cfg.series{i, 2};
    desc = cfg.series{i, 3};
    f = fullfile(cfg.rawDir, [id '.csv']);
    if isempty(url)
        status = "not configured";
    elseif isfile(f)
        status = "cached";
    else
        try
            websave(f, url, opts);
            status = "downloaded";
        catch ME
            status = "FAILED: " + string(ME.message);
            warning('Could not download %s from %s (%s).\nDownload it manually and save it as %s', ...
                id, url, ME.message, f);
        end
    end
    vintage = "";
    if isfile(f)
        info = dir(f);
        vintage = string(datestr(info.datenum, 'yyyy-mm-dd'));
    end
    rows(end+1, :) = {string(id), string(desc), string(url), status, vintage}; %#ok<AGROW>
end
src = cell2table(rows, 'VariableNames', {'id', 'description', 'url', 'status', 'vintage'});
disp(src(:, {'id', 'status', 'vintage'}));
end
