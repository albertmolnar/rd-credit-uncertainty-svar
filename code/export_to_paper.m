%% Copy tables and figures into the paper folder and build an Overleaf zip
% After this, v2/paper is self-contained:
%   paper_v2.tex, references.bib, IEEEtran.cls, tables/*.tex, figures/*.pdf
% and v2/paper_v2_overleaf.zip can be uploaded to Overleaf as a new project
% (New Project > Upload Project).
here = fileparts(mfilename('fullpath'));
out = fullfile(here, 'output');
paper = fullfile(here, '..', 'paper');
for d = {'tables', 'figures'}
    if ~exist(fullfile(paper, d{1}), 'dir'), mkdir(fullfile(paper, d{1})); end
end

t = dir(fullfile(out, 'tables', '*.tex'));
for i = 1:numel(t)
    copyfile(fullfile(t(i).folder, t(i).name), fullfile(paper, 'tables', t(i).name));
end
f = dir(fullfile(out, 'figures', '*.pdf'));
for i = 1:numel(f)
    name = strrep(f(i).name, '.jpg.pdf', '.pdf');     % tolerate a renamed file
    copyfile(fullfile(f(i).folder, f(i).name), fullfile(paper, 'figures', name));
end
fprintf('Copied %d tables and %d figures into %s\n', numel(t), numel(f), paper);

zipfile = fullfile(paper, '..', 'paper_v2_overleaf.zip');
zip(zipfile, {'paper_v2.tex', 'references.bib', 'IEEEtran.cls', 'tables', 'figures'}, paper);
fprintf('Overleaf upload file: %s\n', zipfile);
