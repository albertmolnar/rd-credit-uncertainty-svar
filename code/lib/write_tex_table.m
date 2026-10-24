function write_tex_table(file, header, body, align)
% Writes a booktabs-free tabular (IEEE style uses \hline) for \input in the paper.
fid = fopen(file, 'w');
fprintf(fid, '\\begin{tabular}{%s}\n\\hline\n', align);
fprintf(fid, '%s \\\\\n\\hline\n', strjoin(header, ' & '));
for i = 1:size(body, 1)
    fprintf(fid, '%s \\\\\n', strjoin(body(i, :), ' & '));
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);
end
