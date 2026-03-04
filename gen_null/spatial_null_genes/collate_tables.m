load('./results/gene_names.mat'); % loads variable gene_names (cell array)
if iscell(gene_names)
    gene_names = string(gene_names); % convert to string array
end

n_files = 1000;
null_genes = cell(1, n_files);

for i = 1:n_files
    data = load(sprintf('./results/gene_null/null_genes_%d.mat', i));
    null_slice = data.null_genes; % 66 x 10297
    null_genes{i} = array2table(null_slice, 'VariableNames', cellstr(gene_names));
end

null_genes = single(null_genes);
%save("./results/null_genes.mat", "null_genes", "-v7.3", "-nocompression")
%save("/fs04/scratch2/kg98/oldscratch/gchan/SIR_SCZ/SIR_simulator/data/workspace_S132.mat", "null_genes", "-append", "-nocompression");
save("/fs04/scratch2/kg98/oldscratch/gchan/SIR_SCZ/SIR_simulator/data/workspace_S132_null.mat", "null_genes", "-v7.3", "-nocompression");

% Now tables is a 1x1000 cell array of tables, each with gene names as columns