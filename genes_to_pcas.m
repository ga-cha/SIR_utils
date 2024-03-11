% run pca
gene_arr = table2array(genes);
[coeff, score, latent, ~, explained] = pca(gene_arr);

% choose n PCs to retain 95% variance
n_pcs = find(cumsum(explained) >= 95, 1);

% take the most heavily weighted columns in each PC
n_cols = 15;
w_cols = cell(n_pcs, 1);
for i = 1:n_pcs
    [~, sorted] = sort(abs(coeff(:,i)),'descend');
    w_cols{i} = sorted(1:n_cols);
end

% store PCs
gene_pcs_arr = score(:, 1:n_pcs);
gene_pcs = array2table(gene_pcs_arr, 'VariableNames',cellstr(arrayfun(@(x) sprintf('PC%d', x), 1:n_pcs, 'UniformOutput', false)));