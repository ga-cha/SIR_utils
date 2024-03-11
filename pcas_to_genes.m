% given PCs of interest and pca matrix, returns top contributing genes


% something like take every gene and correlate with PC
% then take like, top 100 or above some arbitrary corr threshold

n_hits = 100;
corrs = corr(table2array(genes), table2array(gene_pcs(:,13)));
[max_corrs, max_i] = maxk(corrs, n_hits);
corr_genes = string(genes.Properties.VariableNames(max_i));
corr_genes = reshape(corr_genes, [n_hits,1]);
corr_genes = cat(2, corr_genes, max_corrs);