import pandas as pd
import numpy as np
from brainsmash.mapgen.base import Base
import scipy.io
import sys

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python gene_brainsmash.py gene_id")
        sys.exit(1)
    i = int(sys.argv[1])
    input_file = sys.argv[2]

    genes = pd.read_csv(input_file)                                 # from abagen
    distmat = np.loadtxt('./data/distmat_ifod.csv', delimiter=',')  # ifod2 sift2
    gene_names = genes.columns.tolist()

    surrogates = pd.DataFrame(np.zeros(genes.shape), columns=gene_names)

    for g in range(genes.shape[1]):
        gene_map = genes.iloc[:, g].to_numpy()

        gen = Base(gene_map, distmat)
        surrogates[gene_names[g]] = gen(n=1)

    scipy.io.savemat(f'./results/gene_null/null_genes_{i}.mat', mdict={'null_genes': surrogates.to_numpy()})