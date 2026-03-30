import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


# Step 1: Calculate the Enrichment Score (ES)
def calculate_es(gene_list, gene_set):
    N = len(gene_list)
    Nh = len(gene_set)
    hits = np.isin(gene_list, list(gene_set)).astype(int)
    no_hits = 1 - hits

    # Calculate the running sum
    running_sum = np.cumsum(hits / Nh - no_hits / (N - Nh))
    es = running_sum.max()
    return es


# Step 2: Estimate the Statistical Significance (Nominal P value)
def permute_and_calculate_es(gene_list, gene_set, n_permutations=1000, seed=None):
    if seed is not None:
        np.random.seed(seed)
    permuted_es = np.zeros(n_permutations)
    for i in tqdm(range(n_permutations)):
        permuted_list = np.random.permutation(gene_list)
        permuted_es[i] = calculate_es(permuted_list, gene_set)
    return permuted_es


def run_gsea(gene_list, gene_set, n_permutations=1000):
    observed_es = calculate_es(gene_list, gene_set)
    permuted_es = permute_and_calculate_es(gene_list, gene_set, n_permutations)
    p_value = np.mean(permuted_es >= observed_es)

    print(f"Observed ES: {observed_es}")
    print(f"Empirical P value: {p_value}")

    return observed_es, permuted_es, p_value


# New: run GSEA across labels and components and return DataFrame like compute_null_p
def run_gsea_df(
    weights_df: pd.DataFrame,
    gene_labels: pd.DataFrame,
    n_permutations: int = 1000,
    adjust: str = "fdr_bh",
    seed: int = None,
):
    """
    weights_df: DataFrame indexed by gene names, columns = components (e.g., C1, C2, ...)
    gene_labels: DataFrame with columns ['label','gene'] specifying gene sets
    n_permutations: number of null permutations to run per label-component
    seed: optional int seed for reproducibility (varied per label/component)
    returns: DataFrame similar in form to compute_null_p with columns including
             ['label','C','pct','true_mean','null_mean','null_std','z','pos','p','q','sig']
    """
    genes = list(weights_df.index)
    comps = list(weights_df.columns)

    rows = []
    for label in gene_labels["label"].unique():
        gene_set = set(gene_labels.query("label == @label")["gene"].values)
        if len(gene_set) == 0:
            continue
        for c_idx, c in enumerate(comps):
            scores = weights_df[c].values
            order = np.array(genes)[np.argsort(-scores)]
            obs_es = calculate_es(order, gene_set)

            # create a varied seed per label/component if seed provided
            seed_variant = None
            if seed is not None:
                seed_variant = int((seed + hash((label, c_idx))) % (2**32))

            # compute null ES distribution via permuting the ordered gene list
            null_es = permute_and_calculate_es(
                order, gene_set, n_permutations=n_permutations, seed=seed_variant
            )

            pct = percentileofscore(null_es, obs_es) / 100.0
            null_mean = null_es.mean()
            null_std = null_es.std(ddof=0) if null_es.std(ddof=0) > 0 else np.nan
            z = (
                (obs_es - null_mean) / null_std
                if not np.isnan(null_std) and null_std != 0
                else np.nan
            )
            pos = pct > 0.5
            p = (1 - pct) * 2 if pct > 0.5 else pct * 2

            rows.append(
                {
                    "label": label,
                    "C": c,
                    "pct": pct,
                    "true_mean": obs_es,
                    "null_mean": null_mean,
                    "null_std": null_std,
                    "z": z,
                    "pos": pos,
                    "p": p,
                }
            )

    df = pd.DataFrame(rows)

    # adjust p -> q and sig
    if df.shape[0] > 0:
        q = multipletests(df["p"].values, method=adjust)[1]
        df["q"] = q
        df["sig"] = df["q"] < 0.05
    else:
        df["q"] = []
        df["sig"] = []

    # preserve ordering similar to compute_null_p
    df["label"] = pd.Categorical(
        df["label"], categories=gene_labels["label"].unique(), ordered=True
    )
    df = df.sort_values("label").reset_index(drop=True)

    return df
