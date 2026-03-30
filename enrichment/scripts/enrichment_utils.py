# from https://github.com/richardajdear

import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm


def fit_weights(gene_expr, pc_scores, n_components=3):
    """
    Get gene weights by correlating expression with scores
    """
    x = gene_expr.values
    y = pc_scores.values / np.std(pc_scores.values, axis=0)
    xv = x - x.mean(axis=0)
    yv = y - y.mean(axis=0)
    xvss = (xv * xv).sum(axis=0)
    yvss = (yv * yv).sum(axis=0)
    result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
    # bound the values to -1 to 1 in the event of precision issues
    result = np.maximum(np.minimum(result, 1.0), -1.0)

    weights = pd.DataFrame(
        result,
        index=gene_expr.columns,
        columns=["C" + str(i + 1) for i in range(n_components)],
    ).iloc[:, :n_components]

    return weights


def shuffle_gene_weights(weights, n=100):
    """
    Make null model by randomizing gene weights / ranks
    """
    # reshape to 2D if weights is a Series
    weights = weights.values if hasattr(weights, "values") else np.asarray(weights)
    if weights.ndim == 1:
        weights = weights[:, np.newaxis]

    n_components = weights.shape[1]
    null_weights = np.repeat(weights[:, :n_components, np.newaxis], n, axis=2)
    # null_weights = np.take_along_axis(null_weights, np.random.randn(*null_weights.shape).argsort(axis=0), axis=0)
    for c in range(n_components):
        for i in range(n):
            np.random.shuffle(null_weights[:, c, i])

    return null_weights


def match_genes(gene_labels, weights):
    """
    Make mask of which genes from each gene list are in the PCs
    """
    genes = weights.index
    gene_masks = {}
    gene_counts = {}
    for l in gene_labels["label"].unique():
        genes_in_label = gene_labels.query("label == @l")["gene"]
        matches = np.isin(genes, genes_in_label)
        if sum(matches) > 0:
            gene_masks[l] = matches

        gene_counts[l] = pd.Series(
            {"n_genes": len(genes_in_label), "n_matches": sum(matches)}
        )
    gene_counts = pd.concat(gene_counts).unstack()

    return gene_masks, gene_counts


def compute_enrichments(
    weights, null_weights, gene_labels, how="mean", norm=False, posneg=None
):
    """
    Compute scores for each gene label, either mean, or median rank
    """
    weights_array = weights.values if hasattr(weights, "values") else np.asarray(weights)
    if weights_array.ndim == 1:
        weights_array = weights_array[:, np.newaxis]

    n_components = weights_array.shape[1]
    axis_names = (
        list(weights.columns)
        if hasattr(weights, "columns")
        else [f"C{i + 1}" for i in range(n_components)]
    )
    gene_index = (
        weights.index
        if hasattr(weights, "index")
        else pd.RangeIndex(start=0, stop=weights_array.shape[0], step=1)
    )
    weights_df = pd.DataFrame(weights_array, index=gene_index, columns=axis_names)
    gene_masks, gene_counts = match_genes(gene_labels, weights_df)

    weights = weights_array.copy()
    nulls = null_weights.copy()
    # Take absolute values of standardized weights
    if norm:
        weights = StandardScaler().fit_transform(weights)
        for i in range(nulls.shape[2]):
            nulls[:, :, i] = StandardScaler().fit_transform(nulls[:, :, i])

    if posneg == "abs":
        weights = np.abs(weights)
        nulls = np.abs(nulls)
    elif posneg == "pos":
        weights = np.where(weights < 0, np.nan, weights)
        nulls = np.where(nulls < 0, np.nan, nulls)
    elif posneg == "neg":
        weights = np.where(weights > 0, np.nan, weights)
        nulls = np.where(nulls > 0, np.nan, nulls)

    true_enrichments = {}
    null_enrichments = {}

    for label, mask in gene_masks.items():
        if how == "mean":
            true_enrichments[label] = pd.Series(np.nanmean(weights[mask, :], axis=0))
            null_enrichments[label] = pd.DataFrame(
                np.nanmean(nulls[mask, :, :], axis=0)
            ).T
        elif how == "median":  #### not working
            true_ranks = weights.argsort(0).argsort(0)
            true_enrichments[label] = pd.Series(
                np.nanmedian(true_ranks[mask, :], axis=0)
            )
            null_enrichments[label] = pd.DataFrame(
                np.nanmedian(nulls[mask, :, :], axis=0)
            ).T

    true_enrichments = (
        pd.concat(true_enrichments).unstack(1).set_axis(axis_names, axis=1)
    )
    null_enrichments = (
        pd.concat(null_enrichments)
        .set_axis(axis_names, axis=1)
        .reset_index(level=0)
        .rename({"level_0": "label"}, axis=1)
    )

    return true_enrichments, null_enrichments, gene_counts


def compute_null_p(
    true_enrichments,
    null_enrichments,
    gene_counts=None,
    adjust="fdr_bh",
    adjust_by_label=False,
    order=None,
):
    """
    Compute null p values
    """
    null_pct = np.zeros(true_enrichments.shape)
    for m, label in enumerate(true_enrichments.index):
        for i in range(true_enrichments.shape[1]):
            nulls_ = null_enrichments.set_index("label").loc[label].iloc[:, i]
            true_ = true_enrichments.iloc[m, i]
            pct = percentileofscore(nulls_, true_) / 100
            null_pct[m, i] = pct

    true_mean = true_enrichments.stack().rename("true_mean")

    null_mean = (
        null_enrichments.groupby("label")
        .agg(["mean", "std"])
        .stack(0)
        .rename_axis([None, None])
        .set_axis(["null_mean", "null_std"], axis=1)
    )

    null_p = (
        pd.DataFrame(
            null_pct, index=true_enrichments.index, columns=true_enrichments.columns
        )
        .stack()
        .rename("pct")
        .to_frame()
        .join(true_mean)
        .join(null_mean)
        .assign(z=lambda x: (x["true_mean"] - x["null_mean"]) / x["null_std"])
        .assign(pos=lambda x: [pct > 0.5 for pct in x["pct"]])
        .assign(
            p=lambda x: [(1 - pct) * 2 if pct > 0.5 else pct * 2 for pct in x["pct"]]
        )  # x2 for two-sided
    )

    # Apply multiple comparisons
    if adjust is not None:
        # Adjust across axes only (not by label)?
        if adjust_by_label:
            null_p = (
                null_p.assign(
                    q=lambda x: x.groupby(level=0)
                    .apply(
                        lambda y: pd.Series(
                            multipletests(y["p"], method=adjust)[1], index=y.index
                        )
                    )
                    .reset_index(0, drop=True)  # make index match
                ).assign(sig=lambda x: x["q"] < 0.05)
                # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
            )
        else:
            null_p = (
                null_p.assign(
                    q=lambda x: multipletests(x["p"], method=adjust)[1]
                ).assign(sig=lambda x: x["q"] < 0.05)
                # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
            )
    else:
        null_p = null_p.assign(q=lambda x: x["p"]).assign(sig=lambda x: x["q"] < 0.05)

    null_p = null_p.reset_index().rename({"level_0": "label", "level_1": "C"}, axis=1)

    if gene_counts is not None:
        null_p = null_p.join(gene_counts, on="label")

    # Fix order of gene labels
    if order is None:
        order = true_enrichments.index
    null_p = null_p.assign(
        label=lambda x: pd.Categorical(x["label"], ordered=True, categories=order)
    ).sort_values("label")

    return null_p
