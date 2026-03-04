from typing import Optional
import pandas as pd
import numpy as np
import mygene


def find_aliases(gene_df: pd.DataFrame, query: str = "abagen symbol") -> pd.DataFrame:
    """
    finds possible aliases for genes names in gene_df using mygene
    returns gene_df with additional columns: entrezgene, ensembl, symbol
    """
    mg = mygene.MyGeneInfo()
    mg_df = mg.querymany(
        gene_df[query],
        scopes="symbol, alias",
        fields="entrezgene, symbol, ensembl.gene",
        species="human",
        as_dataframe=True,
    )
    mg_df["ensembl"] = mg_df.apply(
        lambda row: (
            row["ensembl.gene"]
            if pd.notna(row["ensembl.gene"])
            else (
                row["ensembl"][0]["gene"]
                if isinstance(row["ensembl"], list) and len(row["ensembl"]) > 0
                else None
            )
        ),
        axis=1,
    )
    mg_df = mg_df[["entrezgene", "symbol", "ensembl"]]

    # remove duplicates
    dup = mg_df.index.duplicated(keep=False)  # find all duplicate indices
    mask = (~dup) | (  # keep all non-duplicates
        (dup) & (mg_df["symbol"] == mg_df.index)
    )  # keep duplicate where symbol matches index
    mg_df = mg_df[mask]

    # merge aliases into gene_df
    og_cols = gene_df.columns.tolist()
    gene_df = gene_df.merge(mg_df, left_on=query, right_index=True, how="left")
    gene_df = gene_df[["entrezgene", "ensembl", "symbol"] + og_cols]
    return gene_df


def pool_genes(gene_df: pd.DataFrame, gene_type: str) -> pd.DataFrame:
    """
    Pools gene_df for each gene in gene_type
    gene_type should be either 'risk gene' or 'clearance gene'
    """

    def get_corr(group):
        # group["correlation"] = group["correlation"] * group["pass_all"]
        group.loc[group["correlation"] == 0, "correlation"] = np.nan
        group["correlation"] = np.arctanh(group["correlation"])  # Fisher z-transform
        mean = group["correlation"].mean()
        return mean if not np.isnan(mean) else 0

    def get_seed(group):
        if group["correlation"].isna().all():
            return np.nan
        else:
            return group.loc[group["correlation"].idxmax(), "seed"]

    def normalise(data):
        return (data - data.min()) / (data.max() - data.min())

    pooled_df = (
        gene_df.groupby(gene_type)
        .apply(
            lambda group: pd.Series(
                {
                    "pass_all": group["pass_all"].any().astype(int),
                    "num_pass": group["pass_all"].sum(),
                    "correlation": get_corr(group),
                    "seed": get_seed(group),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )
    # reverse fisher z-transform
    pooled_df["correlation"] = np.tanh(pooled_df["correlation"])

    # normalise data between 0 and 1
    pooled_df.loc[pooled_df["correlation"] != 0, "correlation"] = normalise(
        pooled_df.loc[pooled_df["correlation"] != 0, "correlation"]
    )
    pooled_df["num_pass"] = normalise(pooled_df["num_pass"])

    pooled_df["score"] = pooled_df["correlation"] * pooled_df["num_pass"]
    pooled_df.rename(columns={gene_type: "abagen symbol"}, inplace=True)
    return pooled_df


if __name__ == "__main__":
    # TODO: basic functionality here
    pass
