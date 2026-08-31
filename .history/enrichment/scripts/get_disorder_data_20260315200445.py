import pandas as pd

replace_dict = {
    "01-Mar": "MARCH1",
    "02-Mar": "MARCH2",
    "03-Mar": "MARCH3",
    "04-Mar": "MARCH4",
    "05-Mar": "MARCH5",
    "06-Mar": "MARCH6",
    "07-Mar": "MARCH7",
    "08-Mar": "MARCH8",
    "09-Mar": "MARCH9",
    "10-Mar": "MARCH10",
    "11-Mar": "MARCH11",
    "01-Sep": "SEPT1",
    "02-Sep": "SEPT2",
    "03-Sep": "SEPT3",
    "04-Sep": "SEPT4",
    "05-Sep": "SEPT5",
    "06-Sep": "SEPT6",
    "07-Sep": "SEPT7",
    "08-Sep": "SEPT8",
    "09-Sep": "SEPT9",
    "10-Sep": "SEPT10",
    "11-Sep": "SEPT11",
    "12-Sep": "SEPT12",
    "13-Sep": "SEPT13",
    "14-Sep": "SEPT14",
    "15-Sep": "SEPT15",
    "01-Dec": "DECR1",
    "02-Dec": "DECR2",
}


def get_deg_combined():
    deg_dict = {
        "SCZ--Gandal 2018": pd.read_csv("../data/deg/gandal_genes_rnaseq.csv").loc[
            lambda x: x["SCZ.fdr"] <= 0.05, ["gene_name", "SCZ.log2FC"]
        ],
        "SCZ--Fromer 2016": pd.read_csv(
            "../data/deg/fromer2016_tableS3.csv", header=1
        ).loc[lambda x: x["FDR estimate"] <= 0.05, ["Gene Symbol", "logFC"]],
        "SCZ--Collado-Torres 2019": pd.read_csv(
            "../data/deg/colladotorres2019_tableS11.csv"
        )
        .query("region=='DLPFC'")
        .loc[lambda x: x["adj.P.Val"] <= 0.05, ["Symbol", "logFC"]],
        "SCZ--Jaffe 2018": pd.read_csv(
            "../data/deg/jaffe2018_tableS9.csv", header=1
        ).loc[lambda x: x["fdr_qsva"] <= 0.05, ["Symbol", "log2FC_qsva"]],
    }
    deg_dict = {
        name: deg.set_axis(["gene", "log2FC"], axis=1) for name, deg in deg_dict.items()
    }

    deg_genes = (
        pd.concat(deg_dict)
        .reset_index(0)
        .rename({"level_0": "label"}, axis=1)
        .replace({"gene": replace_dict})
    )

    # drop logFC to cleanly remove duplicates
    deg_genes = deg_genes.drop("log2FC", axis=1)

    deg_genes = deg_genes.drop_duplicates().dropna()
    return deg_genes


def get_deg_consensus():
    deg_all = get_deg_combined()

    grouping = ["label", "gene"]
    agg_dict = {"study": "nunique"}

    deg_consensus = (
        deg_all.assign(
            study=lambda x: x["label"].str.split("--", expand=True)[1],
            label=lambda x: x["label"].str.split("--", expand=True)[0],
        )
        .groupby(grouping, as_index=False)
        .agg(agg_dict)
        .loc[
            lambda x: (x["study"] >= 2) | (x["label"] == "MDD")
        ]  # either 2+ studies, or MDD
    )

    return deg_consensus


def get_gwas_combined():
    # SCZ data from https://figshare.com/articles/dataset/scz2022/19426775?file=35775617
    trubetskoy = (
        pd.read_csv(f"../data/gwas/trubetskoy2022_extended.csv")
        .loc[lambda x: x["Extended.GWAS"] == "YES", "Symbol.ID"]
        #   .rename('gene')
        #   .loc[lambda x: x["Prioritised"]==1, :]
        #   .loc[lambda x: x["Extended.GWAS"]=='YES', 'Ensembl.ID']
        #   .pipe(ensembl_id_to_gene_symbol)
        .rename("gene")
        .reset_index(drop=True)
    )
    gwas_dict = {
        "SCZ": trubetskoy,
    }
    df = (
        pd.concat(gwas_dict)
        .reset_index(0)
        .rename({"level_0": "label"}, axis=1)
        .assign(
            gene=lambda x: x["gene"].str.replace("\\..*", "", regex=True)
        )  # drop variants
        .replace({"gene": replace_dict})
        .drop_duplicates()
        .dropna()
    )
    return df
