import warnings
from pathlib import Path
import abagen
from abagen import io
from abagen.samples_ import _get_struct
from abagen.correct import keep_stable_genes

import pandas as pd
import numpy as np
import nibabel as nib


def get_labels_hcp():
    """
    Get HCP atlas labels from source file and format for ggseg
    """
    hcp_info = (
        pd.read_csv("../data/parcellations/HCP-MMP1_UniqueRegionList.txt")
        .rename(columns={"LR": "hemisphere", "regionID": "id", "region": "label"})
        .assign(structure="cortex")
        .assign(id=lambda x: [id - 20 if id > 180 else id for id in x["id"]])
        # Aurina's HCP images code regions as 1-360, not 1-180,201-380
        # so recode to match the image files
    )
    labels_hcp = hcp_info.set_index("id")["label"]
    return labels_hcp


def get_labels_s100ts2():
    """
    Get labels
    """
    s100ts2_info = pd.read_csv(
        "/fs03/kg98/gchan/Atlases/Tian/Schaefer_Tian/reordered/"
        "Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S2_label.csv"
    )
    labels_s100ts2 = s100ts2_info.set_index("id")["label"]

    return labels_s100ts2


def fetch_hcp(native=True, only_left=True):
    """
    Get HCP atlas
    """
    hcp_info = pd.read_csv("../data/parcellations/HCP-MMP1+TianS2_ordered.csv")

    path = "/fs03/kg98/gchan/Atlases/AHBA_individual/Schaefer_GC/hcp-mmp_xfm_AA"
    fname = "parcs/HCPMMP1_acpc_uncorr+TianS2.nii"
    donors = ["9861", "10021", "12876", "14380", "15496", "15697"]
    subjs = [
        "H0351_2001",
        "H0351_2002",
        "H0351_1009",
        "H0351_1012",
        "H0351_1015",
        "H0351_1016",
    ]

    if native:
        atlas_hcp = {
            "image": {
                donor: nib.load(f"{path}/{subj}/{fname}")
                for subj, donor in zip(subjs, donors)
            },
            "info": hcp_info,
        }
    else:
        atlas_hcp = {
            "image": nib.load("../data/parcellations/HCP-MMP_1mm.nii.gz"),
            "info": pd.read_csv("../data/parcellations/HCP-MMP1_UniqueRegionList.txt")
            .rename(columns={"LR": "hemisphere", "regionID": "id", "region": "label"})
            .assign(structure="cortex"),
        }

    def remove_right(img):
        img_cortex = np.where((img.dataobj[:] > 196), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)

    if only_left:
        for donor, image in atlas_hcp["image"].items():
            atlas_hcp["image"][donor] = remove_right(image)

    return atlas_hcp


def fetch_s100ts2(native=True, only_left=True):
    """
    Get Schaefer 100 + Tian S2 atlas
    """

    path = "/fs03/kg98/gchan/Atlases/AHBA_individual/Schaefer_GC/Schaefer_xfm_AA/"
    atlas_info = pd.read_csv(
        path + "Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S2_label.csv"
    )

    fname = "parcs/Schaefer100_7net+TianS2.nii"
    donors = ["9861", "10021", "12876", "14380", "15496", "15697"]
    subjs = [
        "H0351_2001",
        "H0351_2002",
        "H0351_1009",
        "H0351_1012",
        "H0351_1015",
        "H0351_1016",
    ]

    if native:
        atlas_s100ts2 = {
            "image": {
                donor: nib.load(f"{path}/{subj}/{fname}")
                for subj, donor in zip(subjs, donors)
            },
            "info": atlas_info,
        }
    else:
        assert (
            native
        ), "Non-native version of Schaefer100+TianS2 not currently implemented"

    def remove_right(img):
        img_cortex = np.where((img.dataobj[:] > 66), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)

    if only_left:
        for donor, image in atlas_s100ts2["image"].items():
            atlas_s100ts2["image"][donor] = remove_right(image)

    return atlas_s100ts2


# Wrapper function around abagen.allen.get_expression_data that adds gene and region filtering
def get_abagen_by_roi(
    atlas,
    data_dir="../data/abagen-data",
    save_name=None,
    donors_threshold=0,
    gene_stability_threshold=0.0,
    return_donors=False,
    return_counts=False,
    return_stability=False,
    return_report=False,
    verbose=1,
    lr_mirror="bidirectional",
    **kwargs,
):
    """
    Helper function to get filtered AHBA expression matrix from abagen

    Parameters
    ----------
    atlas : dict
        A dictionary including 'image' and 'info', i.e. in the format returned by
        abagen.get_desikan_killiany(). For details see abagen.get_expression_data.
    data_dir : os.PathLike, optional
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded. If not specified will use the current
        directory.
    save_name : str, optional
        If provided will save the final cleaned expression matrix to
        '{data_dir}/expression/{save_name}.csv'
    region_donors_threshold : int, optional
        Minimum number of donors (1-6) that must have a sampled matched to a region
        For example, 3 means keep only regions with samples from at least 3 of 6 donors
    gene_stability_threshold : float, optional
        Quantile threshold of differential stability to retain genes
        For example, 0.9 means keep only the top 10% most stable genes
    return_donors : bool
        Whether to return separate expression matrices for each donor
    return_counts : bool
        Whether to return table of samples matched per region per donor
    return_stability : bool
        Whether to return differential stability of each gene
        NB: gene stability is computed after filtering for regions
    return report : bool
        Whether to return report of processing for use in Methods
    lr_mirror : {None, 'bidirectional', 'leftright', 'rightleft'}, optional
        NB: Overwrite abagen default of None to 'bidirectional'
        Whether to mirror microarray expression samples across hemispheres to
        increase spatial coverage. Using 'bidirectional' will mirror samples
        across both hemispheres, 'leftright' will mirror samples in the left
        hemisphere to the right, and 'rightleft' will mirror the right to the
        left.

    Returns
    -------
    expression : (R, G) pandas.DataFrame
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors, where the index corresponds to the unique
        integer IDs of `atlas` and the columns are gene names. Cleaned
    counts : (R, D) pandas.DataFrame
        Number of samples assigned to each of `R` regions in `atlas` for each
        of `D` donors (if multiple donors were specified); only returned if
        ``return_counts=True``.
    stability : (G,) pandas.Series
        Differential stability (i.e. mean correlation between donor pairs)
        of each gene across retained regions; only returned if
        ``return_stability=True``.
    report : str
        Methods describing processing procedures implemented to generate
        `expression`, suitable to be used in a manuscript Methods section. Only
        returned if ``return_report=True``.
    """
    # Step 1: Use abagen to get expression data with provided atlas
    # return_counts=True to get sample counts per region
    # return_donors=True to keep expression separate by donors
    # Overwrite default for lr_mirror to 'bidirectional'
    # Suppress FutureWarning for frame.append method
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expression, counts, report = abagen.allen.get_expression_data(
            # expression, counts, report = abagen_allen_tweaked.get_expression_data(
            # expression = abagen.get_expression_data(
            atlas=atlas["image"],
            atlas_info=atlas["info"],
            data_dir=data_dir + "/microarray",
            n_proc=6,
            return_donors=True,
            return_counts=True,
            return_report=True,
            verbose=verbose,
            lr_mirror=lr_mirror,
            **kwargs,
        )
        # abagen return_donors returns a dict, but keep_stable_genes wants a list
        expression = list(expression.values())

    # Step 2: Identify regions to keep and filter each donor's expression matrix
    if donors_threshold is not None:
        # Find regions with at least 1 sample from at least X donors
        region_filter = (counts > 0).sum(axis=1) >= donors_threshold
        # Filter all donor expression matrices for those regions
        expression = [e.loc[region_filter, :] for e in expression]
        if verbose > 0:
            print(
                f"{region_filter.sum()} / {len(region_filter)} regions remain "
                f"after filtering for regions with samples from >= {donors_threshold} donors"
            )

    # Step 2: Identify stable genes and filter each donor's expression matrix
    if gene_stability_threshold is not None:
        # Remember gene labels to use as index
        gene_labels = expression[0].columns
        # Filter all donor expression matrices for stable genes
        expression, stability = keep_stable_genes(
            expression, threshold=gene_stability_threshold, return_stability=True
        )
        stability = pd.Series(stability, index=gene_labels)
        if verbose > 0:
            print(
                f"{expression[0].shape[1]} / {len(stability)} genes remain "
                f"after filtering for top {round(1-gene_stability_threshold,2)} stability"
            )

    # Step 3: Aggregate expression across donors
    if not return_donors:
        expression = pd.concat(expression).groupby("label").mean()

    # Optionally save combined expression data
    if save_name is not None:
        Path(data_dir + "/expression").mkdir(parents=True, exist_ok=True)
        if not return_donors:
            save_path = f"{data_dir}/expression/{save_name}.csv"
            expression.to_csv(save_path)
            print(f"Expression matrix saved to {save_path}")
        else:
            for i, donor_expression in enumerate(expression):
                save_path = f"{data_dir}/expression/{save_name}_{i}.csv"
                donor_expression.to_csv(save_path)
                print(f"Donor expression matrix {i} saved to {save_path}")

    # Pack outputs
    out = (expression,)
    if return_counts:
        out += (counts,)
    if return_stability:
        out += (stability,)
    if return_report:
        out += (report,)
    if len(out) == 1:
        out = out[0]

    return out
