# from https://github.com/richardajdear

"""
Version class for analyzing gradients of AHBA
"""
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSCanonical, PLSRegression
from sklearn.preprocessing import StandardScaler
import brainspace
from brainsmash.mapgen.base import Base
from neuromaps.images import annot_to_gifti
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels


### Monkey Patch compute_affinity function to not zero-out negative values
def compute_affinity_new(
    x, kernel=None, sparsity=0.9, pre_sparsify=True, non_negative=False, gamma=None
):

    # run original function with different output
    return brainspace.gradient.kernels.compute_affinity(
        x,
        kernel=kernel,
        sparsity=sparsity,
        pre_sparsify=pre_sparsify,
        non_negative=False,
        gamma=gamma,
    )


# NB: must patch in 'gradient.kernels' namespace, not 'gradient.gradient'
brainspace.gradient.gradient.compute_affinity = compute_affinity_new


class gradientVersion:

    def __init__(
        self,
        approach="dm",
        n_components=5,
        sparsity=0,
        kernel=None,
        #  marker_genes=['NEFL', 'LGALS1', 'SYT6'],
        marker_genes=["NEFL", "LGALS1", "RFTN1"],
        random_state=0,
        **kwargs,
    ):
        """
        Initialize
        """
        self.n_components = n_components
        self.marker_genes = marker_genes

        self.approach = approach
        self.sparsity = sparsity

        if approach == "dm":
            kwargs["alpha"] = kwargs.get(
                "alpha", 1
            )  # set alpha=1 as default, but only for approach = 'dm'
            kernel = "normalized_angle" if kernel is None else kernel
        self.params = kwargs  # set embedding-specific parameters
        self.kernel = kernel

        self.expression = None
        self.scores = None
        self.gradients = brainspace.gradient.gradient.GradientMaps(
            n_components=n_components,
            approach=approach,
            kernel=kernel,
            random_state=random_state,
        )

    def fit(
        self,
        expression,
        scale=False,
        message=True,
        data_dir="../data/abagen-data/expression/",
    ):
        """
        Fit to data
        Marker genes is a list of genes to define the gradient direction, if gradient n is inversely aligned to marker n, the gradient will be flipped
        """
        # If expression is given as a string name, read the file
        if isinstance(expression, str):
            X = pd.read_csv(data_dir + expression + ".csv", index_col=0)
        else:
            X = expression
            # expression = name # for printing output

        # Clean data: drop regions with all nulls, and genes with any nulls
        X = X.dropna(axis=0, how="all").dropna(axis=1, how="any")
        # Optionally z-score expression
        if scale:
            X = X.apply(lambda x: (x - np.mean(x)) / np.std(x))
        self.expression = X

        # Fit gradients
        self.gradients.fit(X.values, sparsity=self.sparsity, **self.params)

        # Align gradients to marker genes
        # if no marker genes, results will be arbitrarily oriented but valid
        scores = pd.DataFrame(self.gradients.gradients_, index=X.index)
        for i, marker in enumerate(self.marker_genes):
            if marker in self.expression.columns:
                r_ = X.loc[:, marker].corr(scores.loc[:, i])
                if r_ < 0:
                    scores.loc[:, i] *= -1

        self.scores = scores
        self.affinity = compute_affinity_new(
            x=X.values, kernel=self.kernel, sparsity=self.sparsity
        )
        self.eigenvalues = self.gradients.lambdas_
        self.weights = self.fit_weights(n_components=self.n_components)

        if message:
            data_name = expression if isinstance(expression, str) else "(data given)"
            print(
                f"New gradients version: method={self.approach}, kernel={self.kernel}, sparsity={self.sparsity}, data={data_name}"
            )

        return self

    def clean_scores(self, scores=None, norm=True, n_components=3):
        """
        Normalize C1-3 scores, add labels
        """
        if scores is None:
            scores = self.scores

        if scores.shape[0] >= 120:
            labels = get_labels_hcp()
        elif scores.shape[0] <= 34:
            labels = get_labels_dk()
        else:
            labels = get_labels_dx().drop_duplicates()

        scores = (
            scores.iloc[:, :n_components]
            .set_axis(["C" + str(i + 1) for i in range(n_components)], axis=1)
            .rename_axis("id")
        )

        if norm:
            scores = scores.apply(lambda x: (x - np.mean(x)) / np.std(x))

        return scores.join(labels)

    def fit_weights(self, sort=False, n_components=3):
        """
        Get gene weights by correlating expression with scores
        """
        x = self.expression.values
        # y = self.scores.values
        y = self.scores.values / np.std(self.scores.values, axis=0)
        xv = x - x.mean(axis=0)
        yv = y - y.mean(axis=0)
        xvss = (xv * xv).sum(axis=0)
        yvss = (yv * yv).sum(axis=0)
        result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
        # bound the values to -1 to 1 in the event of precision issues
        result = np.maximum(np.minimum(result, 1.0), -1.0)

        weights = pd.DataFrame(
            result,
            index=self.expression.columns,
            columns=["C" + str(i + 1) for i in range(n_components)],
        ).iloc[:, :n_components]

        # Output sorted lists, or unsorted dataframe
        if sort:
            return self.sort_weights(weights)
        else:
            return weights

    def sort_weights(self, gene_weights):
        """
        Return genes with each column independently sorted by weight
        """
        gene_ranks = {}
        for g in range(gene_weights.shape[1]):
            gene_ranks[g] = (
                gene_weights.iloc[:, g]
                .sort_values(ascending=False)
                .reset_index()
                .set_axis(["gene", "weight"], axis=1)
            )

        return pd.concat(gene_ranks, axis=1)

    def match_components(self, df_corr, n_components=5):
        """
        Matching logic for a correlation df
        """
        _matches = [None] * n_components
        _corrs = [None] * n_components

        for i in range(n_components):
            # Find columns and rows already matched
            xs = [m[0] for m in _matches if m != None]
            ys = [m[1] for m in _matches if m != None]
            # Take correlation matrix and remove columns and rows already matched
            df_remain = df_corr.drop(xs, axis=0).drop(ys, axis=1)
            # Find the next highest correlation across the remaining values
            new_match = df_remain.stack().index[np.argmax(np.abs(df_remain.values))]
            # Add the new match into the list of matches
            _matches[i] = new_match
            # Add the new correlation into the list of correlations
            _corrs[i] = df_corr.loc[new_match]

        return _matches, _corrs

    def match_and_sort(self, other, df_corr):
        """
        Do correlation matching for genes
        """
        matches, corrs = self.match_components(df_corr)

        xs = [m[0] for m in matches]
        ys = [m[1] for m in matches]
        x_vars = [self.eigenvalues[i] for i in xs]
        y_vars = [other.eigenvalues[i] for i in ys]
        mean_vars = [(x + y) / 2 for x, y in zip(x_vars, y_vars)]
        sort_idx = np.argsort(mean_vars)[::-1]

        vars_sort = [mean_vars[i] for i in sort_idx]
        corrs_sort = [corrs[i] for i in sort_idx]
        matches_sort = [matches[i] for i in sort_idx]

        return (
            pd.DataFrame([matches_sort, corrs_sort, vars_sort])
            .T.set_axis(["match", "corr", "var_expl"], axis=1)
            .astype(dtype={"corr": float, "var_expl": float})
        )

    def corr_scores(self, other, match=False):
        """
        Correlate region scores with another version, with matching
        """

        df_corr = pd.concat([self.scores, other.scores], axis=1).corr().iloc[:5, 5:]

        if match:
            return self.match_and_sort(other, df_corr)
        else:
            return df_corr

    def corr_weights(self, other, match=False):
        """
        Correlate gene weights with another version
        Optionally with matching
        """

        df_corr = pd.concat([self.weights, other.weights], axis=1).corr().iloc[:5, 5:]

        if match:
            return self.match_and_sort(other, df_corr)
        else:
            return df_corr

    def score_from(self, other, clean=True):
        """
        Score from other weights
        """
        other_weights = other.weights.set_axis(range(other.weights.shape[1]), axis=1)
        gene_intersection = np.intersect1d(self.expression.columns, other_weights.index)
        _expression = self.expression.loc[:, gene_intersection]
        _weights = other_weights.loc[gene_intersection, :]

        scores = _expression @ _weights

        if clean:
            scores = self.clean_scores(scores=scores)

        return scores

    def fill_missing_scores(self, other, clean=True):
        """
        Fill in NA regions in this set of scores using scores from other version
        """
        filled_scores = self.scores.reindex(range(1, 181)).fillna(other.scores).dropna()

        if clean:
            filled_scores = self.clean_scores(scores=filled_scores)

        return filled_scores
