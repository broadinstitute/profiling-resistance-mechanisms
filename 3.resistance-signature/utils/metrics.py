import numpy as np
import pandas as pd

from sklearn.metrics import accuracy_score, average_precision_score


def get_metrics(df):
    """A helper function to output various performance metrics

    Parameters
    ----------
    df : pandas.DataFrame
        a data frame storing metadata and predictions, must include two columns:
        ("Metadata_clone_type_indicator", and "y_pred")

    Returns
    -------
    dict
        A bundle of predefined performance metrics
    """
    acc = accuracy_score(df.Metadata_clone_type_indicator, df.y_pred)

    avg_prec = average_precision_score(
        df.Metadata_clone_type_indicator, df.y_pred, average="samples"
    )

    metric_dict = {"accuracy": acc, "avg_precision": avg_prec}

    return pd.Series(metric_dict)


def get_metric_pipeline(
    df, metric_comparisons, datasets, shuffle=False, threshold=0, signature=True
):
    """A helper function to output performance metric for a given dataframe and options

    Paramaters
    ----------
    df : pandas.DataFrame
        a data frame storing metadata and predictions, must include the columns:
        ("Metadata_clone_type_indicator", "y_pred", "signature", "TotalScore", "dataset")
    metric_comparisons : dict
        a dictionary with metadata splits indicating how to track performance
    datasets : list
        list of strings indicating which "dataset" to use in calculation
        (subsets dataset column)
    shuffle : bool, optional
        Whether or not to shuffle the actual signature scores before computing metrics.
        Defaults to False.
    threshold : float, optional
        How to distinguish positive from negative classes. Defaults to 0.
    signature : bool, optional
        In cases with multiple datasets and predictions made across datasets, only
        track performance for the plates used in the given dataset. Defaults to True.

    Returns
    -------
    dict
        All predefined performance metrics (as described in :py:func:`get_metrics`).
    """
    metric_results = {}
    for metric_compare in metric_comparisons:
        metadata_groups = metric_comparisons[metric_compare]
        metric_results[metric_compare] = {}
        for dataset in datasets:
            result_subset_df = df.query("dataset == @dataset").assign(y_pred=0)

            if signature:
                result_subset_df = result_subset_df.query("signature == @dataset")

            if shuffle:
                score = np.random.permutation(result_subset_df.TotalScore)
            else:
                score = result_subset_df.TotalScore

            result_subset_df.loc[score > threshold, "y_pred"] = 1

            # Now, get metrics
            metric_results[metric_compare][dataset] = (
                result_subset_df.groupby(metadata_groups)
                .apply(get_metrics)
                .reset_index()
                .melt(
                    id_vars=metadata_groups,
                    value_vars=["accuracy", "avg_precision"],
                    var_name="metric",
                    value_name="metric_value",
                )
                .assign(dataset=dataset, shuffle=shuffle)
            )

        # Combine results into metric specific dataframes
        metric_results[metric_compare] = pd.concat(
            metric_results[metric_compare]
        ).reset_index(drop=True)

    return metric_results
