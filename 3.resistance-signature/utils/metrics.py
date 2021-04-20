import numpy as np
import pandas as pd

from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    roc_curve,
    roc_auc_score,
)


def apply_shuffle(score):
    """A helper function to randomly permute scores

    Parameters
    ----------
    score : list
        A list like object that will be randomly shuffled

    Returns
    -------
    A randomly shuffled list
    """

    shuff_score = np.random.permutation(score)
    return shuff_score


def assign_pred_score(df, threshold=0, shuffle=False):
    """A helper function to determine status of input samples

    Parameters
    ----------
    df : pandas.DataFrame
        a data frame storing metadata and predictions, must include two columns:
        ("Metadata_clone_type_indicator", and "y_pred")
    threshold : float, optional
        How to distinguish positive from negative classes. Defaults to 0.
    shuffle : bool, optional
        Whether or not to shuffle the actual signature scores before computing metrics.
        Defaults to False.

    Returns
    -------
    df
        A copy of the input data frame but with a y_pred column
    """

    df = df.copy().assign(y_pred=0)

    if shuffle:
        score = apply_shuffle(df.TotalScore)
    else:
        score = df.TotalScore

    df.loc[score > threshold, "y_pred"] = 1
    return df


def get_metrics(df, return_roc_curve=False, threshold=0, shuffle=False):
    """A helper function to output various performance metrics

    Parameters
    ----------
    df : pandas.DataFrame
        a data frame storing metadata and predictions, must include two columns:
        ("Metadata_clone_type_indicator", and "y_pred")
    return_roc_curve : bool, optional
        determine if you should output the full roc curve info only, defaults to False
    threshold : float, optional
        How to distinguish positive from negative classes. Defaults to 0.
    shuffle : bool, optional
        Whether or not to shuffle the actual signature scores before computing metrics.
        Defaults to False.

    Returns
    -------
    dict
        A bundle of predefined performance metrics
    """
    if return_roc_curve:
        if shuffle:
            score = apply_shuffle(df.TotalScore)
        else:
            score = df.TotalScore

        roc_curve_results = roc_curve(
            y_true=df.Metadata_clone_type_indicator,
            y_score=score,
            drop_intermediate=False,
        )

        roc_curve_results = pd.DataFrame(roc_curve_results).transpose()
        roc_curve_results.columns = ["fpr", "tpr", "treshold"]

        df = assign_pred_score(df, threshold=threshold, shuffle=shuffle)
        roc_auc = roc_auc_score(
            y_true=df.Metadata_clone_type_indicator, y_score=df.y_pred
        )

        return roc_auc, roc_curve_results

    else:
        acc = accuracy_score(df.Metadata_clone_type_indicator, df.y_pred)

        avg_prec = average_precision_score(
            df.Metadata_clone_type_indicator, df.y_pred, average="samples",
        )

        metric_results = {"accuracy": acc, "avg_precision": avg_prec}

        return pd.Series(metric_results)


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
            result_subset_df = df.query("dataset == @dataset")

            if signature:
                result_subset_df = result_subset_df.query("signature == @dataset")

            result_subset_df = assign_pred_score(
                result_subset_df, threshold=threshold, shuffle=shuffle
            )

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
