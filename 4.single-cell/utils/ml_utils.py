import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    average_precision_score
)


def get_threshold_metrics(y_true, y_pred, drop_intermediate=False):
    """
    Retrieve true/false positive rates and auroc/aupr for class predictions
    Arguments:
    y_true - an array of gold standard mutation status
    y_pred - an array of predicted mutation status
    disease - a string that includes the corresponding TCGA study acronym
    Output:
    dict of AUROC, AUPR, pandas dataframes of ROC and PR data, and cancer-type
    """

    roc_columns = ['fpr', 'tpr', 'threshold']
    roc_items = zip(
        roc_columns,
        roc_curve(y_true, y_pred, drop_intermediate=drop_intermediate)
     )
    roc_df = pd.DataFrame.from_dict(dict(roc_items))

    pr_columns = ['precision', 'recall', 'threshold']
    prec, rec, thresh = precision_recall_curve(y_true, y_pred)
    pr_df = pd.DataFrame.from_records([prec, rec]).transpose()
    pr_df = pd.concat([pr_df, pd.Series(thresh)], ignore_index=True, axis=1)
    pr_df.columns = pr_columns

    auroc = roc_auc_score(y_true, y_pred, average='weighted')
    avg_precision = average_precision_score(y_true, y_pred, average='weighted')

    metric_dict = {
        'auroc': auroc,
        'average_precision': avg_precision,
        'roc_df': roc_df,
        'pr_df': pr_df,
    }
    return metric_dict


def shuffle_columns(feature):
    """
    To be used in an `apply` pandas func to shuffle columns around a datafame
    Import only
    """
    return np.random.permutation(feature.tolist())


def cross_validation_performance(trained_pipeline, output_file):
    # Cross-validated performance heatmap
    cv_results = (
        pd.concat([
            pd.DataFrame(trained_pipeline.cv_results_).drop('params', axis=1),
            pd.DataFrame.from_records(trained_pipeline.cv_results_['params'])
        ], axis=1)
    )

    cv_score_mat = pd.pivot_table(
        cv_results,
        values='mean_test_score',
        index='classify__l1_ratio',
        columns='classify__C'
    )

    ax = sns.heatmap(cv_score_mat, annot=True, fmt='.1%')
    ax.set_xlabel('Regularization strength multiplier (C)')
    ax.set_ylabel('L1 Ratio')
    ax.set_title("Multiclass model predictions")
    plt.tight_layout()
    plt.savefig(cv_heatmap_file, dpi=600, bbox_inches='tight')

    return cv_score_mat


def output_coefficients(trained_pipeline, x_df, column_recode, output_file):
    coef_df = pd.DataFrame(
        trained_pipeline
        .best_estimator_
        .named_steps['classify']
        .coef_
    ).transpose()

    coef_df = coef_df.rename(column_recode, axis="columns")
    coef_df.index = x_df.columns
    coef_df.index.name = 'feature'

    coef_df.to_csv(output_file, index=True, sep="\t")
    return coef_df


def model_apply(model, x_df, meta_df, y_recode, data_fit, shuffled, predict_proba=True):
    if predict_proba:
        scores_df = model.predict_proba(x_df)
    else:
        scores_df = model.predict(x_df)
    scores_df = (
        pd.DataFrame(
            scores_df, index=x_df.index
        )
        .rename(y_recode, axis="columns")
        .merge(meta_df, left_index=True, right_index=True)
    ).assign(data_fit=data_fit, shuffled=shuffled)
    return scores_df
