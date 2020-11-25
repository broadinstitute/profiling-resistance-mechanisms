#!/usr/bin/env python
# coding: utf-8

# ## Quantify signature performance
# 
# **Gregory Way, 2020**
# 
# Here, we collect the following performance metrics for the singscore predictions for each dataset:
# 
# 1. Accuracy
# 2. Average precision
# 
# We calculate these metrics for each model split (training, test, and validation) as well as stratified by the following:
# 
# * Total
# * Per plate
# * Per sample

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from sklearn.metrics import accuracy_score, average_precision_score

import plotnine as gg


# In[2]:


def get_metrics(df):
    acc = accuracy_score(
        df.Metadata_clone_type_indicator,
        df.y_pred
    )
    
    avg_prec = average_precision_score(
        df.Metadata_clone_type_indicator,
        df.y_pred,
        average="samples"
    )
    
    metric_dict = {
        "accuracy": acc,
        "avg_precision": avg_prec
    }
    
    return pd.Series(metric_dict)


def get_metric_pipeline(df, metric_comparisons, datasets, shuffle=False, threshold=0):
    metric_results = {}
    for metric_compare in metric_comparisons:
        metadata_groups = metric_comparisons[metric_compare]
        metric_results[metric_compare] = {}
        for dataset in datasets:
            result_subset_df = (   
                df
                .query("dataset == @dataset")
                .query("signature == @dataset")
                .assign(y_pred=0)
            )

            if shuffle:
                score = np.random.permutation(result_subset_df.TotalScore)
            else:
                score = result_subset_df.TotalScore

            result_subset_df.loc[score > threshold, "y_pred"] = 1

            # Now, get metrics
            metric_results[metric_compare][dataset] = (
                result_subset_df
                    .groupby(metadata_groups)
                    .apply(get_metrics)
                    .reset_index()
                    .melt(
                        id_vars=metadata_groups,
                        value_vars=["accuracy", "avg_precision"],
                        var_name="metric",
                        value_name="metric_value"
                    )
                    .assign(dataset=dataset, shuffle=shuffle)
            )

        # Combine results into metric specific dataframes
        metric_results[metric_compare] = (
            pd.concat(metric_results[metric_compare]).reset_index(drop=True)
        )

    return metric_results


# In[3]:


np.random.seed(5678)


# In[4]:


# Set constants
sig_dir = pathlib.Path("results", "singscore")
results_file = pathlib.Path(sig_dir, "full_bulk_signature_singscore_results.tsv.gz")

output_dir = pathlib.Path("results", "performance")

num_permutations = 25
threshold = 0
datasets = ["cloneAE", "cb5083", "ixazomib"]

metric_comparisons = {
    "total": ["Metadata_model_split"],
    "plate": ["Metadata_model_split", "Metadata_Plate"],
    "sample": ["Metadata_model_split", "Metadata_clone_number"]
}


# In[5]:


# Load data
sig_dir = pathlib.Path("results", "singscore")
results_file = pathlib.Path(sig_dir, "full_bulk_signature_singscore_results.tsv.gz")

output_dir = pathlib.Path("results", "performance")

results_df = pd.read_csv(results_file, sep="\t").query("Metadata_model_split != 'perturbation'")

print(results_df.shape)
results_df.head()


# ## Get performance metrics

# In[6]:


# Using real predictions
real_metric_results = get_metric_pipeline(
    results_df,
    metric_comparisons,
    datasets,
    shuffle=False,
    threshold=threshold
)


# In[7]:


# Using shuffled predictions
all_shuffle_results = {compare: [] for compare in metric_comparisons}
for i in range(0, num_permutations):
    np.random.seed(i)
    shuffle_metric_results = get_metric_pipeline(
        results_df,
        metric_comparisons,
        datasets,
        shuffle=True,
        threshold=threshold
    )
    for compare in metric_comparisons:
        metric_df = shuffle_metric_results[compare].assign(permutation=i)
        all_shuffle_results[compare].append(metric_df)


# ## Write to file

# In[8]:


for compare in metric_comparisons:
    full_results_df = real_metric_results[compare]
    shuffle_results_df = pd.concat(all_shuffle_results[compare]).reset_index(drop=True)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_metric_performance.tsv")
    full_results_df.to_csv(output_file, sep="\t", index=False)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_shuffle_metric_performance.tsv")
    shuffle_results_df.to_csv(output_file, sep="\t", index=False)

