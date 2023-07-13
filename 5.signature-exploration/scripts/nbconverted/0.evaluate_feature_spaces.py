#!/usr/bin/env python
# coding: utf-8

# # Evaluate feature spaces
# 
# Perform two analyses to evaluate four different feature spaces.
# 
# The two analyses are:
# 
# 1. UMAP
# 2. Clustering with Silhouette and enrichment evaluations
# 
# The feature spaces include:
# 
# 1. All features
# 2. Feature selected features by traditional methods
# 3. Bortezomib signature features
# 4. All features except bortezomib signature features

# In[1]:


import umap
import pathlib
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact

from typing import List, Union

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


np.random.seed(1234)


# In[3]:


def process_umap(
        data_df: pd.DataFrame,
        umap_category: str,
        features: Union[str, List[str]] = "infer"
        ) -> pd.DataFrame:
    """
    Apply UMAP to a dataframe and return a dataframe with UMAP coordinates and metadata.

    Parameters
    ----------
    data_df : pd.DataFrame
        The input dataframe. Should include feature columns and metadata columns.
    umap_category : str
        A label to add to the output dataframe under the column 'Metadata_umap_category'.
    features : str or list of str, default 'infer'
        The names of the feature columns to use for UMAP. If 'infer', all non-metadata columns are used.

    Returns
    -------
    pd.DataFrame
        A dataframe with the UMAP coordinates ('umap_0' and 'umap_1'), the metadata, and the UMAP category.

    Notes
    -----
    This function assumes that metadata columns can be identified by the function 'infer_cp_features' with 'metadata=True'.
    """
    # Prepare UMAP input by removing metadata columns
    metadata_cols = infer_cp_features(data_df, metadata=True)
    metadata_df = data_df.loc[:, metadata_cols]

    if features == "infer":
        umap_data_df = data_df.drop(metadata_cols, axis="columns")
    else:
        umap_data_df = data_df.loc[:, features]

    # Apply UMAP
    reducer = umap.UMAP(random_state=123)
    embedding = reducer.fit_transform(umap_data_df)
    
    # Setup plotting logic
    embedding_df = pd.DataFrame(embedding, columns=['umap_0', 'umap_1'])
    embedding_df = embedding_df.merge(metadata_df, left_index=True, right_index=True)
    embedding_df = embedding_df.assign(Metadata_umap_category=umap_category)
    
    return embedding_df


def perform_clustering(
        df: pd.DataFrame,
        features: Union[str, List[str]],
        pca_n_components: int,
        class_column: str,
        positive_class: Union[int, str],
        feature_category: str,
        low_k: int, 
        high_k: int
        ) -> pd.DataFrame:
    """
    Perform KMeans clustering on a subset of a dataframe and calculate:
        1) Silhouette width of clustering solution
        2) Enrichment of a class within clusters

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe. Should include feature columns and a class column.
    features : list of str
        The names of the feature columns to use for clustering.
    pca_n_components : int
        The number of principal coponents to use for PCA transformation
    class_column : str
        The name of the class column.
    positive_class : int or str
        The value in the class column that represents the "positive" class.
    all_features : str
        A label to add to the results dataframe.
    low_k : int
        The lower bound for the range of k for KMeans clustering.
    high_k : int
        The upper bound for the range of k for KMeans clustering.

    Returns
    -------
    pd.DataFrame
        A dataframe with the results.
        Columns are 'k', 'silhouette_width', 'average_enrichment', 'average_p_value', and 'feature_category'.

    """
    # Prepare output dataframe
    results = pd.DataFrame(
        columns=[
            "k",
            "silhouette_width",
            "average_enrichment",
            "maximum_enrichment",
            "average_p_value",
            "feature_category"
        ]
    )

    # Subset dataframe
    if features == "infer":
        metadata_cols = infer_cp_features(df, metadata=True)
        subset_df = df.drop(metadata_cols, axis="columns")
    else:
        subset_df = df.loc[:, features]

    # Transform data into PCA space to account for differential feature numbers influencing distance metrics
    pca = PCA(n_components=pca_n_components)
    subset_df = pca.fit_transform(subset_df)

    # Get boolean series indicating whether each sample belongs to the positive class
    is_in_class = df.loc[:, class_column] == positive_class

    # Perform kmeans clustering and calculate metrics
    for k in range(low_k, high_k+1):
        # KMeans clustering
        kmeans = KMeans(n_clusters=k, random_state=0, n_init=20).fit(subset_df)

        # Calculate silhouette width
        silhouette_width = silhouette_score(subset_df, kmeans.labels_)

        # Calculate enrichment for each cluster
        enrichments = []
        p_values = []
        for cluster_id in np.unique(kmeans.labels_):
            # Get boolean series indicating whether each sample is in the current cluster
            is_in_cluster = kmeans.labels_ == cluster_id

            # Calculate contingency table
            table = np.zeros((2, 2))
            table[0, 0] = np.sum(is_in_class & is_in_cluster)   # in class and in cluster
            table[0, 1] = np.sum(is_in_class & ~is_in_cluster)  # in class and not in cluster
            table[1, 0] = np.sum(~is_in_class & is_in_cluster)  # not in class and in cluster
            table[1, 1] = np.sum(~is_in_class & ~is_in_cluster) # not in class and not in cluster

            # Calculate Fisher's exact test
            odds_ratio, p_value = fisher_exact(table)
            enrichments.append(odds_ratio)
            p_values.append(p_value)

        # Calculate average enrichment
        average_enrichment = sum(enrichments) / len(enrichments)
        maximum_enrichment = np.max(enrichments)
        average_p_value = sum(p_values) / len(p_values)

       # Append results
        results_row = pd.DataFrame({
            "k": [k],
            "silhouette_width": [silhouette_width],
            "average_enrichment": [average_enrichment],
            "maximum_enrichment": [maximum_enrichment],
            "average_p_value": [average_p_value],
            "feature_category": [feature_category]
        })
        results = pd.concat([results, results_row], ignore_index=True)

    return results.reset_index(drop=True)


# ## Define paths

# In[4]:


data_dir = pathlib.Path("..", "2.describe-data", "data", "merged")
signature_dir = pathlib.Path("..", "3.resistance-signature")

profile_file = pathlib.Path(f"{data_dir}/all_merged_profiles_before_feature_selection.csv.gz")
profile_feature_selection_file = pathlib.Path(f"{data_dir}/all_merged_profiles.csv.gz")

bz_signature_file = pathlib.Path(f"{signature_dir}/results/signatures/signature_summary_bortezomib_signature.tsv.gz")

output_umap_file = pathlib.Path("results", "umap_feature_summary.tsv.gz")
output_cluster_file = pathlib.Path("results", "clustering_feature_summary.tsv.gz")


# ## Load profiles

# In[5]:


profile_df = pd.read_csv(profile_file, low_memory=False)

print(profile_df.shape)
profile_df.head()


# In[6]:


# Load feature selected dataframe
profile_feature_select_df = pd.read_csv(profile_feature_selection_file)

print(profile_feature_select_df.shape)
profile_feature_select_df.head()


# In[7]:


# Load bortezomib signature features
bz_sig_df = pd.read_csv(bz_signature_file, sep="\t")

bz_sig_features = bz_sig_df.query("final_signature").features.to_list()

print(bz_sig_df.shape)
print(len(bz_sig_features))
bz_sig_df.head()


# In[8]:


# Get all features except those in the signature
all_features_except_bz_sig_features = infer_cp_features(profile_df)
all_features_except_bz_sig_features = [x for x in all_features_except_bz_sig_features if x not in bz_sig_features]
len(all_features_except_bz_sig_features)


# ## Step 1: Calculate UMAP coordinates using four different feature spaces:

# In[9]:


# 1) All feature umap
umap_all_feature_df = process_umap(
    profile_df,
    umap_category="all_features"
)

# 2) Feature selected umap
umap_fs_feature_df = process_umap(
    profile_feature_select_df,
    umap_category="feature_selected"
)

# 3) Bortezomib signature features
umap_bz_sig_df = process_umap(
    profile_df,
    features=bz_sig_features,
    umap_category="bortezomib_signature_features"
)

# 4) All features except bortezomib signature
umap_non_bz_sig_df = process_umap(
    profile_df,
    features=all_features_except_bz_sig_features,
    umap_category="all_except_bortezomib_signature_features"
)


# In[10]:


# Output umap summary
umap_df = pd.concat(
    [
        umap_all_feature_df,
        umap_fs_feature_df,
        umap_bz_sig_df,
        umap_non_bz_sig_df
    ]
).reset_index(drop=True)

umap_df.to_csv(output_umap_file, sep="\t", index=False)

print(umap_df.shape)
umap_df.head()


# ## Step 2: Perform clustering and calculate enrichment and silhouette scores for the four feature spaces

# In[11]:


low_k = 2
high_k = 14
pca_n_components = 30

# 1) All feature clustering
all_feature_cluster_df = perform_clustering(
    df = profile_df,
    features = "infer",
    pca_n_components=pca_n_components,
    class_column = "Metadata_clone_type",
    positive_class = "wildtype",
    feature_category = "all_features",
    low_k = low_k, 
    high_k = high_k
)

# 2) Feature selected clustering
fs_feature_cluster_df = perform_clustering(
    df = profile_feature_select_df,
    features = "infer",
    pca_n_components=pca_n_components,
    class_column = "Metadata_clone_type",
    positive_class = "wildtype",
    feature_category = "feature_selected",
    low_k = low_k, 
    high_k = high_k
)

# 3) Bortezomib signature features clustering
bz_feature_cluster_df = perform_clustering(
    df = profile_df,
    features = bz_sig_features,
    pca_n_components=pca_n_components,
    class_column = "Metadata_clone_type",
    positive_class = "wildtype",
    feature_category = "bortezomib_signature_features",
    low_k = low_k, 
    high_k = high_k
)

# 4) All features except bortezomib signature clustering
non_bz_feature_cluster_df = perform_clustering(
    df = profile_df,
    features = all_features_except_bz_sig_features,
    pca_n_components=pca_n_components,
    class_column = "Metadata_clone_type",
    positive_class = "wildtype",
    feature_category = "all_except_bortezomib_signature_features",
    low_k = low_k, 
    high_k = high_k
)


# In[12]:


# Output clustering summary
clustering_df = pd.concat(
    [
        all_feature_cluster_df,
        fs_feature_cluster_df,
        bz_feature_cluster_df,
        non_bz_feature_cluster_df
    ]
).reset_index(drop=True)

clustering_df.to_csv(output_cluster_file, sep="\t", index=False)

print(clustering_df.shape)
clustering_df.head()

