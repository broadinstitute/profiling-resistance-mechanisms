#!/usr/bin/env python
# coding: utf-8

# # Exploring Batch 3 Replicate Correlation
# 
# **Gregory Way, 2019**

# In[1]:


import os
import numpy as np
import pandas as pd
import plotnine as gg


# In[2]:


def get_pairwise_corr(x_df):
    x_cor_df = x_df.loc[:, ~x_df.columns.str.startswith("Metadata")].transpose().corr()
    x_cor_df = x_cor_df.where(pd.np.tril(pd.np.ones(x_cor_df.shape), k=-1).astype(bool))
    x_cor_df.index = x_df.Metadata_clone_number
    x_cor_df.columns = x_df.Metadata_clone_number
    x_cor_df.index.name = "{}_b".format(x_cor_df.index.name)
    x_cor_df = x_cor_df.stack().reset_index().rename({0: "correlation"}, axis='columns')
    x_cor_df = x_cor_df.assign(replicate=x_cor_df.iloc[:, 0] == x_cor_df.iloc[:, 1])
    return x_cor_df

def plot_replicate_corr(x_df):
    plot_obj = gg.ggplot(x_df, gg.aes(x="correlation", fill="replicate")) +         gg.geom_density(alpha = 0.3) +         gg.scale_fill_manual(name="Replicate",
                             labels={"True": "True",
                                     "False": "False"},
                             values=["#B99638", "#2DB898"]) + \
        gg.xlab("Pearson Correlation") + \
        gg.ylab("Density") + \
        gg.theme_light()
    return plot_obj

def plot_individual_replicate_corr(x_df):
    # Process input data to ensure all sample correlations are captured in
    # one column, but that replicate correlation is not captured
    left_df = (
        x_df
        .loc[:, ["Metadata_clone_number_b", "correlation", "replicate"]]
        .reset_index(drop=True)
        .rename({"Metadata_clone_number_b": "Metadata_clone_number"}, axis='columns')
    )
    right_df = (
        x_df
        .loc[:, ["Metadata_clone_number", "correlation", "replicate"]]
        .query("replicate == False")
        .reset_index(drop=True)
    )

    x_df = pd.concat([left_df, right_df], axis='rows', ignore_index=True)

    # Now plot
    plot_obj = gg.ggplot(x_df, gg.aes(y="correlation", x="replicate", fill="replicate")) +         gg.geom_jitter(shape = ".", size=0.5, alpha=0.3) +         gg.geom_boxplot(alpha=0.3, outlier_alpha=0) +         gg.scale_fill_manual(name="Replicate",
                             labels={"True": "True",
                                     "False": "False"},
                             values=["#B99638", "#2DB898"]) + \
        gg.xlab("Replicates") + \
        gg.ylab("Pearson Correlation") + \
        gg.theme_light() + \
        gg.theme(subplots_adjust={'wspace': 0.2},
                 axis_text=gg.element_text(size=7),
                 axis_title=gg.element_text(size=9),
                 strip_text=gg.element_text(size=6, color="black"),
                 strip_background=gg.element_rect(colour="black",
                                                  fill="#fdfff4")) + \
        gg.facet_wrap("~Metadata_clone_number")
    return plot_obj


# In[3]:


batch = "2019_06_25_Batch3"
backend_dir = os.path.join("..", "..", "backend", batch)

backend_folders = os.listdir(backend_dir)
mut_file, wt_file = [
    os.path.join(backend_dir, x, "{}_normalized_variable_selected.csv".format(x))
    for x in backend_folders
]

print(mut_file)
print(wt_file)


# ## Read in Files

# In[4]:


mut_df = pd.read_csv(mut_file)

print(mut_df.shape)
mut_df.head(2)


# In[5]:


wt_df = pd.read_csv(wt_file)

print(wt_df.shape)
wt_df.head(2)


# ## How Many Replicates?

# In[6]:


mut_df.Metadata_clone_number.value_counts()


# In[7]:


wt_df.Metadata_clone_number.value_counts()


# ## Get Pairwise Correlations

# In[8]:


mut_cor_df = get_pairwise_corr(mut_df)

print(mut_cor_df.shape)
mut_cor_df.head()


# In[9]:


wt_cor_df = get_pairwise_corr(wt_df)

print(wt_cor_df.shape)
wt_cor_df.head()


# ## Plot Replicate Correlation

# In[10]:


wt_cor_gg = plot_replicate_corr(wt_cor_df)

file = os.path.join("figures", "wt_correlation_{}.png".format(batch))
wt_cor_gg.save(filename=file, height=4, width=5, dpi=300)

wt_cor_gg


# In[11]:


mut_cor_gg = plot_replicate_corr(mut_cor_df)

file = os.path.join("figures", "mut_correlation_{}.png".format(batch))
mut_cor_gg.save(filename=file, height=4, width=5, dpi=300)

mut_cor_gg


# In[12]:


wt_cor_individual_gg = plot_individual_replicate_corr(wt_cor_df)

file = os.path.join("figures", "wt_correlation_individual_{}.png".format(batch))
wt_cor_individual_gg.save(filename=file, height=4, width=5, dpi=300)

wt_cor_individual_gg


# In[13]:


mut_cor_individual_gg = plot_individual_replicate_corr(mut_cor_df)

file = os.path.join("figures", "mut_correlation_individual_{}.png".format(batch))
mut_cor_individual_gg.save(filename=file, height=4, width=5, dpi=300)

mut_cor_individual_gg


# ## Output Median Replicate Correlation

# In[14]:


wt_median_cor_df = (
    wt_cor_df
    .groupby(["Metadata_clone_number", "replicate"])
    .median()
    .reset_index()
)

mut_median_cor_df = (
    mut_cor_df
    .groupby(["Metadata_clone_number", "replicate"])
    .median()
    .reset_index()
)

# Save Files
file = os.path.join("results", "WT_median_replicate_correlation.tsv")
wt_median_cor_df.to_csv(file, sep='\t', index=False)

file = os.path.join("results", "MUT_median_replicate_correlation.tsv")
mut_median_cor_df.to_csv(file, sep='\t', index=False)

