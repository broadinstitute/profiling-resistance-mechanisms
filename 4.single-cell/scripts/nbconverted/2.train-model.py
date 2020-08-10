#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from joblib import dump
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, f1_score

from pycytominer.cyto_utils import infer_cp_features

from utils.data_utils import load_data
from utils.ml_utils import (
    shuffle_columns,
    model_apply,
    cross_validation_performance,
    output_coefficients,
)


# In[2]:


np.random.seed(12345)


# ## Set parameters for grid search

# In[3]:


cs = [0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 10]
l1_ratios = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 1]

n_folds = 5


# ## Load data

# In[4]:


data_dict = load_data(
    return_meta=True,
    shuffle_row_order=True,
    holdout=True,
    othertreatment=True
)

print(data_dict["train"]["x"].shape)
print(data_dict["test"]["x"].shape)
print(data_dict["holdout"]["x"].shape)
print(data_dict["othertreatment"]["x"].shape)

data_dict["test"]["x"].head(3)


# ## Set Y

# In[5]:


y_column = "Metadata_clone_number"
y_recode = {"WT parental": 0, "Clone A": 1, "Clone E": 2}
y_recode_reverse = {y: x for x, y in y_recode.items()}

y_train = data_dict["train"]["meta"].loc[:, y_column].replace(y_recode)
y_test = data_dict["test"]["meta"].loc[:, y_column].replace(y_recode)

y_train.head()


# ## Setup Pipeline and Fit

# In[6]:


clf_parameters = {
    'classify__C': cs,
    'classify__l1_ratio': l1_ratios
}

estimator = Pipeline(
    steps=[(
        'classify',
        LogisticRegression(
            multi_class="ovr",
            penalty='elasticnet',
            random_state=0,
            class_weight='balanced',
            solver="saga",
            max_iter=50,
            tol=1e-3
        )
    )]
)

# Custom scorer that optimizes f1 score weighted by class proportion
weighted_f1_scorer = make_scorer(f1_score, average='weighted')

cv_pipeline = GridSearchCV(
    estimator=estimator,
    param_grid=clf_parameters,
    n_jobs=-1,
    cv=n_folds,
    scoring=weighted_f1_scorer,
    return_train_score=True
)
 
shuffle_cv_pipeline = GridSearchCV(
    estimator=estimator,
    param_grid=clf_parameters,
    n_jobs=-1,
    cv=n_folds,
    scoring=weighted_f1_scorer,
    return_train_score=True
)


# In[7]:


get_ipython().run_cell_magic('time', '', 'cv_pipeline.fit(X=data_dict["train"]["x"], y=y_train.values)')


# In[8]:


get_ipython().run_cell_magic('time', '', 'x_train_shuffled_df = data_dict["train"]["x"].apply(shuffle_columns, axis=0, result_type="broadcast")\n\nshuffle_cv_pipeline.fit(X=x_train_shuffled_df, y=y_train.values)')


# ## Visualize Cross Validation Results

# In[9]:


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


# In[10]:


cv_heatmap_file = pathlib.Path("figures", "cross_validation", "single_cell_multiclass_train.png")
cross_validation_performance(cv_pipeline, cv_heatmap_file)


# In[11]:


cv_heatmap_file = pathlib.Path("figures", "cross_validation", "single_cell_multiclass_train_shuffled.png")
cross_validation_performance(shuffle_cv_pipeline, cv_heatmap_file)


# ## Output model coefficients

# In[12]:


coef_file = pathlib.Path("coefficients", "single_cell_multiclass_coefficients.tsv")
coef_df = output_coefficients(
    cv_pipeline, data_dict["train"]["x"], y_recode_reverse, coef_file
)

coef_df.head()


# In[13]:


shuff_coef_file = pathlib.Path("coefficients", "single_cell_multiclass_coefficients_shuffled.tsv")
shuff_coef_df = output_coefficients(
    shuffle_cv_pipeline, data_dict["train"]["x"], y_recode_reverse, shuff_coef_file
)

shuff_coef_df.head()


# ## Apply models to data and save output scores

# In[14]:


scores_dict = {}
for data_fit in data_dict:
    scores_dict[data_fit] = {}
    for shuffled in [True, False]:
        if shuffled:
            trained_pipeline = shuffle_cv_pipeline
        else:
            trained_pipeline = cv_pipeline

        scores_dict[data_fit][f"shuffled_{shuffled}"] = model_apply(
            model=trained_pipeline.best_estimator,
            x_df=data_dict[data_fit]["x"],
            meta_df=data_dict[data_fit]["meta"],
            y_recode=y_recode_reverse,
            data_fit=data_fit,
            shuffled=shuffled
        )


# In[15]:


scores_df = []
for data_fit in data_dict:
    for shuffled in [False, True]:
        scores_df.append(scores_dict[data_fit][f"shuffled_{shuffled}"])

scores_df = pd.concat(scores_df, axis="rows").reset_index(drop=True)

output_file = pathlib.Path("scores", "all_single_cell_scores.tsv.gz")
scores_df.to_csv(output_file, sep="\t", compression="gzip", index=False)

print(scores_df.shape)
scores_df.head()


# ## Output Models to File

# In[16]:


model_file = pathlib.Path("models", "multiclass_cloneAE_wildtype.joblib")
top_model = cv_pipeline.best_estimator_.named_steps["classify"]
dump(top_model, model_file)

shuffle_model_file = pathlib.Path("models", "multiclass_cloneAE_wildtype_shuffled.joblib")
top_shuffle_model = shuffle_cv_pipeline.best_estimator_.named_steps["classify"]
dump(top_shuffle_model, shuffle_model_file)


# In[17]:


top_model

