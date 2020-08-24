#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pathlib
import sqlite3
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features

from utils.single_cell_utils import process_sites, normalize_sc
sys.path.append("../0.generate-profiles")
from scripts.profile_util import load_config


# In[2]:


pd.np.random.seed(1234)


# In[3]:


# Set constants
batch = "2020_07_02_Batch8"
plate = "218360"
cell_line_column = "Metadata_clone_number"
cell_lines = ["Clone A", "Clone E", "WT parental"]

feature_filter = ["Object", "Location", "Count", "Parent"]
test_split_prop = 0.15
scaler_method = "standard"
seed = 123

feature_select_opts = [
    "variance_threshold",
    "drop_na_columns",
    "blacklist",
    "drop_outliers",
]
corr_threshold = 0.8
na_cutoff = 0


# In[4]:


# Load locations of single cell files
config = pathlib.Path("../0.generate-profiles/profile_config.yaml")
pipeline, single_cell_files = load_config(config, append_sql_prefix=False, local=True)


# In[5]:


workspace_dir = pipeline["workspace_dir"]
batch_dir = pathlib.Path(workspace_dir, "backend", batch)
metadata_dir = pathlib.Path("../0.generate-profiles", "metadata", batch)

barcode_plate_map_file = pathlib.Path(metadata_dir, "barcode_platemap.csv")
barcode_plate_map_df = pd.read_csv(barcode_plate_map_file)

barcode_plate_map_df


# In[6]:


plate_map_name = (
    barcode_plate_map_df
    .query("Assay_Plate_Barcode == @plate")
    .Plate_Map_Name
    .values[0]
)

plate_map_file = pathlib.Path(metadata_dir, "platemap", f"{plate_map_name}.txt")
plate_map_df = pd.read_csv(plate_map_file, sep="\t")
plate_map_df.columns = [x if x.startswith("Metadata_") else f"Metadata_{x}" for x in plate_map_df.columns]
plate_map_df.head()


# ## Setup Single Cell Connection

# In[7]:


plate_column = pipeline["aggregate"]["plate_column"]
well_column = pipeline["aggregate"]["well_column"]


# In[8]:


# Establish connection to sqlite file
single_cell_sqlite = single_cell_files[batch]["plates"][plate]
conn = sqlite3.connect(single_cell_sqlite)


# In[9]:


image_cols = f"TableNumber, ImageNumber, {plate_column}, {well_column}"
image_query = f"select {image_cols} from image"
image_df = (
    pd.read_sql_query(image_query, conn)
    .merge(
        plate_map_df,
        left_on=well_column,
        right_on="Metadata_well_position"
    )
    .drop(["Metadata_well_position"], axis="columns")
)

print(image_df.shape)
image_df.head()


# ## Identify Representative Wells

# In[10]:


# Assert that image number is unique
assert len(image_df.ImageNumber.unique()) == image_df.shape[0]


# In[11]:


# How many wells were collected per treatment
replicate_info_df = (
    image_df.loc[:, ["Metadata_Well", cell_line_column, "Metadata_treatment"]]
    .drop_duplicates()
)

pd.crosstab(replicate_info_df.loc[:, cell_line_column], replicate_info_df.Metadata_treatment)


# ## Identify wells to use for training and holdout sets
# 
# There are three wells per replicate cell line and treatment.
# We will select two at random to use in training and use the remaining one as a holdout set.

# In[12]:


untreated_wells = []
imagenumber_dict = {}
for cell_line in cell_lines:
    imagenumber_dict[cell_line] = {}
    wells = (
        image_df
        .query(f"{cell_line_column} == @cell_line")
        .query("Metadata_treatment == '0.1% DMSO'")
    ).Metadata_Well.unique()
    
    train_wells = pd.np.random.choice(wells, size=3, replace=False)
    holdout_wells = [x for x in wells if x not in train_wells]

    untreated_wells.extend(train_wells)
    untreated_wells.extend(holdout_wells)
    
    imagenumber_dict[cell_line]["train"] = (
        image_df
        .query("Metadata_Well in @train_wells")
        .ImageNumber
        .tolist()
    )
    imagenumber_dict[cell_line]["holdout"] = (
        image_df
        .query("Metadata_Well in @holdout_wells")
        .ImageNumber
        .tolist()
    )


# In[13]:


other_wells = [x for x in image_df.Metadata_Well.unique() if x not in untreated_wells]

for cell_line in cell_lines:
    imagenumber_dict[cell_line]["other"] = (
        image_df
        .query("Metadata_clone_number == @cell_line")
        .query("Metadata_Well in @other_wells")
        .ImageNumber
        .tolist()
    )


# ## Load Single Cell Data

# In[14]:


training_dict_df = {}
holdout_dict_df = {}
other_dict_df = {}
for clone_type, clone_info_dict in imagenumber_dict.items():
    for data_split, clone_imagenumbers in clone_info_dict.items():
        print(f"Now loading... {clone_type}, {data_split}")
        sc_df = process_sites(
            connection=conn,
            imagenumbers=clone_imagenumbers,
            image_df=image_df,
            feature_filter=feature_filter,
            seed=seed,
            normalize=False
        )
        if data_split == "holdout":
            holdout_dict_df[clone_type] = sc_df.reset_index(drop=True)
        elif data_split == "train":
            training_dict_df[clone_type] = sc_df.reset_index(drop=True)
        elif data_split == "other":
            other_dict_df[clone_type] = sc_df.reset_index(drop=True)


# ## Normalize, split, and shuffle row order

# In[15]:


# Training and testing sets
train_df = pd.concat(training_dict_df).sample(frac=1).reset_index(drop=True)
train_df = normalize_sc(train_df, scaler_method=scaler_method)

train_df, test_df = train_test_split(
    train_df,
    test_size=test_split_prop,
    stratify=train_df.Metadata_clone_number,
    random_state=seed
)

print(train_df.shape)
print(test_df.shape)


# In[16]:


# Holdout set
holdout_df = pd.concat(holdout_dict_df).sample(frac=1).reset_index(drop=True)
holdout_df = normalize_sc(holdout_df, scaler_method=scaler_method)

print(holdout_df.shape)


# In[17]:


# Other data
other_df = pd.concat(other_dict_df).sample(frac=1).reset_index(drop=True)
other_df = normalize_sc(other_df, scaler_method=scaler_method)

print(other_df.shape)


# ## Apply Feature Selection

# In[18]:


meta_features = infer_cp_features(train_df, metadata=True)
meta_features


# In[19]:


train_df = feature_select(
    train_df,
    operation=feature_select_opts,
    na_cutoff=na_cutoff,
    corr_threshold=corr_threshold
)

selected_features = infer_cp_features(train_df)
reindex_features = meta_features + selected_features

test_df = test_df.reindex(reindex_features, axis="columns")
train_df = train_df.reindex(reindex_features, axis="columns")
holdout_df = holdout_df.reindex(reindex_features, axis="columns")
other_df = other_df.reindex(reindex_features, axis="columns")


# In[20]:


# Shapes after feature selection
print(train_df.shape)
print(test_df.shape)
print(holdout_df.shape)
print(other_df.shape)


# ## Output Files

# In[21]:


out_file = pathlib.Path("data", "single_cell_train.tsv.gz")
train_df.to_csv(out_file, sep="\t", compression="gzip", index=False)

out_file = pathlib.Path("data", "single_cell_test.tsv.gz")
test_df.to_csv(out_file, sep="\t", compression="gzip", index=False)

out_file = pathlib.Path("data", "single_cell_holdout.tsv.gz")
holdout_df.to_csv(out_file, sep="\t", compression="gzip", index=False)

out_file = pathlib.Path("data", "single_cell_othertreatment.tsv.gz")
other_df.to_csv(out_file, sep="\t", compression="gzip", index=False)

