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
batch = "2019_03_20_Batch2"
plate = "207106_exposure320"

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
pipeline, single_cell_files = load_config(config, append_sql_prefix=False)


# In[5]:


workspace_dir = pipeline["workspace_dir"]
batch_dir = pathlib.Path(workspace_dir, "backend", batch)
metadata_dir = pathlib.Path(workspace_dir, "metadata", batch)

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


image_df.Metadata_CellLine.value_counts()


# In[12]:


image_df.Metadata_Well.value_counts()


# In[13]:


clone_e_wells = pd.np.random.choice(
    (
        image_df
        .query("Metadata_CellLine == 'CloneE'")
        .query("Metadata_Dosage == 0")
    )
    .Metadata_Well.unique(), size=2, replace=False
)

wt_wells = pd.np.random.choice(
    (
        image_df
        .query("Metadata_CellLine == 'WT'")
        .query("Metadata_Dosage == 0")
    ).Metadata_Well.unique(), size=2, replace=False
)

clone_e_holdout_wells = pd.np.random.choice(
    (
        image_df
        .query("Metadata_CellLine == 'CloneE'")
        .query("Metadata_Dosage == 0")
        .query("Metadata_Well not in @clone_e_wells")
    )
    .Metadata_Well.unique(), size=1, replace=False
)

wt_holdout_wells = pd.np.random.choice(
    (
        image_df
        .query("Metadata_CellLine == 'WT'")
        .query("Metadata_Dosage == 0")
        .query("Metadata_Well not in @wt_wells")
    ).Metadata_Well.unique(), size=1, replace=False
)

clone_a_wells = pd.np.random.choice(
    (
        image_df
        .query("Metadata_CellLine == 'CloneA'")
        .query("Metadata_Dosage == 0")
    )
    .Metadata_Well.unique(), size=1, replace=False
)

print(
    f"Clone E Wells: {clone_e_wells}",
    f"\nWT Wells: {wt_wells}",
    f"Clone E Holdout Wells: {clone_e_holdout_wells}",
    f"\nWT Holdout Wells: {wt_holdout_wells}",
    f"\nClone A Wells: {clone_a_wells}"
)


# # Load Cells

# In[14]:


imagenumber_dict = {}
imagenumber_dict["clone_e"] = image_df.query("Metadata_Well in @clone_e_wells").ImageNumber.tolist()
imagenumber_dict["wt"] = image_df.query("Metadata_Well in @wt_wells").ImageNumber.tolist()
imagenumber_dict["clone_a"] = image_df.query("Metadata_Well in @clone_a_wells").ImageNumber.tolist()
imagenumber_dict["clone_e_holdout"] = image_df.query("Metadata_Well in @clone_e_holdout_wells").ImageNumber.tolist()
imagenumber_dict["wt_holdout"] = image_df.query("Metadata_Well in @wt_holdout_wells").ImageNumber.tolist()

imagenumber_dict


# In[15]:


train_dict_df = {}
test_dict_df = {}
holdout_dict_df = {}
for clone_type, clone_imagenumbers in imagenumber_dict.items():
    print(f"Now processing clone: {clone_type}")
    train_df, test_df = process_sites(
        connection=conn,
        imagenumbers=clone_imagenumbers,
        image_df=image_df,
        feature_filter=feature_filter,
        scaler_method=scaler_method,
        seed=seed,
        test_split_prop=test_split_prop,
        normalize=False 
    )
    print(train_df.shape)
    print(test_df.shape)
    
    if clone_type in ["clone_e", "wt"]:
        train_dict_df[clone_type] = train_df.reset_index(drop=True)
        test_dict_df[clone_type] = test_df.reset_index(drop=True)
    else:
        holdout_dict_df[clone_type] = train_df


# In[16]:


# Normalize and shuffle row order
train_df = normalize_sc(pd.concat(train_dict_df).reset_index(drop=True), scaler_method=scaler_method)
test_df = normalize_sc(pd.concat(test_dict_df).reset_index(drop=True), scaler_method=scaler_method)

train_df = train_df.sample(frac=1).reset_index(drop=True)
test_df = test_df.sample(frac=1).reset_index(drop=True)


# ## Apply Feature Selection

# In[17]:


# Original shapes
print(train_df.shape)
print(test_df.shape)


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


# In[20]:


# Shapes after feature selection
print(train_df.shape)
print(test_df.shape)


# ## Output Files

# In[21]:


out_file = pathlib.Path("data", "example_train.tsv.gz")
train_df.to_csv(out_file, sep="\t", compression="gzip", index=False)

out_file = pathlib.Path("data", "example_test.tsv.gz")
test_df.to_csv(out_file, sep="\t", compression="gzip", index=False)


# In[22]:


# Normalize, shuffle row order, and output for holdout sets
for clone_type in holdout_dict_df:
    df = normalize_sc(holdout_dict_df[clone_type].reset_index(drop=True), scaler_method=scaler_method)
    df = df.sample(frac=1).reset_index(drop=True).reindex(reindex_features, axis="columns")
    print(clone_type)
    print(df.shape)
    out_file = pathlib.Path("data", f"example_holdout_{clone_type}.tsv.gz")
    df.to_csv(out_file, sep="\t", compression="gzip", index=False)

