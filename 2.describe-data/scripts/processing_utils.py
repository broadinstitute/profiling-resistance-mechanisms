"""
Various processing utility functions

Usage:
import only
"""

import os
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


def load_data(
    batch,
    profile_dir="profiles",
    suffix="normalized_feature_selected.csv.gz",
    combine_dfs=False,
):
    batch_dir = os.path.join(profile_dir, batch)

    plate_folders = os.listdir(batch_dir)

    plate_files = [
        os.path.join(batch_dir, x, "{}_{}".format(x, suffix))
        for x in plate_folders
        if ".DS_Store" not in x
    ]

    plate_data = {}
    for plate_idx in range(0, len(plate_files)):
        plate = plate_folders[plate_idx]
        plate_data[plate] = pd.read_csv(plate_files[plate_idx]).assign(
            Metadata_batch=batch
        )

    if combine_dfs:
        plate_data = convert_data(plate_data)

    return plate_data


def convert_data(df_dict):
    df = pd.concat(df_dict.values(), ignore_index=True, sort=True).reset_index(
        drop=True
    )
    cp_cols = infer_cp_features(df)
    meta_cols = df.drop(cp_cols, axis="columns").columns.tolist()

    return df.reindex(meta_cols + cp_cols, axis="columns")
