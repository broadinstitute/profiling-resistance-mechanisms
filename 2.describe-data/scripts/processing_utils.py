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
    add_cell_count=False,
    cell_count_dir="cell_counts",
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
        df = pd.read_csv(plate_files[plate_idx]).assign(Metadata_batch=batch)
        if add_cell_count:
            df = merge_cell_count(df, batch, cell_count_dir=cell_count_dir)

        plate_data[plate] = df

    if combine_dfs:
        plate_data = convert_data(plate_data)

    return plate_data


def merge_cell_count(df, batch, cell_count_dir="cell_counts"):
    # Load cell counts for the specific plates
    count_files = [
        os.path.join(cell_count_dir, x)
        for x in os.listdir(cell_count_dir)
        if batch in x
    ]
    all_plate_dfs = []
    for count_file in count_files:
        plate = os.path.basename(count_file)
        plate = plate.replace(batch, "").replace("cell_count.tsv", "").strip("_")

        plate_df = pd.read_csv(count_file, sep="\t").rename(
            {"cell_count": "Metadata_cell_count"}, axis="columns"
        )
        all_plate_dfs.append(plate_df)

    # Merge all plates and append cell count information as a metadata feature
    plate_df = pd.concat(all_plate_dfs, sort=True)
    df = plate_df.merge(
        df, on=plate_df.drop("Metadata_cell_count", axis="columns").columns.tolist()
    )
    return df


def convert_data(df_dict):
    df = pd.concat(df_dict.values(), ignore_index=True, sort=True).reset_index(
        drop=True
    )
    cp_cols = infer_cp_features(df)
    meta_cols = df.drop(cp_cols, axis="columns").columns.tolist()

    return df.reindex(meta_cols + cp_cols, axis="columns")
