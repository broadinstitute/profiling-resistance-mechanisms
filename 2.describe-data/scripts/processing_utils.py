"""
Various processing utility functions

Usage:
import only
"""

import os
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


def get_recode_cols():
    return_dict = {}
    return_dict["recode_cols"] = {
        "Metadata_CellLine": "Metadata_clone_number",
        "Metadata_Dosage": "Metadata_treatment",
    }

    return_dict["recode_sample"] = {
        "Clone A": "CloneA",
        "Clone E": "CloneE",
        "WT": "WT_parental",
        "WT parental": "WT_parental",
    }

    return_dict["recode_treatment"] = {
        "0.0": "0.1% DMSO",
        "DMSO": "0.1% DMSO",
        "0.7": "2.1 nM bortezomib",
        "7.0": "21 nM bortezomib",
        "70.0": "210 nM bortezomib",
        "bortezomib": "21 nM bortezomib",
    }

    return return_dict


def load_data(
    batch,
    plates="all",
    profile_dir="profiles",
    suffix="normalized_feature_selected.csv.gz",
    combine_dfs=False,
    add_cell_count=False,
    harmonize_cols=False,
    cell_count_dir="cell_counts",
):
    batch_dir = os.path.join(profile_dir, batch)

    plate_folders = [x for x in os.listdir(batch_dir) if ".DS_Store" not in x]

    plate_files = [
        os.path.join(batch_dir, x, f"{x}_{suffix}")
        for x in plate_folders
        if ".DS_Store" not in x
    ]

    plate_data = {}
    for plate_idx in range(0, len(plate_files)):
        plate = plate_folders[plate_idx]
        if plates != "all":
            if plate not in plates:
                continue

        df = pd.read_csv(plate_files[plate_idx]).assign(Metadata_batch=batch)

        if add_cell_count:
            df = merge_cell_count(df, batch, cell_count_dir=cell_count_dir)

        if harmonize_cols:
            recode_cols = get_recode_cols()

            # Update columns and specific entries
            if batch == "2019_06_25_Batch3":
                df = df.assign(Metadata_treatment="Untreated")

            df = df.rename(recode_cols["recode_cols"], axis="columns")

            df.Metadata_clone_number = df.Metadata_clone_number.astype(str)
            df.Metadata_treatment = df.Metadata_treatment.astype(str)

            df.Metadata_clone_number = df.Metadata_clone_number.replace(
                recode_cols["recode_sample"]
            )
            df.Metadata_treatment = df.Metadata_treatment.replace(
                recode_cols["recode_treatment"]
            )

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
