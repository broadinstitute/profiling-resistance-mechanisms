"""
Helper functions to enable profile generation and processing
"""

import os
import yaml
import pandas as pd

from pycytominer.aggregate import AggregateProfiles
from pycytominer import (
    annotate,
    normalize,
    feature_select,
    audit
)


def process_pipeline(pipeline, option):
    if option == "compression":
        if option in pipeline.keys():
            output = pipeline["compression"]
        else:
            output = "None"

    if option == "samples":
        if option in pipeline.keys():
            output = pipeline["samples"]
        else:
            output = "all"

    return output


def process_profile(sql_file, batch, plate, pipeline):
    """
    Given batch details and a pipeline, process morphology profiles
    """
    assert batch in sql_file, "batch {} not recognized in sql file {}".format(
        batch, sql_file
    )
    assert plate in sql_file, "plate {} not recognized in sql file {}".format(
        plate, sql_file
    )

    # Set output directory information
    pipeline_output = pipeline["output_dir"]
    output_dir = os.path.join(pipeline_output, batch, plate)
    os.makedirs(output_dir, exist_ok=True)

    # Set output file information
    aggregate_out_file = os.path.join(output_dir, "{}.csv.gz".format(plate))
    annotate_out_file = os.path.join(output_dir, "{}_augmented.csv.gz".format(plate))
    normalize_out_file = os.path.join(output_dir, "{}_normalized.csv.gz".format(plate))
    feature_out_file = os.path.join(
        output_dir, "{}_normalized_feature_selected.csv.gz".format(plate)
    )

    # Load pipeline options
    compression = process_pipeline(pipeline["options"], option="compression")
    samples = process_pipeline(pipeline["options"], option="samples")

    # Load and setup platemap info
    workspace_dir = pipeline["workspace_dir"]
    batch_dir = os.path.join(workspace_dir, "backend", batch)
    metadata_dir = os.path.join(workspace_dir, "metadata", batch)

    barcode_plate_map_file = os.path.join(
        metadata_dir, sorted(os.listdir(metadata_dir))[0]
    )
    barcode_plate_map_df = pd.read_csv(barcode_plate_map_file)
    plate_map_name = barcode_plate_map_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    plate_map_file = os.path.join(
        metadata_dir, "platemap", "{}.txt".format(plate_map_name)
    )
    plate_map_df = pd.read_csv(plate_map_file, sep="\t")
    plate_map_df.columns = [
        "Metadata_{}".format(x) if not x.startswith("Metadata_") else x
        for x in plate_map_df.columns
    ]
    platemap_well_column = pipeline["platemap_well_column"]

    # Step 1: Aggregate
    aggregate_steps = pipeline["aggregate"]
    if aggregate_steps["perform"]:
        aggregate_features = aggregate_steps["features"]
        aggregate_operation = aggregate_steps["method"]
        aggregate_plate_column = aggregate_steps["plate_column"]
        aggregate_well_column = aggregate_steps["well_column"]

        ap = AggregateProfiles(
            sql_file,
            strata=[aggregate_plate_column, aggregate_well_column],
            features=aggregate_features,
            operation=aggregate_operation,
        )

        ap.aggregate_profiles(output_file=aggregate_out_file, compression=compression)

        if pipeline["count"]["perform"]:
            count_dir = pipeline["count"]["output_dir"]

            cell_count_file = os.path.join(
                count_dir, "{}_{}_cell_count.tsv".format(batch, plate)
            )

            cell_count_df = ap.count_cells()

            cell_count_df = cell_count_df.merge(
                plate_map_df,
                left_on=aggregate_well_column,
                right_on=platemap_well_column,
            ).drop(platemap_well_column, axis="columns")

    # Annotate Profiles
    annotate_steps = pipeline["annotate"]
    if annotate_steps["perform"]:
        annotate_well_column = annotate_steps["well_column"]
        annotate(
            profiles=aggregate_out_file,
            platemap=plate_map_df,
            join_on=[platemap_well_column, annotate_well_column],
            output_file=annotate_out_file,
            compression=compression,
        )

    # Normalize Profiles
    normalize_steps = pipeline["normalize"]
    if normalize_steps["perform"]:
        norm_features = normalize_steps["features"]
        norm_method = normalize_steps["method"]
        normalize(
            profiles=annotate_out_file,
            features=norm_features,
            samples=samples,
            method=norm_method,
            output_file=normalize_out_file,
            compression=compression,
        )

    # Apply feature selection
    feature_select_steps = pipeline["feature_select"]
    if feature_select_steps["perform"]:
        feature_select_operations = feature_select_steps["operations"]
        feature_select_features = feature_select_steps["features"]
        feature_select(
            profiles=normalize_out_file,
            features=feature_select_features,
            samples="none",
            operation=feature_select_operations,
            output_file=feature_out_file,
            compression=compression,
            corr_threshold=0.9,
            corr_method="pearson",
        )
