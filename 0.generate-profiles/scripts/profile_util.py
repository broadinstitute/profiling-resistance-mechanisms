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
)
from pycytominer.cyto_utils import output


def load_config(config_file, append_sql_prefix=True):
    # Load configuration file info
    if append_sql_prefix:
        sql_prefix = "sqlite:////"
    else:
        sql_prefix = "/"
    profile_config = {}
    with open(config_file, "r") as stream:
        for data in yaml.load_all(stream, Loader=yaml.FullLoader):
            if "pipeline" in data.keys():
                pipeline = data
            else:
                process = data["process"]
                if not process:
                    continue
                batch = data["batch"]
                plates = [str(x) for x in data["plates"]]
                profile_config[batch] = {}
                profile_config[batch]["plates"] = {
                    x: "{}{}".format(
                        sql_prefix,
                        os.path.join(
                            pipeline["workspace_dir"],
                            "backend",
                            batch,
                            x,
                            "{}.sqlite".format(x),
                        ),
                    )
                    for x in plates
                }

    return pipeline, profile_config


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
    
    # Process Bulk profiles
    # Step 1: Aggregate
    aggregate_steps = pipeline["aggregate"]
    if aggregate_steps["perform"]:
        aggregate_features = aggregate_steps["features"]
        aggregate_operation = aggregate_steps["method"]
        aggregate_plate_column = aggregate_steps["plate_column"]
        aggregate_well_column = aggregate_steps["well_column"]

        strata = [aggregate_plate_column, aggregate_well_column]

        if "site_column" in aggregate_steps:
            aggregate_site_column = aggregate_steps["site_column"]
            strata += [aggregate_site_column]

        ap = AggregateProfiles(
            sql_file,
            strata=strata,
            features=aggregate_features,
            operation=aggregate_operation,
        )

        ap.aggregate_profiles(output_file=aggregate_out_file, compression=compression)

        if pipeline["count"]["perform"]:
            count_dir = pipeline["count"]["output_dir"]
            os.makedirs(count_dir, exist_ok=True)

            cell_count_file = os.path.join(
                count_dir, "{}_{}_cell_count.tsv".format(batch, plate)
            )

            cell_count_df = ap.count_cells()

            cell_count_df = cell_count_df.merge(
                plate_map_df,
                left_on=aggregate_well_column,
                right_on=platemap_well_column,
            ).drop(platemap_well_column, axis="columns")

            cell_count_df.to_csv(cell_count_file, sep="\t", index=False)

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
    norm_features = normalize_steps["features"]
    norm_method = normalize_steps["method"]
    if normalize_steps["perform"]:
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
    feature_select_operations = feature_select_steps["operations"]
    feature_select_features = feature_select_steps["features"]
    if feature_select_steps["perform"]:
        feature_select(
            profiles=normalize_out_file,
            features=feature_select_features,
            samples=samples,
            operation=feature_select_operations,
            output_file=feature_out_file,
            compression=compression,
            corr_threshold=0.9,
            corr_method="pearson",
        )
        
    sc_steps = pipeline["single_cell"]
    if sc_steps["perform"]:
        if not aggregate_steps["perform"]:
            ap = AggregateProfiles(
                sql_file,
                strata=strata,
                features=aggregate_features,
                operation=aggregate_operation,
            )
        
        # Load cells
        query = "select * from cells"
        cell_df = pd.read_sql(sql=query, con=ap.conn)
        
        # Load cytoplasm
        query = "select * from cytoplasm"
        cytoplasm_df = pd.read_sql(sql=query, con=ap.conn)
        
        # Load nuclei
        query = "select * from nuclei"
        nuclei_df = pd.read_sql(sql=query, con=ap.conn)
        
        # Merge single cells together
        sc_merged_df = cell_df.merge(
            cytoplasm_df.drop("ObjectNumber", axis="columns"),
            left_on=["TableNumber", "ImageNumber", "ObjectNumber"],
            right_on=["TableNumber", "ImageNumber", "Cytoplasm_Parent_Cells"],
            how="inner"
        ).drop("ObjectNumber", axis="columns").merge(
            nuclei_df,
            left_on=["TableNumber", "ImageNumber", "Cytoplasm_Parent_Nuclei"],
            right_on=["TableNumber", "ImageNumber", "ObjectNumber"],
            how="inner"
        )
        
        # Merge image data info
        sc_merged_df = ap.image_df.merge(
            sc_merged_df, how="right", on=ap.merge_cols
        )
        
        # Make sure column names are correctly prefixed
        prefix = ["Metadata", "Cells", "Cytoplasm", "Nuclei"]
        cols = []
        for col in sc_merged_df.columns:
            if any([col.startswith(x) for x in prefix]):
                cols.append(col)
            else:
                cols.append(f"Metadata_{col}")

        sc_merged_df.columns = cols
        
        sc_merged_df = annotate(
            profiles=sc_merged_df,
            platemap=plate_map_df,
            join_on=[platemap_well_column, annotate_well_column],
            output_file="none"
        )
        
        if sc_steps["normalize"]:
            sc_merged_df = normalize(
                profiles=sc_merged_df,
                features=norm_features,
                samples=samples,
                method=norm_method,
                output_file="none"
            )
            
        if sc_steps["feature_select"]:
            sc_merged_df = feature_select(
                profiles=sc_merged_df,
                features=feature_select_features,
                samples=samples,
                operation=feature_select_operations,
                output_file="none",
                corr_threshold=0.9,
                corr_method="pearson",
            )
           
        sc_pipeline_output = pipeline["sc_output_dir"]
        sc_output_dir = os.path.join(sc_pipeline_output, batch, plate)
        os.makedirs(sc_output_dir, exist_ok=True)

        # Set output file information
        sc_out_file = os.path.join(sc_output_dir, "{}_single_cell.csv.gz".format(plate))
        output(
            df=sc_merged_df,
            output_filename=sc_out_file,
            compression="gzip",
            float_format=float_format
        )
