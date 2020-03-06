import os
import yaml
import argparse
import numpy as np
import pandas as pd

from pycytominer import audit

from scripts.viz_utils import plot_replicate_correlation, plot_replicate_density

parser = argparse.ArgumentParser()
parser.add_argument("--config", help="configuration yaml file for batch information")
parser.add_argument(
    "--profile_dir",
    help="directory storing profiles",
    default="../0.generate-profiles/profiles",
)
parser.add_argument(
    "--output_dir", help="directory where to save audit results", default="results"
)
parser.add_argument(
    "--figure_dir", help="directory where to save audit figures", default="figures"
)
args = parser.parse_args()

config = args.config
profile_dir = args.profile_dir
output_dir = args.output_dir
figure_dir = args.figure_dir

np.random.seed(1234)

audit_config = {}
stream = open(config, "r")
for data in yaml.load_all(stream, Loader=yaml.FullLoader):
    batch = data["batch"]
    audit_level = data["auditlevel"]
    plates = [str(x) for x in data["plates"]]
    audit_config[batch] = {}
    audit_config[batch]["plates"] = plates
    audit_config[batch]["auditcols"] = data["auditcols"]
    audit_config[batch]["plate_files"] = {
        x: os.path.join(profile_dir, batch, x, "{}_{}.csv.gz".format(x, audit_level))
        for x in plates
    }

for batch in audit_config:
    batch_dict = audit_config[batch]
    audit_cols = batch_dict["auditcols"]
    plate_files = batch_dict["plate_files"]
    plates = batch_dict["plates"]
    for plate in plates:
        print("Now auditing... Batch: {}; Plate: {}".format(batch, plate))
        audit_output_dir = os.path.join(output_dir, batch, plate)
        os.makedirs(audit_output_dir, exist_ok=True)

        figure_output_dir = os.path.join(figure_dir, batch, plate)
        os.makedirs(figure_output_dir, exist_ok=True)

        audit_output_file = os.path.join(audit_output_dir, "{}_audit.csv".format(plate))
        df = pd.read_csv(plate_files[plate])

        audit(
            df,
            audit_groups=audit_cols,
            audit_resolution="full",
            output_file=audit_output_file,
        )

        audit_df = pd.read_csv(audit_output_file)
        pair_a_df = audit_df.loc[:, ["{}_pair_a".format(x) for x in audit_cols]]
        pair_a_df.columns = audit_cols
        pair_b_df = audit_df.loc[:, ["{}_pair_b".format(x) for x in audit_cols]]
        pair_b_df.columns = audit_cols
        audit_df = audit_df.assign(
            replicate_info=(pair_a_df == pair_b_df).all(axis="columns")
        )

        # Build a dataframe that does not skip any sample combinations
        pair_a_df = pair_a_df.assign(
            replicate_info=audit_df.replicate_info,
            pairwise_correlation=audit_df.pairwise_correlation,
        )
        pair_b_df = pair_b_df.assign(
            replicate_info=audit_df.replicate_info,
            pairwise_correlation=audit_df.pairwise_correlation,
        )
        full_audit_df = pd.concat([pair_a_df, pair_b_df], axis="rows").drop_duplicates()

        grid_string = "~{}".format("+".join(audit_cols))

        # Visualize the audit - output two plots for each plate
        output_base = os.path.join(
            figure_output_dir, "{}_{}_replicate_correlation".format(batch, plate)
        )
        _ = plot_replicate_correlation(
            full_audit_df,
            batch,
            plate,
            grid_string,
            dpi=400,
            split_samples=True,
            output_file_base=output_base,
        )

        output_base = os.path.join(
            figure_output_dir, "{}_{}_density".format(batch, plate)
        )
        _ = plot_replicate_density(
            audit_df, batch, plate, dpi=400, output_file_base=output_base
        )
