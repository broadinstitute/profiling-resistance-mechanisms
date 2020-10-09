import os
import yaml
import argparse
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates

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
output_file_extensions = [".png"]

audit_config = {}
stream = open(config, "r")
for data in yaml.load_all(stream, Loader=yaml.FullLoader):
    batch = data["batch"]
    audit_level = data["auditlevel"]
    plates = [str(x) for x in data["plates"]]
    audit_config[batch] = {}
    audit_config[batch]["plates"] = plates
    audit_config[batch]["auditcols"] = data["auditcols"]
    audit_config[batch]["process"] = data["process"]
    audit_config[batch]["plate_files"] = {
        x: os.path.join(profile_dir, batch, x, "{}_{}.csv.gz".format(x, audit_level))
        for x in plates
    }

for batch in audit_config:
    batch_dict = audit_config[batch]
    process = batch_dict["process"]
    if not process:
        continue
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
        
        # Determine feature class
        features = infer_cp_features(df)
        meta_features = infer_cp_features(df, metadata=True)
        
        # Calculate and process pairwise similarity matrix
        audit_df = metric_melt(
            df=df,
            features=features,
            metadata_features=meta_features,
            similarity_metric="pearson",
            eval_metric="percent_strong"
        )

        audit_df = assign_replicates(
            similarity_melted_df=audit_df, replicate_groups=audit_cols
        )
        # What is 95% of the non replicate null distribution
        cutoff = audit_df.query("not group_replicate").similarity_metric.quantile(0.95)

        # Calculate a single number for percent strong
        percent_strong = evaluate(
            profiles=df,
            features=features,
            meta_features=meta_features,
            replicate_groups=audit_cols,
            operation="percent_strong",
            similarity_metric="pearson",
            percent_strong_quantile=0.95
        )
        
        grid_string = "~{}".format("+".join([f"{x}_pair_a" for x in audit_cols]))
        
        # Visualize the audit - output two plots for each plate
        output_base = os.path.join(
            figure_output_dir, "{}_{}_replicate_correlation".format(batch, plate)
        )
        plot_replicate_correlation(
            df=audit_df,
            batch=batch,
            plate=plate,
            facet_string=grid_string,
            dpi=500,
            split_samples=True,
            output_file_base=output_base,
            output_file_extensions=output_file_extensions,
        )

        output_base = os.path.join(
            figure_output_dir, "{}_{}_density".format(batch, plate)
        )
        plot_replicate_density(
            df=audit_df,
            batch=batch,
            plate=plate,
            cutoff=cutoff,
            percent_strong=percent_strong,
            dpi=500,
            output_file_base=output_base,
            output_file_extensions=output_file_extensions,
        )