import os
import yaml
import argparse
import numpy as np
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import drop_outlier_features
from pycytominer import audit

parser = argparse.ArgumentParser()
parser.add_argument("--config", help="configuration yaml file for batch information")
parser.add_argument(
    "--profile_dir",
    help="directory storing profiles",
    default="1.process-profiles/profiles",
)
parser.add_argument(
    "--output_dir",
    help="directory where to save audit results",
    default="2.profiling-audit/results",
)
args = parser.parse_args()

config = args.config
output_dir = args.output_dir
profile_dir = args.profile_dir

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
        x: os.path.join(profile_dir, batch, x, "{}_{}.csv".format(x, audit_level))
        for x in plates
    }

for batch in audit_config:
    batch_dict = audit_config[batch]
    audit_cols = batch_dict["auditcols"]
    plate_files = batch_dict["plate_files"]
    plates = batch_dict["plates"]
    for plate in plates:
        audit_output_dir = os.path.join(output_dir, batch, plate)
        os.makedirs(audit_output_dir, exist_ok=True)

        audit_output_file = os.path.join(audit_output_dir, "{}_audit.csv".format(plate))
        df = pd.read_csv(plate_files[plate])
        drop_columns = drop_outlier_features(population_df=df)

        df = df.drop(drop_columns, axis="columns")

        audit(
            df,
            audit_groups=audit_cols,
            audit_resolution="full",
            output_file=audit_output_file,
        )
