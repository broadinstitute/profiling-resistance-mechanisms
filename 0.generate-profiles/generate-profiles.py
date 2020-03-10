"""
Script to generate profiles given a configuration file
"""

import os
import yaml
import argparse
import numpy as np
import pandas as pd

from scripts.profile_util import process_profile

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config", help="configuration yaml file for pipeline and batch information"
)
args = parser.parse_args()
config = args.config

# Load configuration file info
profile_config = {}
with open(config, "r") as stream:
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
                x: "sqlite:////{}".format(
                    os.path.join(
                        pipeline["workspace_dir"],
                        "backend",
                        batch,
                        x,
                        "{}.sqlite".format(x),
                    )
                )
                for x in plates
            }


for batch in profile_config:
    for plate in profile_config[batch]["plates"]:
        sql_file = profile_config[batch]["plates"][plate]
        print("Now processing... batch: {}, plate: {}".format(batch, plate))
        process_profile(sql_file=sql_file, batch=batch, plate=plate, pipeline=pipeline)
