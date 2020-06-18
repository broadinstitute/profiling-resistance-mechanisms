"""
Script to generate profiles given a configuration file
"""

import os
import yaml
import argparse
import numpy as np
import pandas as pd

from scripts.profile_util import process_profile, load_config

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config", help="configuration yaml file for pipeline and batch information"
)
args = parser.parse_args()
config = args.config

# Load configuration file info
pipeline, profile_config = load_config(config)

for batch in profile_config:
    for plate in profile_config[batch]["plates"]:
        sql_file = profile_config[batch]["plates"][plate]
        print("Now processing... batch: {}, plate: {}".format(batch, plate))
        process_profile(sql_file=sql_file, batch=batch, plate=plate, pipeline=pipeline)
