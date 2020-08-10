import pathlib
import pandas as pd
from pycytominer.cyto_utils import infer_cp_features


def load_data(return_meta=False, shuffle_row_order=False, holdout=False, othertreatment=False):
    output_data_dict = {"train": {}, "test": {}}
    train_file = pathlib.Path("data", "single_cell_train.tsv.gz")
    train_df = pd.read_csv(train_file, sep="\t")

    test_file = pathlib.Path("data", "single_cell_test.tsv.gz")
    test_df = pd.read_csv(test_file, sep="\t")

    if shuffle_row_order:
        train_df = train_df.sample(frac=1).reset_index(drop=True)
        test_df = test_df.sample(frac=1).reset_index(drop=True)

    cp_features = infer_cp_features(train_df)
    output_data_dict["train"]["x"] = train_df.loc[:, cp_features]
    output_data_dict["test"]["x"] = test_df.loc[:, cp_features]

    if holdout:
        output_data_dict["holdout"] = {}
        holdout_file = pathlib.Path("data", "single_cell_holdout.tsv.gz")
        holdout_df = pd.read_csv(holdout_file, sep="\t")
        if shuffle_row_order:
            holdout_df = holdout_df.sample(frac=1).reset_index(drop=True)

        output_data_dict["holdout"]["x"] = holdout_df.loc[:, cp_features]

    if othertreatment:
        output_data_dict["othertreatment"] = {}
        other_file = pathlib.Path("data", "single_cell_othertreatment.tsv.gz")
        other_df = pd.read_csv(other_file, sep="\t")
        if shuffle_row_order:
            other_df = other_df.sample(frac=1).reset_index(drop=True)

        output_data_dict["othertreatment"]["x"] = other_df.loc[:, cp_features]

    if return_meta:
        output_data_dict["train"]["meta"] = train_df.drop(cp_features, axis="columns")
        output_data_dict["test"]["meta"] = test_df.drop(cp_features, axis="columns")

        if holdout:
            output_data_dict["holdout"]["meta"] = holdout_df.drop(cp_features, axis="columns")

        if othertreatment:
            output_data_dict["othertreatment"]["meta"] = other_df.drop(cp_features, axis="columns")

    return output_data_dict
