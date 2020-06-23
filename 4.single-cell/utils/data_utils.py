import pathlib
import pandas as pd
from pycytominer.cyto_utils import infer_cp_features


def load_data(
    y_col="Metadata_CellLine", wt_col="WT", return_meta=False, shuffle_row_order=False
):
    train_file = pathlib.Path("data", "example_train.tsv.gz")
    train_df = pd.read_csv(train_file, sep="\t")

    test_file = pathlib.Path("data", "example_test.tsv.gz")
    test_df = pd.read_csv(test_file, sep="\t")

    if shuffle_row_order:
        train_df = train_df.sample(frac=1).reset_index(drop=True)
        test_df = test_df.sample(frac=1).reset_index(drop=True)

    y_train_df = pd.DataFrame(train_df.loc[:, y_col]).assign(status=0)
    y_train_df.loc[y_train_df.loc[:, y_col] != wt_col, "status"] = 1

    y_test_df = pd.DataFrame(test_df.loc[:, y_col]).assign(status=0)
    y_test_df.loc[y_test_df.loc[:, y_col] != wt_col, "status"] = 1

    cp_features = infer_cp_features(train_df)
    x_train_df = train_df.loc[:, cp_features]
    x_test_df = test_df.loc[:, cp_features]

    if return_meta:
        meta_train_df = train_df.drop(cp_features, axis="columns")
        meta_test_df = test_df.drop(cp_features, axis="columns")
        return x_train_df, y_train_df, meta_train_df, x_test_df, y_test_df, meta_test_df
    else:
        return x_train_df, y_train_df, x_test_df, y_test_df
