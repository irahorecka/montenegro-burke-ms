import os

import seaborn as sns

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
sns.set_theme()


def isolate_cols_and_transpose_df(df):
    df = ms.get_df_with_cols_to_keep(df, col_substrings=["Compound Name", "Area"])
    df = ms.transpose_and_reset_idx(df)
    return ms.mv_row_as_header(df, row_idx=0)


def mutate_and_relabel_nutrient_data(df):
    ternary_exp = [
        df["Compound Name"].str.split("_").str[1].str[-1].str.isdigit(),
        df["Compound Name"].str.split("_").str[1].str[-1],
        df["Compound Name"].str.split("_").str[1],
    ]
    df = ms.map_ternary_exp_and_replace_values(df, ternary_exp, new_colname="Sample Group")
    # Drop useless Compound Name column, as we have our Sample Group column
    df = df.drop(columns="Compound Name")

    group_id_to_nutrient = {
        "1": "GLC | ASP",
        "2": "GLC | GLN",
        "3": "GLC | AMN",
        "4": "GAL | ASP",
        "5": "GAL | GLN",
        "6": "GAL | AMN",
    }
    return ms.map_dict_and_replace_values(df, colname="Sample Group", dict_map=group_id_to_nutrient)


def aggregate_mean_std_cv_from_nutrient_data(df):
    df = ms.convert_to_numerics(df)
    # Begin aggregation
    df_mean = ms.group_and_agg(df, colname="Sample Group", agg_type="mean")
    df_std = ms.group_and_agg(df, colname="Sample Group", agg_type="std")
    df_cv = ms.convert_to_numerics(df_std.div(df_mean).reset_index())
    return df_mean, df_std, df_cv


def filter_mean_data_from_control_cv_threshold(df_mean, df_cv, cv_threshold=0.15):
    # Non-numeric column - pop and reintroduce later
    sample_group = df_cv.pop("Sample Group")
    df_cv_adj = df_cv[df_cv < cv_threshold]
    # 5 is the index of our control
    df_cv_adj = df_cv_adj.loc[:, ~df_cv_adj.iloc[5].isna()]
    df_cv_adj["Sample Group"] = sample_group
    df_cv_adj = df_cv_adj.set_index("Sample Group")
    return df_mean[list(df_cv_adj)]


def normalize_nutrient_data_to_control(df):
    # Remove blank row and CTRL
    df = df.drop(index=["Blank", "CTRL"])
    # Divide groups 1-6 by group 3 (i.e., our control) to get relative expression diff
    df.loc[:, "2_2_Dimethylsuccinic_acid":] = df.loc[:, "2_2_Dimethylsuccinic_acid":].div(
        df.iloc[3]["2_2_Dimethylsuccinic_acid":]
    )
    # Drop control row (Glucose and ammonia)
    df = df.drop(index="GLC | AMN")
    # Re-apply float data type
    return df.astype(float)


if __name__ == "__main__":
    # Get pandas df from CSV path
    untargeted_yeast_ms_path = os.path.join(
        DATA_PATH, "exportFile_irahorecka_yeast_nutrient_array_150milliminute_retention_time.csv"
    )
    untargeted_yeast_ms_df = read_csv(untargeted_yeast_ms_path)
    # Manipulate df for aggregation
    untargeted_yeast_ms_df = isolate_cols_and_transpose_df(untargeted_yeast_ms_df)
    untargeted_yeast_ms_df = mutate_and_relabel_nutrient_data(untargeted_yeast_ms_df)

    # Aggregate df and find normalized values in respect to the control
    agg_nutrient_mean, agg_nutrient_std, agg_nutrient_cv = aggregate_mean_std_cv_from_nutrient_data(
        untargeted_yeast_ms_df
    )
    # agg_nutrient_mean = filter_mean_data_from_control_cv_threshold(
    #     agg_nutrient_mean, agg_nutrient_cv, cv_threshold=0.25
    # )
    norm_agg_nutrient_mean = normalize_nutrient_data_to_control(agg_nutrient_mean)
    norm_agg_nuttient_mean_log2 = ms.get_log2_df(
        norm_agg_nutrient_mean, downregulated=False, log2_weight=1
    )

    # Plot results as a hierarchical dendrogram
    down_regulated = sns.clustermap(
        norm_agg_nuttient_mean_log2.fillna(0).T, cmap="rocket_r", z_score=0
    )
    down_regulated.savefig("test_zscore.eps", format="eps")
