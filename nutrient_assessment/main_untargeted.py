import os

import seaborn as sns

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")
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
    df.loc[:, df.columns[0] :] = df.loc[:, df.columns[0] :].div(df.iloc[3][df.columns[0] :])
    # Drop control row (Glucose and ammonia)
    df = df.drop(index="GLC | AMN")
    # Re-apply float data type
    return df.astype(float)


if __name__ == "__main__":
    # Get pandas df from CSV path
    untargeted_yeast_ms_path = os.path.join(
        DATA_PATH, "exportFile_irahorecka_yeast_nutrient_array_350milliminute_retention_time.csv"
    )
    untargeted_yeast_ms_df = read_csv(untargeted_yeast_ms_path)

    # Try below two expressions to read 350milliminute retention time file.
    # These samples contain _REF or _MET suffix for values in column "Compound Name".
    if untargeted_yeast_ms_df["Compound Name"].str[-3:].str.contains("MET").any():
        untargeted_yeast_ms_df = ms.drop_rows_with_substring_in_col_value(
            untargeted_yeast_ms_df, "Compound Name", "REF"
        )
        untargeted_yeast_ms_df["Compound Name"] = untargeted_yeast_ms_df["Compound Name"].str[:-4]

    # Manipulate df for aggregation
    untargeted_yeast_ms_df = isolate_cols_and_transpose_df(untargeted_yeast_ms_df)
    untargeted_yeast_ms_df = mutate_and_relabel_nutrient_data(untargeted_yeast_ms_df)
    # Aggregate df and find normalized values in respect to the control
    agg_nutrient_mean, agg_nutrient_std, agg_nutrient_cv = aggregate_mean_std_cv_from_nutrient_data(
        untargeted_yeast_ms_df
    )
    # Only look at metabolites where the CV % for the control is < 15%
    # agg_nutrient_mean = filter_mean_data_from_control_cv_threshold(
    #     agg_nutrient_mean, agg_nutrient_cv, cv_threshold=0.5
    # )

    # Normalize aggregated mean data to mean of control - perform log2 scaling of results
    norm_agg_nutrient_mean = normalize_nutrient_data_to_control(agg_nutrient_mean)
    norm_agg_nutrient_mean_log2 = ms.get_log2_df_directional(
        norm_agg_nutrient_mean, downregulated=True, log2_weight=1
    )

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(os.path.join(DATA_PATH, "log2_nutrient_mean.csv"))

    # Plot results as a hierarchical dendrogram
    # down_regulated = sns.clustermap(
    #     norm_agg_nutrient_mean_log2.fillna(0).T, cmap="rocket_r", z_score=0
    # )
    # down_regulated.savefig(os.path.join(FIGURES_PATH, "test_zscore.eps"), format="eps")