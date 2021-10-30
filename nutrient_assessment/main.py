import os

from dataproc import read_csv
import dataproc.untargeted_ms as ms


BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")

if __name__ == "__main__":
    untargeted_yeast_ms_path = os.path.join(
        DATA_PATH, "exportFile_irahorecka_yeast_nutrient_array_150milliminute_retention_time.csv"
    )
    untargeted_yeast_ms_df = read_csv(untargeted_yeast_ms_path)
    untargeted_yeast_ms_df = ms.get_df_with_cols_to_keep(
        untargeted_yeast_ms_df, col_substrings=["Compound Name", "Area"]
    )
    untargeted_yeast_ms_df = ms.transpose_and_reset_idx(untargeted_yeast_ms_df)
    untargeted_yeast_ms_df = ms.mv_row_as_header(untargeted_yeast_ms_df, row_idx=0)

    ternary_exp = [
        untargeted_yeast_ms_df["Compound Name"].str.split("_").str[1].str[-1].str.isdigit(),
        untargeted_yeast_ms_df["Compound Name"].str.split("_").str[1].str[-1],
        untargeted_yeast_ms_df["Compound Name"].str.split("_").str[1],
    ]
    untargeted_yeast_ms_df = ms.map_ternary_exp_and_replace_values(
        untargeted_yeast_ms_df, ternary_exp, new_colname="Sample Group"
    )
    # Drop useless Compound Name column, as we have our Sample Group column
    untargeted_yeast_ms_df = untargeted_yeast_ms_df.drop(columns="Compound Name")

    group_id_to_nutrient = {
        "1": "GLC | ASP",
        "2": "GLC | GLN",
        "3": "GLC | AMN",
        "4": "GAL | ASP",
        "5": "GAL | GLN",
        "6": "GAL | AMN",
    }
    untargeted_yeast_ms_df = ms.map_dict_and_replace_values(
        untargeted_yeast_ms_df, colname="Sample Group", dict_map=group_id_to_nutrient
    )
    untargeted_yeast_ms_df = ms.convert_to_numerics(untargeted_yeast_ms_df)

    # Begin aggregation
    untargeted_yeast_ms_df_mean = ms.group_and_agg(
        untargeted_yeast_ms_df, colname="Sample Group", agg_type="mean"
    )
    untargeted_yeast_ms_df_std = ms.group_and_agg(
        untargeted_yeast_ms_df, colname="Sample Group", agg_type="std"
    )
    untargeted_yeast_ms_df_cv = untargeted_yeast_ms_df_std.div(
        untargeted_yeast_ms_df_mean
    ).reset_index()

    # Remove blank row and CTRL
    untargeted_yeast_ms_df_mean = untargeted_yeast_ms_df_mean.drop(index=["Blank", "CTRL"])
    # Divide groups 1-6 by group 3 (i.e. our control) to get relative expression diff
    untargeted_yeast_ms_df_mean.loc[
        :, "2_2_Dimethylsuccinic_acid":
    ] = untargeted_yeast_ms_df_mean.loc[:, "2_2_Dimethylsuccinic_acid":].div(
        untargeted_yeast_ms_df_mean.iloc[3]["2_2_Dimethylsuccinic_acid":]
    )
    # Re-apply float data type
    untargeted_yeast_ms_df_mean = untargeted_yeast_ms_df_mean.drop(index="GLC | AMN")
    untargeted_yeast_ms_df_mean = untargeted_yeast_ms_df_mean.astype(float)
    # Drop control row (Glucose and ammonia)

    untargeted_yeast_ms_df_mean_log2 = ms.get_log2_df(
        untargeted_yeast_ms_df_mean, downregulated=False, log2_weight=1
    )


    import seaborn as sns
    sns.set_theme()

    down_regulated = sns.clustermap(
        untargeted_yeast_ms_df_mean_log2.fillna(0).T, cmap="rocket_r", z_score=0
    )
    down_regulated.savefig("test_zscore.eps", format="eps")
