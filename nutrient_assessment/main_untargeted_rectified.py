"""
This file is similar to `main_untargeted`, however the aggregation of data is different.

This file normalizes each biological replicate to the control (GLC | AMN) prior to aggregating the biological replicate.
So we are essentially aggregating the normalized data, whereas in the original file, we are averaging the raw values THEN
aggregating the mean to the mean control value (GLC | AMN).

The differences are not that significant, although with this method, we are able to see more molecules involved in the glutamate
biosynthesis pathway. However, this alone should not be enough reason to begin using this method of aggregation.
"""

import os

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")


def isolate_cols_and_transpose_df(df, col_substrings):
    """Mutates data to keep colname substring and move the metabolite
    compounds as column headers."""
    df = ms.get_df_with_cols_to_keep(df, col_substrings)
    df = ms.transpose_and_reset_idx(df)
    return ms.mv_row_as_header(df, row_idx=0)


def mutate_and_relabel_nutrient_data(df, src_colname, dest_colname="Sample Group"):
    """Renames plate well names to nutrient conditions."""
    ternary_exp = [
        df[src_colname].str.split("_").str[1].str[-1].str.isdigit(),
        df[src_colname].str.split("_").str[1].str[-1],
        df[src_colname].str.split("_").str[1],
    ]
    df = ms.map_ternary_exp_and_replace_values(df, ternary_exp, new_colname=dest_colname)
    # Drop useless source column column, as we have our Sample Group column
    df = df.drop(columns=src_colname)

    group_id_to_nutrient = {
        "1": "GLC | ASP",
        "2": "GLC | GLN",
        "3": "GLC | AMN",
        "4": "GAL | ASP",
        "5": "GAL | GLN",
        "6": "GAL | AMN",
    }
    return ms.map_dict_and_replace_values(df, colname=dest_colname, dict_map=group_id_to_nutrient)


def normalize_metabolite_levels_within_biological_replicate(df, chunk_size):
    """Normalizes metabolite levels to the control (GLC | AMN) within each biological replicate."""
    df = ms.convert_to_numerics(df).set_index("Sample Group")
    df = filter_data_with_more_than_3_reads_among_4_samples(df, "Sample Group").drop(
        index=["Blank", "CTRL"]
    )
    df_concat = ms.get_empty_df_from_df(df)
    list_df = [df[i : i + chunk_size] for i in range(0, df.shape[0], chunk_size)]
    for df_ in list_df:
        df_.loc[:, df_.columns[0] :] = df_.loc[:, df_.columns[0] :].div(
            df_.iloc[3][df_.columns[0] :]
        )
        df_concat = ms.concat_df(df_concat, df_, axis=0)

    return df_concat.reset_index().rename(columns={"index": "Sample Group"})


def aggregate_mean_std_cv_from_nutrient_data(df, agg_colname="Sample Group"):
    """Aggregates dataframe to find mean, std, and cv, grouping by column
    bound to `agg_colname`. Returns all three aggregated dataframes to caller."""
    df = ms.convert_to_numerics(df)
    df = filter_data_with_more_than_3_reads_among_4_samples(df, agg_colname)
    # Begin aggregation
    df_mean = ms.group_and_agg(df, colname=agg_colname, agg_type="mean")
    df_std = ms.group_and_agg(df, colname=agg_colname, agg_type="std")
    df_cv = ms.convert_to_numerics(df_std.div(df_mean).reset_index())
    return df_mean, df_std, df_cv


def filter_data_with_more_than_3_reads_among_4_samples(df, agg_colname):
    """Filters and drops metabolite samples that have more than 1 NaN value in any
    of the nutrient conditions. I.e., there must be greater than 3/4 valid samples in
    all nutrient conditions when assessing a particular metabolite."""
    # Generate aggregated data to count NaN
    df_count_valid_data = ms.group_and_agg(df, colname=agg_colname, agg_type="count")
    try:
        df_count_valid_data = df_count_valid_data.drop(index=["Blank", "CTRL"])
    except:
        pass

    # Find columns to drop that have less than 3/4 reads (i.e., greater than 1 NA value in any nutrient category)
    cols_to_drop = []
    for colname in list(df_count_valid_data.T):
        cols_to_drop.extend(
            list(ms.get_cols_with_less_than_count_in_row(df_count_valid_data.T, colname, 3))
        )
    # Isolate only unique compounds found to have poor resolution
    cols_to_drop_name = list(set(cols_to_drop))
    return df.drop(cols_to_drop_name, axis=1)


def filter_mean_data_from_control_cv_threshold(df_mean, df_cv, cv_threshold=0.15):
    """Returns metabolite samples that have a control CV value less than the cv_threshold.
    I.e., 'GLC | AMN' must have a CV% less than cv_threshold (e.g. 15)."""
    # Non-numeric column - pop and reintroduce later
    sample_group = df_cv.pop("Sample Group")
    df_cv_adj = df_cv[df_cv < cv_threshold]
    # 5 is the index of our control
    df_cv_adj = df_cv_adj.loc[:, ~df_cv_adj.iloc[5].isna()]
    df_cv_adj["Sample Group"] = sample_group
    df_cv_adj = df_cv_adj.set_index("Sample Group")
    return df_mean[list(df_cv_adj)]


def normalize_nutrient_data_to_control(df):
    """Normalizes nutrient data value to the control (GLC | AMN). Drops the control row
    after normalization and returns dataframe to caller."""
    # Remove blank row and CTRL
    try:
        df = df.drop(index=["Blank", "CTRL"])
    except:
        pass
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
    untargeted_yeast_ms_df = isolate_cols_and_transpose_df(
        untargeted_yeast_ms_df, ["Compound Name", "Area"]
    )
    untargeted_yeast_ms_df = mutate_and_relabel_nutrient_data(
        untargeted_yeast_ms_df, src_colname="Compound Name"
    )

    untargeted_yeast_ms_df = normalize_metabolite_levels_within_biological_replicate(
        untargeted_yeast_ms_df, 6
    )

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
        norm_agg_nutrient_mean, downregulated=False, log2_weight=1
    )
    norm_agg_nutrient_mean_log2 = ms.get_df_values_within_range(norm_agg_nutrient_mean_log2, -5, 5)

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(os.path.join(DATA_PATH, "log2_nutrient_mean.csv"))
