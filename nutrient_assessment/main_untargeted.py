"""
montenegro-burke-ms/nutrient_assessment/main_untargeted.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A standalone script to generate a filtered and normalized
CSV file showing log2 transformed expression differences among
nutrient conditions to the control nutriend condition (GLC | AMN).

The functions written here are hard-coded with values that are
found in 'exportFile_irahorecka_yeast_nutrient_array_350milliminute_retention_time.csv'
I know... not the best practice, but it'll do for now.

The filtering conducted here are:
    - Remove metabolite compounds with less than 3/4 reads in any
    nutrient categories.
    - [OPTIONAL] Don't add metabolite samples that showed less than
    x% (e.g. 15%) CV in the control group (GLC | AMN)
    - [VARYING] Keep only metabolite samples that showed at least
    2-fold increase in expression in any nutrient condition in respect
    to the control (GLC | AMN)

Output data can be found here:
montenegro-burke-ms/nutrient_assessment/data/log2_nutrient_mean.csv
"""

import os

import seaborn as sns

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")
sns.set_theme()


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
    df_count_valid_data = ms.group_and_agg(
        df, colname=agg_colname, agg_type=["count", "size"]
    ).drop(index=["Blank", "CTRL"])

    # Find columns to drop that have less than 3/4 reads (i.e., greater than 1 NA value in any nutrient category)
    cols_to_drop = []
    for colname in list(df_count_valid_data.T):
        cols_to_drop.extend(
            list(ms.get_cols_with_less_than_count_in_row(df_count_valid_data.T, colname, 3))
        )
    # Isolate only unique compounds found to have poor resolution
    cols_to_drop_name = list({group[0] for group in cols_to_drop})
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
    untargeted_yeast_ms_df = isolate_cols_and_transpose_df(
        untargeted_yeast_ms_df, ["Compound Name", "Area"]
    )
    untargeted_yeast_ms_df = mutate_and_relabel_nutrient_data(
        untargeted_yeast_ms_df, src_colname="Compound Name"
    )

    # Begin aggregation
    df = ms.convert_to_numerics(untargeted_yeast_ms_df)
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

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(os.path.join(DATA_PATH, "log2_nutrient_mean.csv"))
