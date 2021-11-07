import os

import seaborn as sns

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")
sns.set_theme()


def isolate_cols_and_transpose_df(df):
    """Mutates data to keep colname substring and move the metabolite
    compounds as column headers."""
    df = ms.get_df_with_cols_to_keep(df, col_substrings=["DetectedMass", "Area"])
    df = ms.transpose_and_reset_idx(df)
    return ms.mv_row_as_header(df, row_idx=0)


def mutate_and_relabel_nutrient_data(df):
    """Renames plate well names to nutrient conditions."""
    ternary_exp = [
        df["DetectedMass"].str.split("_").str[1].str[-1].str.isdigit(),
        df["DetectedMass"].str.split("_").str[1].str[-1],
        df["DetectedMass"].str.split("_").str[1],
    ]
    df = ms.map_ternary_exp_and_replace_values(df, ternary_exp, new_colname="Sample Group")
    # Drop useless DetectedMass column, as we have our Sample Group column
    df = df.drop(columns="DetectedMass")

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
    """Aggregates dataframe to find mean, std, and cv, grouping by column
    'Sample Group'. Returns all three aggregated dataframes to caller."""
    df = ms.convert_to_numerics(df)
    df = filter_data_with_more_than_3_reads_among_4_samples(df, "Sample Group")
    # Begin aggregation
    df_mean = ms.group_and_agg(df, colname="Sample Group", agg_type="mean")
    df_std = ms.group_and_agg(df, colname="Sample Group", agg_type="std")
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
    small_molecule_path = os.path.join(
        DATA_PATH,
        "exportFile_irahorecka_yeast_nutrient_array_batch_recursive_small_molecule_350milliminute_retention_time.csv",
    )
    small_molecule_df = read_csv(small_molecule_path).rename(columns={"Mass": "DetectedMass"})
    small_molecule_df["DetectedMass"] = (
        small_molecule_df["DetectedMass"].map(str) + "_" + small_molecule_df["RT"].map(str)
    )

    small_molecule_df.drop_duplicates(subset="DetectedMass", keep="first", inplace=True)
    small_molecule_df.drop(columns=["Compound Name", "Formula", "CAS ID"], inplace=True)
    small_molecule_df.dropna(inplace=True)

    small_molecule_df = isolate_cols_and_transpose_df(small_molecule_df)
    small_molecule_df = mutate_and_relabel_nutrient_data(small_molecule_df)
    agg_nutrient_mean, agg_nutrient_std, agg_nutrient_cv = aggregate_mean_std_cv_from_nutrient_data(
        small_molecule_df
    )

    # Normalize aggregated mean data to mean of control - perform log2 scaling of results
    norm_agg_nutrient_mean = normalize_nutrient_data_to_control(agg_nutrient_mean)
    norm_agg_nutrient_mean_log2 = ms.get_log2_df(norm_agg_nutrient_mean, log2_weight=1)

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(
        os.path.join(DATA_PATH, "log2_nutrient_mean_small_molecule.csv")
    )

    # SNIPPET TO ASSESS SPEARMAN CORRELATION BETWEEN ABS VALUES OF LOG2 MEAN VS CV
    # I.E., DOES A LARGER EXPRESSION PATTERN CORRELATE CLOSELY WITH LARGER CV?
    # OFF THE BAT, WE CAN SEE A SPEARMAN CORRELATION VALUE OF 0.47 - MODERATE
    # -------------------------------------
    # mean_list = norm_agg_nutrient_mean_log2.abs().values.tolist()
    # cv_list = agg_nutrient_cv.set_index("Sample Group").drop(index=['Blank', 'CTRL', 'GLC | AMN'])[list(norm_agg_nutrient_mean_log2)].values.tolist()
    # mean_mod = [val for l in mean_list for val in l]
    # cv_mod = [val for l in cv_list for val in l]
    # from scipy.stats.stats import spearmanr
    # print(spearmanr(mean_mod, cv_mod))

    # SNIPPET TO EXPORT CV DF ONLY WHERE SIGNIFICANT HITS WERE DETECTED
    # -------------------------------------
    # col_to_keep = list(norm_agg_nutrient_mean_log2)
    # col_to_keep.insert(0, "Sample Group")
    # agg_nutrient_cv[col_to_keep].to_csv('test_sig.csv', index=False)
