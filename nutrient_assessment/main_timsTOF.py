import os

from dataproc import read_csv
import dataproc.untargeted_ms as ms

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")


def aggregate_mean_std_cv_from_nutrient_data(df, agg_colname="Sample Group"):
    """Aggregates dataframe to find mean, std, and cv, grouping by column
    bound to `agg_colname`. Returns all three aggregated dataframes to caller."""
    df = ms.convert_to_numerics(df)
    df = filter_data_with_more_than_n_reads_among_4_samples(df, agg_colname, n=2)
    # Begin aggregation
    df_mean = ms.group_and_agg(df, colname=agg_colname, agg_type="mean")
    df_std = ms.group_and_agg(df, colname=agg_colname, agg_type="std")
    df_cv = ms.convert_to_numerics(df_std.div(df_mean).reset_index())
    return df_mean, df_std, df_cv


def filter_data_with_more_than_n_reads_among_4_samples(df, agg_colname, n=0):
    """Filters and drops metabolite samples that have more than (4 - n) NaN value in any
    of the nutrient conditions. I.e., there must be greater than or equal to n/4 valid samples in
    all nutrient conditions when assessing a particular metabolite."""
    # Generate aggregated data to count NaN
    df_count_valid_data = ms.group_and_agg(df, colname=agg_colname, agg_type="count").drop(
        index=["BLANK", "CTRL"]
    )

    # Find columns to drop that have less than n/4 reads (i.e., greater than (4 - n) NA value in any nutrient category)
    cols_to_drop = []
    for colname in list(df_count_valid_data.T):
        cols_to_drop.extend(
            list(ms.get_cols_with_less_than_count_in_row(df_count_valid_data.T, colname, n))
        )
    # Isolate only unique compounds found to have poor resolution
    cols_to_drop_name = list(set(cols_to_drop))
    return df.drop(cols_to_drop_name, axis=1)


def normalize_nutrient_data_to_control(df):
    """Normalizes nutrient data value to the control (GLC | AMN). Drops the control row
    after normalization and returns dataframe to caller."""
    # Remove blank row and CTRL
    df = df.drop(index=["BLANK", "CTRL"])
    # Divide groups 1-6 by group 3 (i.e., our control) to get relative expression diff
    df.loc[:, df.columns[0] :] = df.loc[:, df.columns[0] :].div(df.iloc[3][df.columns[0] :])
    # Drop control row (Glucose and ammonia)
    df = df.drop(index="GLC | AMN")
    # Re-apply float data type
    return df.astype(float)


if __name__ == "__main__":
    # Get pandas df from CSV path
    tims_path = os.path.join(DATA_PATH, "20211104_IH_timsTOF_Experiment.csv")
    tims_df = read_csv(tims_path)
    tims_df = ms.mv_row_as_header(tims_df, row_idx=0)
    tims_df["Sample Group"] = tims_df.iloc[:, 2].astype(str) + "_" + tims_df.iloc[:, 1].astype(str)

    # Trim, transpose, and move metabolite ID as header for timsTOF data
    tims_df_trunc = ms.transpose_and_reset_idx(tims_df.iloc[:, 5:])

    # Map and change index values
    group_sample_to_nutrient = {
        "GLC _ ASP": "GLC | ASP",
        "GLC _ GLN": "GLC | GLN",
        "GLC _ AMN": "GLC | AMN",
        "GAL _ ASP": "GAL | ASP",
        "GAL _ GLN": "GAL | GLN",
        "GAL _ AMN": "GAL | AMN",
    }
    tims_df_trunc = ms.map_dict_and_replace_values(
        tims_df_trunc, colname="index", dict_map=group_sample_to_nutrient
    )
    tims_df_trunc = ms.mv_row_as_header(tims_df_trunc, row_idx=28).set_index("Sample Group")

    # Convert "0.0" readings to NaN
    tims_df_trunc = ms.convert_value_to_nan(tims_df_trunc, "0.0")

    # Begin aggregating and normalizing aggregated results to control value (GLC | AMN)
    agg_nutrient_mean, _, _ = aggregate_mean_std_cv_from_nutrient_data(tims_df_trunc)
    norm_agg_nutrient_mean = normalize_nutrient_data_to_control(agg_nutrient_mean).dropna(axis=1)
    norm_agg_nutrient_mean_log2 = ms.get_log2_df_directional(
        norm_agg_nutrient_mean, downregulated=False, log2_weight=3
    )

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(os.path.join(DATA_PATH, "log2_nutrient_mean_timsTOF.csv"))
