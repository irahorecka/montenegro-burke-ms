import numpy as np
import pandas as pd


def get_df_values_within_range(df, min, max):
    df = df[df.fillna(0) < max]
    return df[df.fillna(0) > min].dropna(axis=1)


def get_empty_df_from_df(df):
    return pd.DataFrame(columns=list(df))


def concat_df(*args, axis=0):
    return pd.concat(list(args), axis=axis)


def get_df_with_cols_to_keep(df, col_substrings=None):
    if not col_substrings:
        return df
    cols_to_keep = []
    for col_substring in col_substrings:
        cols_to_keep.extend([colname for colname in list(df) if col_substring in colname])
    return df[cols_to_keep]


def convert_value_to_nan(df, value):
    return df.replace(value, np.nan, inplace=False)


def drop_rows_with_substring_in_col_value(df, colname, substring):
    return df[~df[colname].str.contains(substring)]


def get_cols_with_less_than_count_in_row(df, colname, count):
    # Drop rows containing values less than count.
    return df[(df[colname] < count)].T


def transpose_and_reset_idx(df):
    return df.T.reset_index()


def mv_row_as_header(df, row_idx):
    """Moves row as dataframe header of df indicated by `row_idx`.
    Drops row name row once complete."""
    return df.rename(columns=df.iloc[row_idx]).drop(df.index[row_idx])


def convert_to_numerics(df):
    """Converts object dtype to numeric, if possible."""
    return df.apply(pd.to_numeric, errors="ignore")


def map_dict_and_replace_values(df, colname, dict_map):
    return df.replace({colname: dict_map})


def map_ternary_exp_and_replace_values(df, ternary_cond, new_colname):
    """Assesses ternary pattern (`ternary_cond`) and maps to `np.where`.
    Assigns output column to `new_colname`."""
    df[new_colname] = np.where(*ternary_cond)
    return df


def group_and_agg(df, colname, agg_type):
    agg_type_map = {
        "mean": np.nanmean,
        "std": np.nanstd,
    }
    if isinstance(agg_type, list):
        return df.groupby([colname]).agg(agg_type)
    return df.groupby([colname]).agg(agg_type_map.get(agg_type, agg_type))


def get_log2_df_directional(df, downregulated=False, log2_weight=1):
    """Get df where any column values greater than log2_weight will be kept
    if downregulated=False. Less than log2_weight will be kept if downregulated=True."""
    # Convert values to log2 - get a sense of upregulation and downregulation
    df_log = np.log2(df)
    # Remove inf / -inf and replace NaN with 0
    df_log = convert_value_to_nan(df_log, [np.inf, -np.inf])

    # Refine values to query - improves resolution of large differences
    if downregulated:
        df_log = df_log[df_log.fillna(0) < -log2_weight]
    else:
        df_log = df_log[df_log.fillna(0) > log2_weight]

    df_log = df_log.dropna(axis=1, how="all")
    header_colnames_to_keep = list(df_log)
    return np.log2(df[header_colnames_to_keep])


def get_log2_df(df, log2_weight=1):
    """Get df where any column values greater than log2_weight or less than
    log2_weight will be kept."""
    # Convert values to log2 - get a sense of upregulation and downregulation
    df_log = np.log2(df)
    # Remove inf / -inf and replace NaN with 0
    df_log = convert_value_to_nan(df_log, [np.inf, -np.inf])

    # Refine values to query - improves resolution of large differences
    df_log_downregulated = df_log[df_log.fillna(0) < -log2_weight]
    df_log_upregulated = df_log[df_log.fillna(0) > log2_weight]
    df_log_downregulated[df_log_downregulated.isnull()] = df_log_upregulated
    df_log = df_log_downregulated

    df_log = df_log.dropna(axis=1, how="all")
    header_colnames_to_keep = list(df_log)
    return np.log2(df[header_colnames_to_keep])
