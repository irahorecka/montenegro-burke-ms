import pandas as pd


def read_csv(path):
    return pd.read_csv(path, delimiter=",")


def rm_empty_cols(df):
    df.dropna(how="all", axis=1, inplace=True)
    return df


def sort_abundance_below_precursor(
    df, abundance_colname, precursor_colname, precursor_value, **kwargs
):
    if not isinstance(precursor_value, (int, float)):
        raise ValueError(
            f"precursor_value requires type int or float. Found {type(precursor_value)}"
        )
    sorted_df = sort_df_by_colname(df, abundance_colname, **kwargs)
    return sorted_df[sorted_df[precursor_colname] < precursor_value]


def sort_df_by_colname(df, colname, ascending=True):
    return df.sort_values(by=[colname], ascending=ascending)
