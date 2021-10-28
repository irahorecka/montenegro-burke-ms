from dataproc import read_csv, rm_empty_cols, sort_abundance_below_precursor


if __name__ == "__main__":
    path = "/Users/irahorecka/Desktop/Harddrive_Desktop/PhD/University of Toronto/Rotations/Montenegro-Burke Lab/ms2-fragpeak-analysis/data/2021-10-27_MSMS_FragAnalysis/ATP.csv"
    df = read_csv(path)
    df_non_nan = rm_empty_cols(df)
    print(
        sort_abundance_below_precursor(
            df_non_nan,
            abundance_colname="Abund",
            precursor_colname="m/z",
            precursor_value=880,
            ascending=False,
        )
    )
