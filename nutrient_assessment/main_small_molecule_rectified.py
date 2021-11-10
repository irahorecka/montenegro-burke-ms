import os

from dataproc import read_csv
import dataproc.untargeted_ms as ms
import main_untargeted_rectified as mu

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
FIGURES_PATH = os.path.join(BASE_PATH, "figures")


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

    small_molecule_df = mu.isolate_cols_and_transpose_df(
        small_molecule_df, ["DetectedMass", "Area"]
    )
    small_molecule_df = mu.mutate_and_relabel_nutrient_data(
        small_molecule_df, src_colname="DetectedMass"
    )
    small_molecule_df = mu.test(small_molecule_df, 6)
    (
        agg_nutrient_mean,
        agg_nutrient_std,
        agg_nutrient_cv,
    ) = mu.aggregate_mean_std_cv_from_nutrient_data(small_molecule_df)

    # Normalize aggregated mean data to mean of control - perform log2 scaling of results
    norm_agg_nutrient_mean = mu.normalize_nutrient_data_to_control(agg_nutrient_mean)
    norm_agg_nutrient_mean_log2 = ms.get_log2_df(norm_agg_nutrient_mean, log2_weight=1)
    norm_agg_nutrient_mean_log2 = ms.get_df_values_within_range(norm_agg_nutrient_mean_log2, -5, 5)

    # Export data as CSV for further analysis
    norm_agg_nutrient_mean_log2.to_csv(
        os.path.join(DATA_PATH, "log2_nutrient_mean_small_molecule.csv")
    )
