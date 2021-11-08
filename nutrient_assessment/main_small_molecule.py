import os

from dataproc import read_csv
import dataproc.untargeted_ms as ms
import main_untargeted as mu

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
    (
        agg_nutrient_mean,
        agg_nutrient_std,
        agg_nutrient_cv,
    ) = mu.aggregate_mean_std_cv_from_nutrient_data(small_molecule_df)

    # Normalize aggregated mean data to mean of control - perform log2 scaling of results
    norm_agg_nutrient_mean = mu.normalize_nutrient_data_to_control(agg_nutrient_mean)
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
