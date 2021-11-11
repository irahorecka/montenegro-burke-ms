# montenegro-burke-ms
Data analysis repository for rotation with the Rafa Montenegro-Burke Lab (Fall 2021).

# Background
## /fragmentation_peaks
Directory that handles extraction and analysis of metabolite fragmentation peaks from an LC-MS QTOF instrument. Of the metabolites found and ran on the instrument, 20 were extracted successfully using the Agilent Qualitative Analysis B.07.00 software.

The goal is to create a simple summary file of the extracted peaks, sorting each fragmentation peak as relative abundance to other peaks within the same sample. If done correctly, this should provide a simple, curated output to make the process of subsequent fragmentation peak database querying and analysis facile.

## /nutrient_assessment
Directory that handles extraction and analysis of metabolite peaks from an S. cerevisiae nutrient array experiment using data acquired from an LC-MS qTOF instrument and timsTOF instrument. For results from the LC-MS qTOF instrument, the output metabolite AUC, AUC quality score, and m/z mass were extracted from Agilent MassHunter Profinder B.06.00 software.

The goal is to create hierarchically clustered dendrograms representing dysregulated metabolites from the control S. cerevisiae nutrient sample (i.e., media with glucose and ammonia). Those metabolites that show confident discrepancies of metabolite levels in different growth media will be isolated and represented in the hierarchical clusters.

View /nutrient_assessment/data for pre- and post-processed data (all CSV format).

View /nutrient_assessment/figures for output figures (currently comprised of hierarchically clustered dendrograms). View in folders /figures/pdf and figures/png, for PDF and PNG files, respectively.

Data processing and manipulation is handled using Numpy and Pandas using the Python3 language. Plotting of the hierarchically clustered dendrogram is accomplished using R.
