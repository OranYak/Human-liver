# A spatial transcriptomics atlas of live donors reveals unique zonation patterns in the healthy human liver.

We performed 10x Visium spatial transcriptomics on 16 adult human liver samples: eight from young, healthy, living donors (marked as 'M') and eight from patients with liver pathology, where we sampled ‘adjacent normal’ tissue (marked as 'P') to construct a spatial expression atlas of the adult human liver. To achieve higher spatial resolution, we also applied high-definition Visium HD measurements.
To compare human liver zonation profiles with those of other mammals, we assembled a spatial transcriptomics dataset from three additional mammalian species with body sizes and metabolic rates more comparable to human: wild boar (n=2. marked as 'non+human_P'), cow (n=2. marked as 'non+human_C'), and domesticated pig (n=3. marked as 'non+human_PD'). 
Spatial transcriptomics data of the human and non-human samples that were generated in this study can be downloaded from Zenodo, containing preprocessed 10X raw count tables, spot barcodes, tissue images and loupe browser files of all the human and non-human samples using the link https://zenodo.org/records/14795740?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijg5NzU1MzRjLTA2ZjYtNDU1ZS04MWZhLTRhNmIwODBkNDk4YSIsImRhdGEiOnt9LCJyYW5kb20iOiJhZDNmNDdjN2JlMjVmMTRhODZlMGU2NmI3NTE2NzA4NCJ9.5dbmnWYhg8UHjzXemALu5n94drsy3YBKgPALlqbwB1ZLrBEoIyw1vPJbU4cjV6IQS0Qd-QeXdaMfYtFdDXNCYg.

# Prerequisites:
Matlab R2022a

# Installation
Download RAW input data files from Zenodo-
https://zenodo.org/records/14795740?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijg5NzU1MzRjLTA2ZjYtNDU1ZS04MWZhLTRhNmIwODBkNDk4YSIsImRhdGEiOnt9LCJyYW5kb20iOiJhZDNmNDdjN2JlMjVmMTRhODZlMGU2NmI3NTE2NzA4NCJ9.5dbmnWYhg8UHjzXemALu5n94drsy3YBKgPALlqbwB1ZLrBEoIyw1vPJbU4cjV6IQS0Qd-QeXdaMfYtFdDXNCYg

Run the the main scripts at the followind order:
(notice that your directory should include all the Zenodo raw data, and the Matlab_functions directory that is also attached here)
a1_import_data_raw_and_filtered_for_github.m
a2_visium_12_patients_initial_QCs_for_github.m
a3_visium_12_patients_zonation_reconstruction_for_github.m
run time 10-12 min.
This will result on outputiny v.mat, which is also uploaded here.

For ferther analysis:

Pseudo-bulk DE genes analysis: Figure1_pseudobulk_analyses.m

Multi-linear regression analysis: multilinear_regression_analysis.m

Computational pipeline for the mapping of Visium HD spots to single segmented cells: Mapping_spots_to_cells.groovy

Computational pipeline for mapping snRNA-seq hepatocytes and NPCs zonation based on Visium HD: parse_snRNAseq_combined_atlas.m
