# A spatial transcriptomics atlas of live donors reveals unique zonation patterns in the healthy human liver.

We performed 10x Visium spatial transcriptomics on 16 adult human liver samples: eight from young, healthy, living donors (marked as 'M') and eight from patients with liver pathology, where we sampled ‘adjacent normal’ tissue (marked as 'P') to construct a spatial expression atlas of the adult human liver. To achieve higher spatial resolution, we also applied high-definition Visium HD measurements to two liver sections from young, healthy, living donors (M2 and M6).
To compare human liver zonation profiles with those of other mammals, we assembled a spatial transcriptomics dataset from three additional mammalian species with body sizes and metabolic rates more comparable to human: wild boar (n=2. marked as 'non+human_P'), cow (n=2. marked as 'non+human_C'), and domesticated pig (n=3. marked as 'non+human_PD'). 
For human Visium samples (n=16. M1, M2, M3, M4, M5, M6, M7, M8, P2, P3, P6, P7, P14, P17, P18 and P21) and non-human Visium samples (n=7. C1, C2, P1, P2, PD1, PD2 and PD3) the following files were uploaded:

Metadata:
Human samples metadata file.
Non-human metadata file.

For each Visium sample individually:
counts_ALL.csv
counts_UTT.csv
scalefactors_json.json
tissue_hires_image.png
tissue_positions_list.csv
cloupe file
filtered_feature_bc_matrix.h5
raw_feature_bc_matrix.h5

For M2 and M6 Visium HD samples:
cloupe file

# Prerequisites:
Matlab R2023a

# Installation
Download RAW input data files from Zenodo-
https://zenodo.org/records/14783760?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjhjOTc0ZTc2LTNiN2YtNDg4OC04NDhmLTk3YTUxODkwYTE2YyIsImRhdGEiOnt9LCJyYW5kb20iOiJhZWQ3YWFlYWU1NGI2YTljZjUyZGY5Yzc0YWM4MTQ3NCJ9.cZvQN0jiLqdD7fqAa_6cRNINc4ywSJOogNcw2lgdUYdekF88Q9xCk7xqFLxOnkFzSJjw7Fw2127DB4N8mTIjuQ

Run the the main scripts at the followind order:
(notice that your directory should include all the Zenodo raw data, and the Matlab_functions directory that is also attached here)
a1_import_data_raw_and_filtered_for_github.m
a2_visium_12_patients_initial_QCs_for_github.m
a3_visium_12_patients_zonation_reconstruction_for_github.m
run time 10-12 min.
