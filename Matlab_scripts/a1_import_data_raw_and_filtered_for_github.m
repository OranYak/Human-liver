%% Info
% This is a  Master script to call all functions needed to align diffrent layer of
% data and produce a visium structure for downstream analysis and figure creation.
% The raw data is available at Zenodo:
% https://zenodo.org/records/14783760?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjhjOTc0ZTc2LTNiN2YtNDg4OC04NDhmLTk3YTUxODkwYTE2YyIsImRhdGEiOnt9LCJyYW5kb20iOiJhZWQ3YWFlYWU1NGI2YTljZjUyZGY5Yzc0YWM4MTQ3NCJ9.cZvQN0jiLqdD7fqAa_6cRNINc4ywSJOogNcw2lgdUYdekF88Q9xCk7xqFLxOnkFzSJjw7Fw2127DB4N8mTIjuQ

%% add the path to all the necessary functions and define the working directory
current_dir = cd; % In this folder you should have all  the Zenodo raw data, and the Matlab_functions folder from Github
addpath([current_dir,'\Matlab_functions\']);

% define the fixed parameters for all samples (same as function defult):

Z_THRESH_MITO = 4;
Z_THRESH_SUM_UMI = -4;
HEMO_FRAC = inf;
BOUNDARY_PRCTILE = 100;
DIST_THRESH = inf;

% set the path for the files:
input_path = ([current_dir,'\human_raw_16_samples\']);
pati = {'P2','P3','P6','P7','P14','P17','P18','P21','M1','M2','M3','M4','M5','M6','M7','M8'};
metadata = readtable([input_path,'metadata.xlsx']);

%% loop over all patients

v=cell(1,(length(pati)));
for i=1:length(v)
    disp(upper(['***** importing visium data for ',pati{i},'*****']));
    v{i} = import_visium_data_function_for_github(input_path,pati{i},...
                                    Z_THRESH_MITO,Z_THRESH_SUM_UMI,HEMO_FRAC,...
                                    BOUNDARY_PRCTILE,DIST_THRESH);
    v{i}.metadata = metadata(i,:);
    v{i}.patient = pati{i};
    %v{i}.coor = v{i}.coor .* optional: v{i}.json.tissue_hires_scalef; % asjust coordianted to the slide image    
end

[Ia,Ib]=ismember(sig.gene_name(ind_markers),v{1}.gene_name);
ind_hep_human=Ib(Ia);

for i=1:length(v)
    disp(upper(['***** computing mat norm over hepatocyte genes for ',pati{i},'*****']));
    % v{i}.mat_norm is already established as v{i}.mat./sum(v{i}.mat)
    v{i}.mat_norm_hep_genes_5_species=v{i}.mat./sum(v{i}.mat(ind_hep_5_species,:)); % ind_hep_5_species is 1344 genes
    v{i}.ind_hep_5_species=ind_hep_5_species;
    v{i}.mat_norm_hep_genes_only_human=v{i}.mat./sum(v{i}.mat(ind_hep_human,:)); % ind_hep_human is 1732 genes
    v{i}.ind_hep_human=ind_hep_human;
end

disp(upper(['***** finished importing visium data for ',num2str(length(v)),' patients *****']));


%% export to R the bg vector:
output_dir = 'X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab\updated_scripts\data\processed\';
for i =1:length(v)
    t = v{i};
    writetable(table(t.gene_name,t.bg_vector),[output_dir,t.patient,'_background_vec.csv']);
end   
