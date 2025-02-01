function t=import_visium_data_funcion_for_github(input_path,patient,Z_THRESH_MITO,Z_THRESH_SUM_UMI,HEMO_FRAC,BOUNDARY_PRCTILE,DIST_THRESH)

if nargin<7
    DIST_THRESH=Inf;
end
if nargin<6
    BOUNDARY_PRCTILE=100; % the spot is removed if it is beyond this distance from the spot center x or y 
end
if nargin<5
    %HEMO_FRAC=0.01;
    HEMO_FRAC=inf ; % we don't want to remove spots with high fraction of hemoglobin genes since they often mark vessels
end
if nargin<4
    Z_THRESH_SUM_UMI=-4;
end

if nargin<3
    Z_THRESH_MITO=4;
end

% read spot position
display('Reading spot positions');
T=readtable([input_path,patient,'_tissue_positions_list.csv']);
barcode_space=table2cell(T(:,1));
coor=table2array(T(:,2:end));

% % import the microscopy image
% im=imread([input_path,patient,'_tissue_hires_image.png']);
% t.slide_image=im;

display('Importing UMI counts');
T=readtable([input_path,patient,'_counts_UTT.csv']);
t.spot_name=T.Properties.VariableNames(2:end);
t.gene_name=table2cell(T(:,1));
t.coor=NaN*ones(length(t.spot_name),2);
t.mat=table2array(T(:,2:end));

% prefrom bg substruction 
bg_vector = compute_background_counts_for_github(input_path,patient,'mean');
if ~isempty(bg_vector)
    t.mat = t.mat-bg_vector;
    t.mat(find(t.mat<0)) = 0;
end
t.bg_vector = bg_vector;
t.mat_norm=t.mat./sum(t.mat);

% assign the coordinates of each spot according to the barcode_space
for i=1:length(t.spot_name)
    str=t.spot_name{i};
    str(findstr(str,'_'))='-';
    indd=find(strcmpi(barcode_space,str));
    if ~isempty(indd)
        t.coor(i,:)=coor(indd,end:-1:end-1);
    end
end
