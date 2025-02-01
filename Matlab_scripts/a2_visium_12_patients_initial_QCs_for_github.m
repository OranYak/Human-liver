addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab\functions')
addpath('X:\Common\Lab_matlab_functions\Violinplot-Matlab-master');
addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab');
addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab\visium analysis without backgroung');
addpath('X:\Common\Lab_matlab_functions');
addpath('X:\Common\Lab_matlab_functions\gseaDotPlot');
addpath('X:\Yotam\Human\Visium\analysis\Matlab\Code\fun');

%% QC - plot sum of umis and mito fraction, filter out spots with low UMI
% Run only once:

num_of_patients=16;
UMI_ZTHRESH=-3;
MITO_FRAC_ZTHRESH= 5;
ttls={'P2-fibrosis','P3-post ALPPS','P6-steatosis','P7-normal','P14-steatosis','P17-severe steatosis','P18-steatosis','P21-fibrosis','M1-normal+steatosis','M2-normal','M3-normal','M4-normal','M5-normal','M6-normal','M7-normal','M8-normal'};
figure;
ind_mito_genes=find(contains(lower(v{1}.gene_name),'mt-'));
for i=1:length(v)
    t=v{i};
    ind_mito_genes=find(contains(lower(v{i}.gene_name),'mt-'));
    % plot sum of umis
    sum_umi=sum(t.mat);
    subplot(4,num_of_patients,i);
    scatter(t.coor(:,1),t.coor(:,2),25,log10(sum_umi),'filled');
    colorbar;
    title([ttls{i} ', sum UMIs']);
    set(gca,'ydir','reverse');
    % plot histogram of Z-scored sum of umis
    Z=(log10(sum_umi+1)-mean(log10(sum_umi+1)))./std(log10(sum_umi+1));
    subplot(4,num_of_patients,i+num_of_patients);
    hist(Z,100);
    hold on;
    line([UMI_ZTHRESH UMI_ZTHRESH],ylim,'color','k');
    colorbar;
    title([ttls{i} ', sum of umis']);
    ind_out_umi=find(Z<UMI_ZTHRESH);
    subplot(4,num_of_patients,i);
    hold on;
    scatter(t.coor(ind_out_umi,1),t.coor(ind_out_umi,2),'r.');

    % plot mitochondrial fraction
    mito_frac=sum(t.mat_norm(ind_mito_genes,:));
    subplot(4,num_of_patients,i+num_of_patients*2);
    scatter(t.coor(:,1),t.coor(:,2),25,mito_frac,'filled');
    colorbar;
    set(gca,'ydir','reverse');
    title([ttls{i} ', mito fraction']);
    % plot histogram of Z-scored mitochondrial fraction
    Z=(log10(mito_frac+1)-nanmean(log10(mito_frac+1)))./nanstd(log10(mito_frac+1));
    subplot(4,num_of_patients,i+num_of_patients*3);
    hist(Z,50);
    hold on;
    line([MITO_FRAC_ZTHRESH MITO_FRAC_ZTHRESH],ylim,'color','k');
    colorbar;
    title([ttls{i} ', mito fraction']);
    ind_out_mito=find(Z>MITO_FRAC_ZTHRESH);
    subplot(4,num_of_patients,i+num_of_patients*2);
    hold on;
    scatter(t.coor(ind_out_mito,1),t.coor(ind_out_mito,2),'r.');

    % discarding the ind_out spots from mat and creating a new mat and a
    % new mat_norm
    ind_out=union(ind_out_mito,ind_out_umi);
    %     ind_out=ind_out_umi;
    ind_out_nan=find(isnan(t.mat_norm(1,:)));
    ind_out=union(ind_out,ind_out_nan);
    t.mat_raw_before_UMI_MITO_filter=t.mat;
    t.ind_out_UMI_MITO_NAN=ind_out;
    t.mat(:,ind_out)=[];
    t.mat_norm=t.mat./sum(t.mat);
    t.spot_name(ind_out)=[];
    t.coor(ind_out,:)=[];
    v{i}=t;
end


%% arranging the v structure in clinical-corresponding order:
% Run only once:

v_temp=v;

v{1}=v_temp{10};
v{2}=v_temp{11};
v{3}=v_temp{12};
v{4}=v_temp{13};
v{5}=v_temp{14};
v{6}=v_temp{15};
v{7}=v_temp{16};
v{8}=v_temp{9};
v{9}=v_temp{4};
v{10}=v_temp{3};
v{11}=v_temp{5};
v{12}=v_temp{7};
v{13}=v_temp{6};
v{14}=v_temp{1};
v{15}=v_temp{8};
v{16}=v_temp{2};

% adding a main_feature for each structure
% ttls={'P2-fibrosis','P3-post ALPPS','P6-steatosis','P7-normal','P14-steatosis','P17-severe steatosis','P18-steatosis','P21-fibrosis','M1-normal+steatosis','M2-normal','M3-normal','M4-normal'};
v{1}.main_feature='M2-normal';
v{2}.main_feature='M3-normal';
v{3}.main_feature='M4-normal';
v{4}.main_feature='M5-normal';
v{5}.main_feature='M6-normal';
v{6}.main_feature='M7-normal';
v{7}.main_feature='M8-normal';
v{8}.main_feature='M1-normal & steatosis';
v{9}.main_feature='P7-minimal steatosis & capsule';
v{10}.main_feature='P6-steatosis';
v{11}.main_feature='P14-steatosis';
v{12}.main_feature='P18-steatosis';
v{13}.main_feature='P17-severe steatosis & fibrosis';
v{14}.main_feature='P2-fibrosis';
v{15}.main_feature='P21-fibrosis';
v{16}.main_feature='P3-post ALPPS';

% str={'P7-minimal steatosis','P6-steatosis','P14-steatosis','P18-steatosis','P17-severe steatosis','P2-fibrosis','P21-fibrosis','P3-post ALPPS','M1-normal+steatosis','M2-normal','M3-normal','M4-normal'};

v_temp=v;

%% find fibrotic zones (for now-P2,P17,P21):
% Run only once:
for i=1:length(v)
    t=v{i};
    fibrotic_spots=[];
    ind_fibrotic_spots=[];
    input_path=['X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Loupe_categories\' v{i}.patient '\fibrotic_zones.csv'];
    if isfile(input_path)
        fibrotic_spots=read_loupe_barcode(input_path);
    else
        continue;
    end
    for j=1:size(fibrotic_spots,1)
        indd=find(strcmpi(t.spot_name,fibrotic_spots(j)));
        if ~isempty(indd)
            ind_fibrotic_spots=[ind_fibrotic_spots;indd];
        end
    end
    figure;
    scatter(t.coor(:,1),t.coor(:,2),20,[0.6 0.6 0.6],'filled');
    hold on;
    scatter(t.coor(ind_fibrotic_spots,1),t.coor(ind_fibrotic_spots,2),'rx');
    title([v{i}.patient ' fibrotic spots from Loupe']);
    set(gca,'ydir','reverse');
    v{i}.ind_fib_spots=ind_fibrotic_spots;
end

% find capsule zones of (for now-P7,M8):
for i=1:length(v)
    t=v{i};
    capsule_spots=[];
    ind_capsule_spots=[];
    input_path=['X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Loupe_categories\' v{i}.patient '\capsule_zones.csv'];
    if isfile(input_path)
        capsule_spots=read_loupe_barcode(input_path);
    else
        continue;
    end
    for j=1:size(capsule_spots,1)
        indd=find(strcmpi(t.spot_name,capsule_spots(j)));
        if ~isempty(indd)
            ind_capsule_spots=[ind_capsule_spots;indd];
        end
    end
    figure; scatter(t.coor(:,1),t.coor(:,2),20,[0.6 0.6 0.6],'filled');
    hold on;
    scatter(t.coor(ind_capsule_spots,1),t.coor(ind_capsule_spots,2),'rx');
    title([v{i}.patient ' capsule spots from Loupe']);
    set(gca,'ydir','reverse');
    v{i}.ind_capsule_spots=ind_capsule_spots;
end

% find lipid zones of (for now-P6,P7,P14,P18,P17,M1):
for i=1:length(v)
    t=v{i};
    lipid_spots=[];
    ind_lipid_spots=[];
    input_path=['X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Loupe_categories\' v{i}.patient '\lipid_zones.csv'];
    if isfile(input_path)
        lipid_spots=read_loupe_barcode(input_path);
    else
        continue;
    end
    for j=1:size(lipid_spots,1)
        indd=find(strcmpi(t.spot_name,lipid_spots(j)));
        if ~isempty(indd)
            ind_lipid_spots=[ind_lipid_spots;indd];
        end
    end
    figure; scatter(t.coor(:,1),t.coor(:,2),20,[0.6 0.6 0.6],'filled');
    hold on;
    scatter(t.coor(ind_lipid_spots,1),t.coor(ind_lipid_spots,2),'rx');
    title([v{i}.patient ' lipid spots from Loupe']);
    set(gca,'ydir','reverse');
    v{i}.ind_lipid_spots=ind_lipid_spots;
end