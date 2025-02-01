addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab\functions')
addpath('X:\Common\Lab_matlab_functions\Violinplot-Matlab-master');
addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab');
addpath('X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Human_Liver_Matlab\visium analysis without backgroung');
addpath('X:\Common\Lab_matlab_functions');
addpath('X:\Common\Lab_matlab_functions\gseaDotPlot');
addpath('X:\Yotam\Human\Visium\analysis\Matlab\Code\fun');


%% Reconstruction strategy: for the MAYO tissue establish consensus 20 portal and 20 central LMs based on correlation or anticorrelation with Cyp2e1

%% Examine correlations with prominent zonated landmark genes - cyp2e1 and Sds

clear corr_mat;
clear corr_mat_all;
central_LM='CYP2e1';
portal_LM='SDS';
EXP_THRESH=10^-5;
corr_mat_all=cell(1,length(v));
for i=1:length(v)
    i
    t=v{i};
    if isfield(t,'ind_fib_spots')
        ind_spots=setdiff(1:length(t.spot_name),t.ind_fib_spots);
    else if isfield(t,'ind_capsule_spots')
            ind_spots=setdiff(1:length(t.spot_name),t.ind_capsule_spots);
    else
        ind_spots=1:length(t.spot_name);
    end
    end
    ind_CM=find(strcmpi(t.gene_name,central_LM));
    ind_PM=find(strcmpi(t.gene_name,portal_LM));
    corr_mat=NaN(length(t.gene_name),2); % 2 columns, first is correlation with central_LM, second correlation with portal_LM
    indin=find(nanmean(t.mat_norm,2)>EXP_THRESH);

    for j=1:length(indin);
        corr_mat(indin(j),1)=corr(t.mat_norm(indin(j),ind_spots)',t.mat_norm(ind_CM,ind_spots)','type','spearman');
        corr_mat(indin(j),2)=corr(t.mat_norm(indin(j),ind_spots)',t.mat_norm(ind_PM,ind_spots)','type','spearman');
    end
    corr_mat_all{i}=corr_mat;
end

%% choose canonical LMs
NUM2TAKE=20;
USE_MAX_NORM=1;
SZ=10;
corr_with_cyp2e1_MAYO=zeros(size(corr_mat_all{1},1),1);

index_mayo=[3 4 5 6 7]; %% in order to include only truly healthy mayo samples, and not the steatotic ones.

for i=1:length(index_mayo)
    corr_with_cyp2e1_MAYO=corr_with_cyp2e1_MAYO+corr_mat_all{index_mayo(i)}(:,1);
end
% Periportal LM
index=find(~isnan(corr_with_cyp2e1_MAYO));
[y,ord]=sort(corr_with_cyp2e1_MAYO(index),'ascend');
LM_PP_ind_corr=index(ord(1:NUM2TAKE));
display('portal LMs:')
sort(v{1}.gene_name(LM_PP_ind_corr))
% Pericentral LM
[y,ord]=sort(corr_with_cyp2e1_MAYO(index),'descend');
LM_PC_ind_corr=index(ord(1:NUM2TAKE));
display('central LMs:')
sort(v{1}.gene_name(LM_PC_ind_corr))

central_LM_genes=v{1}.gene_name(LM_PC_ind_corr);
portal_LM_genes=v{1}.gene_name(LM_PP_ind_corr);

% Compute eta score.
% The eta score can be calculated using the mat_norm normelized to max or
% not. We used the normelized to max.

for i=1:length(v)
    t=v{i};
    v{i}.mat_norm_max=v{i}.mat_norm./max(v{i}.mat_norm,[],2);
    v{i}.mat_norm_hep_genes_only_human_max=v{i}.mat_norm_hep_genes_only_human./max(v{i}.mat_norm_hep_genes_only_human,[],2); % added 23.5.2024 to compute eta also according to hep-normalized mat_norm 

    [Ia,Ib]=ismember(lower(central_LM_genes),lower(t.gene_name));
    v{i}.LM_pc_ind=Ib(Ia);
    [Ia,Ib]=ismember(lower(portal_LM_genes),lower(t.gene_name));
    v{i}.LM_pp_ind=Ib(Ia);
    if USE_MAX_NORM
        sum_pp=sum(v{i}.mat_norm_max(v{i}.LM_pp_ind,:));
        sum_pc=sum(v{i}.mat_norm_max(v{i}.LM_pc_ind,:));
        sum_pp_hep_genes=nansum(v{i}.mat_norm_hep_genes_only_human_max(v{i}.LM_pp_ind,:)); % added 23.5.2024
        sum_pc_hep_genes=nansum(v{i}.mat_norm_hep_genes_only_human_max(v{i}.LM_pc_ind,:)); % added 23.5.2024
    else
        sum_pp=sum(v{i}.mat_norm(v{i}.LM_pp_ind,:));
        sum_pc=sum(v{i}.mat_norm(v{i}.LM_pc_ind,:));
        sum_pp_hep_genes=nansum(v{i}.mat_norm_hep_genes_only_human(v{i}.LM_pp_ind,:)); % added 23.5.2024
        sum_pc_hep_genes=nansum(v{i}.mat_norm_hep_genes_only_human(v{i}.LM_pc_ind,:)); % added 23.5.2024
    end
    
    % v{i}.mat_norm is now v{i}.mat./sum(v{i}.mat), so the eta is computed
    % according to normalization to all genes
    v{i}.eta=sum_pp./(sum_pp+sum_pc);    
    v{i}.eta_orig=v{i}.eta;

end

% plot according to peri-central:
figure;
for i=1:length(v)
    t=v{i};    
    nexttile
    [Ia,Ib]=ismember(lower(central_LM_genes),lower(t.gene_name));
    index=Ib(Ia);
    var=sum(t.mat_norm_max(index,:));
    indin=find(var>0);
    scatter(t.coor(:,1),t.coor(:,2),SZ,repmat(0.7,1,3),'filled');
    hold on;
    scatter(t.coor(indin,1),t.coor(indin,2),SZ,var(indin),'filled'); colorbar;
    set(gca,'ydir','reverse');
    set(gca,'XTick',[], 'YTick', []);
    title([v{i}.main_feature 'central LM']);
    axis square
    axis tight
    box on;
end
set(gcf,'name','central LM');
set(gcf,'position',[469         207        1684         806]);

% plot according to peri-portal:
figure;
for i=1:length(v)
    t=v{i};    
    nexttile
    [Ia,Ib]=ismember(lower(portal_LM_genes),lower(t.gene_name));
    index=Ib(Ia);
    var=sum(t.mat_norm_max(index,:));
    indin=find(var>0);
    scatter(t.coor(:,1),t.coor(:,2),SZ,repmat(0.7,1,3),'filled');
    hold on;
    scatter(t.coor(indin,1),t.coor(indin,2),SZ,var(indin),'filled'); colorbar;
    set(gca,'ydir','reverse');
    set(gca,'XTick',[], 'YTick', []);
    title([v{i}.main_feature 'portal LM']);
    axis square
    axis tight
    box on;
end
set(gcf,'name','portal LM');
set(gcf,'position',[469         207        1684         806]);

% plot according to eta score:
figure;
for i=1:length(v)
    t=v{i};    
    nexttile
    scatter(t.coor(:,1),t.coor(:,2),SZ,repmat(0.7,1,3),'filled');
    hold on;
    var=t.eta;
    indin=find(var>0);
    scatter(t.coor(indin,1),t.coor(indin,2),SZ,var(indin),'filled'); colorbar;
    set(gca,'ydir','reverse');
    set(gca,'XTick',[], 'YTick', []);
    title([v{i}.main_feature ' - eta']);
    axis square
    axis tight
    box on;
end
set(gcf,'name','eta');
set(gcf,'position',[469         207        1684         806]);

%% compute zonation

EXP_THRESH=1*10^(-6);
for i=1:length(v)
    disp(['***** computing zonation struct for ' v{i}.patient '*****'] );
    v{i}.zon_struct=extract_zonation_for_github(v{i},EXP_THRESH,8,1,0);
end

for i=1:length(v)
    v{i}.zon_struct.zone_index_orig=v{i}.zon_struct.zone_index;
end

% add a median filter to avoid salt and pepper noise
for i=1:length(v)
    disp(['***** computing zonation struct for ' v{i}.patient '*****'] );
    v{i}=median_zone_filter_for_github(v{i},0);
    v{i}.zon_struct.zone_index=v{i}.zon_struct.zone_index_med;
    v{i}.zon_struct=extract_zonation_from_zone_index(v{i},EXP_THRESH,8,1,0);
end

figure;
for i=1:length(v)
    t=v{i};    
    nexttile
    var=t.zon_struct.zone_index;
    scatter(t.coor(:,1),t.coor(:,2),SZ,var,'filled'); colorbar;
    set(gca,'ydir','reverse');
    %title([str{i} ', ' gg]);
    title([v{i}.main_feature 'zone index']);
    axis square
    set(gca,'XTick',[], 'YTick', []);
end
set(gcf,'name','zone index');
set(gca,'XTick',[], 'YTick', []);
set(gcf,'position',[469         207        1684         806]);

%% Fixing the "zone 0" problem. converting zone 0 to zone 9, which is the fibrotic areas.

for(i=1:length(v))
    zone_zero_ind=find(v{i}.zon_struct.zone_index==0);
    v{i}.zon_struct.zone_index(zone_zero_ind)=9;
end

%% export to category for Loupe
for j=1:length(v)
    t=v{j};
    filename=['X:\oran\Data\Human_Liver_Project\Human_Liver_Visium\Loupe_categories\' t.patient '\zonation_coors_21_11_2024_LM_genes_M2_TO_M8.csv']
    fid=fopen(filename,'w');
    for i=1:length(t.spot_name),
        str=t.spot_name{i};str(findstr(str,'_'))='-';
        fprintf(fid,'%s,%d\n',str,t.zon_struct.zone_index(i));
    end;
    fclose all
end