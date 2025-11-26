load combined_scRNAseq_atlas_M5M6M7M8 
%% show the map with clusters
SHOW_TEXT=0;
figure;
scatter(t.umap(:,1),t.umap(:,2),10,t.leiden,'filled',...
    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.6);
colormap jet
if SHOW_TEXT
    for i=1:length(t.cluster_unique)
        ind=find(strcmpi(t.cluster_annotations,t.cluster_unique{i}));
        text(median(t.umap(ind,1)),median(t.umap(ind,2)),string(t.cluster_unique(i)),'fontweight','b','color','r','fontsize',16);
    end
end
axis square
axis tight
set(gca,'xtick',[],'ytick',[]);
set(gcf,'position',[680   339   880   656])
box on
axis off

t.cluster_names=t.cluster_unique;

%% show the map with cluster names
index=1:size(t.mat,2);
figure;
SZ=15;
scatter(t.umap(index,1),t.umap(index,2),SZ,t.leiden(index),'filled');
for i=1:length(t.cluster_unique)
    ind=find(strcmpi(t.cluster_annotations,t.cluster_unique{i}));
    ind=intersect(ind,index);
    text(median(t.umap(ind,1)),median(t.umap(ind,2)),t.cluster_names{i},'fontweight','b',...
        'color','r','fontsize',14);
end
colorbar
colormap('jet')
set(gcf,'position',[1000         124        1281        1114]);
box on;
xlabel('UMAP1')
ylabel('UMAP2')
set(gca,'fontsize',18)

t.zone=NaN(size(t.mat,2),1);

%% Extract zonation
% read LM genes from the VisiumHD dataset
cell_types={'LSEC','Kupffer','Fibroblasts','Hepatocytes'};
figure;
t.zone=NaN(size(t.mat,2),1);
COMPUTE_PVALS=1;
MIN_UMIS=500;
t.zone=NaN(size(t.mat,2),1);
ALB_THRESH=0.01;

for ii=1:length(cell_types)
    ind_cells=find(strcmpi(t.cluster_annotations,cell_types{ii}));
    ind_cells=intersect(ind_cells,find(sum(t.mat)>MIN_UMIS));
    if ~strcmpi(cell_types{ii},'Hepatocytes')
        ind_cells=intersect(ind_cells,find(t.mat_norm(find(strcmpi(t.gene_name,'ALB')),:)<ALB_THRESH));
    end
    length(ind_cells)
    if ii==4,
        cell_types{ii}='Hepatocyte';
    end
    PC_LM=readtable([cell_types{ii} '-PC-LM.csv'],'ReadVariableNames',false);
    [Ia,Ib]=ismember(table2cell(PC_LM),t.gene_name);

    t.mat_norm_max=t.mat_norm;%./max(t.mat_norm(:,ind_cells),[],2);

    sum_pc=sum(t.mat_norm_max(Ib(Ia),:));

    PP_LM=readtable([cell_types{ii} '-PP-LM.csv'],'ReadVariableNames',false);
    [Ia,Ib]=ismember(table2cell(PP_LM),t.gene_name);

    sum_pp=sum(t.mat_norm_max(Ib(Ia),:));
    eta=sum_pp./(sum_pp+sum_pc);


    NUM_BINS=8;
    edges=prctile(eta(ind_cells),linspace(0,100,NUM_BINS+1));
    for i=1:NUM_BINS
        index=ind_cells(find(eta(ind_cells)>=edges(i) & eta(ind_cells)<=edges(i+1)));
        t.zone(index)=i;
    end
    
    t.zonation{ii}.mn=zeros(length(t.gene_name),NUM_BINS);
    t.zonation{ii}.se=zeros(length(t.gene_name),NUM_BINS);
    for i=1:NUM_BINS
        ind=ind_cells(find(t.zone(ind_cells)==i));
        t.zonation{ii}.mn(:,i)=nanmean(t.mat_norm(:,ind),2);
        t.zonation{ii}.se(:,i)=nanstd(t.mat_norm(:,ind),[],2)/sqrt(length(ind));
    end

    t.zonation{ii}.mx=nanmax(t.zonation{ii}.mn,[],2);
    t.zonation{ii}.com=nansum(repmat(1:8,length(t.gene_name),1).*t.zonation{ii}.mn,2)./nansum(t.zonation{ii}.mn,2);

    nexttile;
    scatter(log10(t.zonation{ii}.mx),t.zonation{ii}.com,'.');
    ind1=find(t.zonation{ii}.mx>1e-4);
    [y,ord]=sort(t.zonation{ii}.com(ind1),'ascend');
    ind1=ind1(ord(1:20));
    text(log10(t.zonation{ii}.mx(ind1)),t.zonation{ii}.com(ind1),...
        t.gene_name(ind1));
        ind1=find(t.zonation{ii}.mx>1e-4);
    [y,ord]=sort(t.zonation{ii}.com(ind1),'descend');
    ind1=ind1(ord(1:20));
    text(log10(t.zonation{ii}.mx(ind1)),t.zonation{ii}.com(ind1),...
        t.gene_name(ind1));
    xlabel('log10(max)')
    ylabel('COM');
    title(cell_types{ii})
    if COMPUTE_PVALS
        t.zonation{ii}.pval=NaN(length(t.gene_name),1);
        t.zonation{ii}.qval=NaN(length(t.gene_name),1);
        for k=1:length(t.gene_name)
            t.zonation{ii}.pval(k)=kruskalwallis(t.mat_norm(k,ind_cells),...
                t.zone(ind_cells),'off');
        end
        t.zonation{ii}.qval=mafdr( t.zonation{ii}.pval,'BHFDR','true');
    end
    t.zonation{ii}.cell_type=cell_types{ii};
end

%% show zonation figure
figure;
scatter(t.umap(:,1),t.umap(:,2),15,repmat(0.8,1,3),'filled')
hold on;
ind=find(~isnan(t.zone));
scatter(t.umap(ind,1),t.umap(ind,2),15,t.zone(ind),'filled',...
    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.7)
for i=1:length(t.cluster_unique)
    ind=find(strcmpi(t.cluster_annotations,t.cluster_unique{i}));
    %ind=intersect(ind,index);
    text(median(t.umap(ind,1)),median(t.umap(ind,2)),t.cluster_names{i},'fontweight','b',...
        'color','k','fontsize',16);
end
box on;
axis tight
xlabel('UMAP1');
ylabel('UMAP2');
set(gca,'fontsize',14)
colorbar
set(gcf,'position',[ 1          41        2560        1323]);


%% combine
comb.gene_name=t.gene_name;
for i=1:4,
    comb.cell_type{i}=t.zonation{i}.cell_type;
    comb.pval{i}=t.zonation{i}.pval;
    comb.qval{i}=t.zonation{i}.qval;
    mx=max(t.zonation{i}.mn,[],2);    

    comb.mean_mat{i}=t.zonation{i}.mn;
    comb.se_mat{i}=t.zonation{i}.se;

    comb.com{i}=sum(repmat(1:8,size(comb.mean_mat{i},1),1).*comb.mean_mat{i},2)./sum(comb.mean_mat{i},2);

    T=table(comb.gene_name,comb.mean_mat{i}(:,1),comb.mean_mat{i}(:,2),comb.mean_mat{i}(:,3),comb.mean_mat{i}(:,4),...
        comb.mean_mat{i}(:,5),comb.mean_mat{i}(:,6),comb.mean_mat{i}(:,7),comb.mean_mat{i}(:,8),...
        comb.se_mat{i}(:,1),comb.se_mat{i}(:,2),comb.se_mat{i}(:,3),comb.se_mat{i}(:,4),...
        comb.se_mat{i}(:,5),comb.se_mat{i}(:,6),comb.se_mat{i}(:,7),comb.se_mat{i}(:,8),...
        comb.pval{i},comb.qval{i},comb.com{i},max(comb.mean_mat{i},[],2))
    T.Properties.VariableNames{1}='Gene_Name';
    for ii=1:8
        T.Properties.VariableNames{ii+1}=['Mean_Zone_' num2str(ii)];
        T.Properties.VariableNames{9+ii}=['SE_Zone_' num2str(ii)];
    end
    T.Properties.VariableNames{18}='pValue';
    T.Properties.VariableNames{19}='qValue';
    T.Properties.VariableNames{20}='Center_of_Mass';
    T.Properties.VariableNames{21}='Max_expression';
end


%% After running previous block for hepatocytes, kupffer cells, LSEC and fibroblasts, now examine the zonated interactions with them
% - for each of the 16 pairs, compute the correlation of the zonation
% profiles of ligands and matching receptors, record the highest ones

zon_ttls={'LSEC','Kupffer','Fibroblasts','Hepatocytes'};
[A,B] = meshgrid(zon_ttls, zon_ttls);
pairs_ttls = strcat(A(:),'-',B(:));
EXP_THRESH=0;
load("X:\Common\useful_datasets\ramilowsky.mat");
% arrange pairs of ligands and matching receptors in our data
pairs=[];
for i=1:length(ramilowsky.ligands)
    indL=find(strcmpi(comb.gene_name,ramilowsky.ligands{i}));
    indR=find(strcmpi(comb.gene_name,ramilowsky.receptors{i}));
    if ~isempty(indL) & ~isempty(indR)
        pairs=[pairs;indL indR];
    end
end
% compute  interactions
pairs_potential=NaN(size(pairs,1),length(comb.mean_mat).^2);
pairs_potential_pval=NaN(size(pairs,1),length(comb.mean_mat).^2);

counter=1;
for i=1:length(comb.mean_mat)
    mat1=comb.mean_mat{i};
    mx1=nanmax(mat1,[],2);
    for j=1:length(comb.mean_mat)
        [i j]
        mat2=comb.mean_mat{j};
        mx2=nanmax(mat2,[],2);
         for k=1:size(pairs,1)
               [pairs_potential(k,counter),pairs_potential_pval(k,counter)]=...
                corr(mat1(pairs(k,1),:)',mat2(pairs(k,2),:)','type','spearman');           
        end
        counter=counter+1;
    end
end

%% Organize results
comb.cell_type{4}='Hepatocytes';
clear Tstruct;
for i=1:size(pairs_potential,2)
    T=cell2table(comb.gene_name(pairs),'Variablenames',{'Source','Target'});
    poten_mat=pairs_potential(:,i);
    poten_mat_pval=pairs_potential_pval(:,i);
    strs=strsplit(pairs_ttls{i},'-');
    cell_ind1=find(strcmpi(comb.cell_type,strs{1}));
    cell_ind2=find(strcmpi(comb.cell_type,strs{2}));
    com1=comb.com{cell_ind1}(pairs(:,1));
    com2=comb.com{cell_ind2}(pairs(:,2));
    mx1=max(comb.mean_mat{cell_ind1},[],2);
    mx1=mx1(pairs(:,1));
    mx2=max(comb.mean_mat{cell_ind2},[],2);
    mx2=mx2(pairs(:,2));
    poten_mat_qval=mafdr(poten_mat_pval,'bhfdr','true');
    Tnum=table(poten_mat,com1,com2,mx1,mx2,poten_mat_pval,poten_mat_qval,...
        'VariableNames',{'Spearman correlation','COM Source','COM Target',...
        'Max Source','Max Target','p-value','q-value'});
    T=[T Tnum];
    writetable(T,'zonated_interactions.xlsx','Sheet',pairs_ttls{i});
    Tstruct{i}=T;
end


%% Panels 4c-f zonation maps
EXP_THRESH=5e-4;
QTHRESH=0.1;
NUM2SHOW=80;
figure;
for i=1:length(comb.cell_type)
   
    indin=find(max(comb.mean_mat{i},[],2)>EXP_THRESH & ...
        comb.qval{i}<QTHRESH);
    [y,ord]=sort(comb.qval{i}(indin));
    indin=indin(ord(1:NUM2SHOW));
    [y,ord]=sort(comb.com{i}(indin));
    subplot(1,4,i);
    mat=comb.mean_mat{i}(indin(ord),:);
    mat=mat./max(mat,[],2);
    imagesc(mat);
    set(gca,'ytick',1:length(ord),'yticklabel',...
        comb.gene_name(indin(ord)));
    set(gca,'xtick',[1 8],'xticklabel',{'Central','Portal'});
    for j=1:length(ord)
        yline(j-0.5);
    end
    for j=1:8
        xline(j-0.5);
    end
    title(comb.cell_type{i});
    colormap redbluecmap
    set(gca,'fontsize',7)
end
set(gcf,'position',[         348         -11        1238        1006])

%% plot a pair of genes
clear pairs;
clear cell_pairs;
pairs{1}={'WNT2','FZD5'};cell_pairs{1}={'LSEC','Hepatocytes'};
pairs{2}={'WNT2','FZD1'};cell_pairs{2}={'LSEC','Fibroblasts'};
pairs{3}={'RSPO3','LGR5'};cell_pairs{3}={'Fibroblasts','Hepatocytes'};
pairs{4}={'BMP10','ACVR1'};cell_pairs{4}={'Fibroblasts','Hepatocytes'};
pairs{5}={'GDF2','ENG'};cell_pairs{5}={'Fibroblasts','LSEC'};
pairs{6}={'EFNB3','EPHB4'};cell_pairs{6}={'Fibroblasts','LSEC'};
pairs{7}={'LAMA5','SDC1'};cell_pairs{7}={'Fibroblasts','Hepatocytes'};
pairs{8}={'DLL4','NOTCH3'};cell_pairs{8}={'LSEC','Fibroblasts'};

L=length(pairs);
figure;
for i=1:length(pairs)
    subplot(2,L,i);
    ind_source=find(strcmpi(comb.cell_type,cell_pairs{i}{1}));
    ind_target=find(strcmpi(comb.cell_type,cell_pairs{i}{2}));
    ind_gene=find(strcmpi(comb.gene_name,pairs{i}{1}));
    plot_patch(1:8,smooth(comb.mean_mat{ind_source}(ind_gene,:),3),...
        smooth(comb.se_mat{ind_source}(ind_gene,:),3),'b');
    title([pairs{i}{1} ' - ' cell_pairs{i}{1}],'fontsize',10);
    axis tight;
    set(gca,'xtick',[1 8],'xticklabel',{'Central','Portal'});
    box on;
    set(gca,'fontsize',7)
    %ylim([0 max(ylim)]);
    subplot(2,L,i+L);
    ind_gene=find(strcmpi(comb.gene_name,pairs{i}{2}));
    plot_patch(1:8,smooth(comb.mean_mat{ind_target}(ind_gene,:),3),...
        smooth(comb.se_mat{ind_source}(ind_gene,:),3),'b');
    title([pairs{i}{2} ' - ' cell_pairs{i}{2} ],'fontsize',10);
    box on;
    axis tight;
    set(gca,'xtick',[1 8],'xticklabel',{'Central','Portal'});
    set(gcf,'position',[2         348        1642         526]);
    %ylim([0 max(ylim)]);
    set(gca,'fontsize',7)
end





