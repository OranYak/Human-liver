load v
load zon_struct_all_full

%% Choose normalization for mat_norm (normalization over all genes):
for i=1:length(v)
    i
    v{i}.mat_norm=v{i}.mat_norm_all_genes; 
end

% Choose zon_struct to work on:
for i=1:length(v)
    v{i}.zon_struct=v{i}.zon_struct_not_hep_norm; 
end

%% clustergram- figure 1c:
EXP_THRESH=5e-4;
clear ttl;
zone=1:8;
analyzed_samples=length(v);
analyzed_samples=16;  

for i=1:analyzed_samples
    ttl{i}=v{i}.main_feature;
end

% find the common genes:
common_genes=v{1}.gene_name;
if analyzed_samples>16
    for i=2:length(v)
        common_genes=intersect(common_genes,v{i}.gene_name);
    end
end

mat_all=zeros(length(v{1}.gene_name),analyzed_samples);

for i=1:analyzed_samples
    indin=1:length(v{i}.spot_name);
    if isfield(v{i},'ind_capsule_spots')
        indin=setdiff(indin,v{i}.ind_capsule_spots);
    end
    if isfield(v{i},'ind_fib_spots')
        indin=setdiff(indin,v{i}.ind_fib_spots);
    end

    [Ia,Ib]=ismember(common_genes,v{i}.gene_name);
    ind_common_genes=Ib(Ia);
 
    mat_all(:,i)=mean(v{i}.mat_norm(ind_common_genes,indin),2);
end

% clustergram
ind_genes=find(max(mat_all,[],2)>EXP_THRESH);

ind_genes=setdiff(ind_genes,find(startsWith(v{1}.gene_name,'RPL')));
ind_genes=setdiff(ind_genes,find(startsWith(v{1}.gene_name,'RPS')));
ind_genes=setdiff(ind_genes,find(startsWith(v{1}.gene_name,'MT-')));

mat=mat_all./sum(mat_all(ind_genes,:));
mat=log10(1e-6+mat);

Zmat=(mat-mean(mat,2))./std(mat,[],2);
rowLabels = v{1}.gene_name(ind_genes);
cp=clustergram(Zmat(ind_genes,:),'rowpdist','spearman','columnpdist','spearman','rowlabels',rowLabels,'columnlabels',ttl,'colormap','redbluecmap');
cm = struct('GroupNumber',{12,14},'Annotation',{'Time1','Time2'},'Color',{'r','k'});
set(cp,'ColumnGroupMarker',cm);
set(cp,'ColumnLabelsRotate', 30);

% In the matlab clustergram function: If the number of row labels is 200 or more, the labels do not appear in the clustergram plot.
rowLabels = v{1}.gene_name(ind_genes);
showed_gns={'krt','timp1','s100a8','s100a9','pck1','cps','alb','cyp27a1'};
showed_gns=upper(showed_gns);
[ia,ib]=ismember(showed_gns',v{i}.gene_name(ind_genes));
labelIndices=ib(ia);
v{i}.gene_name(ind_genes(labelIndices))
rowLabels = repmat({''}, size(v{1}.gene_name(ind_genes)));
rowLabels(labelIndices) = v{i}.gene_name(ind_genes(labelIndices));


%% PCA- figure 1d:
SZ=150;
[coeff,score,latent]=pca(Zmat(ind_genes,:)');
sample_indicator=zeros(1,length(ttl));
sample_indicator(find(contains(ttl,'fib')))=2;
sample_indicator(find(contains(ttl,'cap')))=2;
sample_indicator(find(contains(ttl,'steato')))=1;
sample_indicator(find(contains(ttl,'ALPPS')))=3;
if analyzed_samples>16
    sample_indicator(find(contains(ttl,'Guilliams')))=4;
    sample_indicator(find(contains(ttl,'Macparland')))=4;
end
ind_LHD=find(startsWith(ttl,'M'));

figure;
scatter(score(:,1),score(:,2),SZ,sample_indicator,'filled');
hold on
scatter(score(:,1),score(:,2),SZ,'k','LineWidth',2);
scatter(score(ind_LHD,1),score(ind_LHD,2),SZ,'r','LineWidth',2);
hold on;
text(score(:,1),score(:,2),ttl);
hold on;
line([0.7 -0.7], ylim, 'LineStyle', '--', 'Color', 'k');

hold on
xlabel(['PC1 (' num2str(100*latent(1)/sum(latent),'%.1f') '%)']);
ylabel(['PC2 (' num2str(100*latent(2)/sum(latent),'%.1f') '%)']);
box on;


%% Figure 1e

EXP_THRESH=5e-5;
QTHRESH=0.25;
RATIO_THRESH=1.5;
ind_LHD=1:8;
ind_adjacent=9:16;
gns=v{1}.gene_name;
ind1=ind_LHD;
ind2=ind_adjacent;
m1=median(mat_all(:,ind1),2);
m2=median(mat_all(:,ind2),2);

mx=max([m1 m2],[],2);
vec=[m1 m2];
vec=vec(:);
pn=min(vec(vec>0));
ratio=(pn+m1)./(pn+m2);

%add pval and qval
pval=NaN(length(gns),1);
qval=NaN(length(gns),1);
indin=find(mx>EXP_THRESH);
for i=1:length(indin)
    [~,pval(indin(i))]=ttest2(mat_all(indin(i),ind1),mat_all(indin(i),ind2));
end
qval(indin)=mafdr(pval(indin),'bhfdr','true');
ind_signif=find(qval<QTHRESH & (ratio>RATIO_THRESH | ratio<1/RATIO_THRESH));

figure;
scatter(log10(mx),log2(ratio),'.');
hold on;
scatter(log10(mx(ind_signif)),log2(ratio(ind_signif)),'r.');
xlabel('Max of means of healthy and Adjacent normal');
ylabel('log2(Healthy/Adjacent normal)');
xlim([-4.5 -0.5])   
ylim([-4.8 2.5])
hold on
box on;

inter_genes_nams={'SDHB','C1QA','C1QB','S100A8','S100A9','S100A10','S100A11','S100A13','IL6R','CRP','SLC5A6','IGKC','IGFBP5','SRPRB','TNFSF10','ALB','HMGCS1','HLA-C','HLA-B','GSTA2','GSTA1','COX7A2','SOD2','CXCL12','CYP2C19','CYP2C8','OAT','CYP2E1','SAA4','SAA2','SAA1','LDHA','SDHD','NNMT','GAPDH','CD4','KRT8','KRT18','GLS2','DCN','PCK2','IGHA1','IGHG1','IGHG3','TNFRSF12A','IL32','IL4R','LDLR','CYP4F3','IGLC1','IGLC2','IGLC3','PPARA','GK'};
[ia,ib]=ismember(inter_genes_nams',gns);
inter_genes=ib(ia);

gns(inter_genes)
hold on;
scatter(log10(mx(inter_genes)),log2(ratio(inter_genes)),'o','m');
text(log10(mx(inter_genes)),log2(ratio(inter_genes)),gns(inter_genes));
ylim([-5.5 2.7])

