T=readtable("supplementary_table_1.xlsx");
clear t;
t.gene_name=T.Var1;
t.samples=T.Properties.VariableNames(2:end);
for i=1:length(t.samples)
    str=strsplit(t.samples{i},'_');
    t.samples{i}=str{2};
end
t.mat_norm=table2array(T(:,2:end));

t.type={'LHD','LHD','LHD','LHD','LHD','LHD','LHD','LHD','Adj','Adj','Adj','Adj','Adj','Adj','Adj','Adj'};
t.sex={'M','F','F','M','M','F','F','F','F','F','F','F','F','M','M','M'};
t.age=[32 45 32 21 27 47 30 57 70 40 57 69 69 42 77 48];
t.histology={'Steatosis','Steatosis','Non steatotic','Non steatotic',...
    'Non steatotic','Non steatotic','Non steatotic','Steatosis',...
    'Steatosis','Steatosis','Steatosis','Steatosis','Severe Steatosis','Fibrosis','Fibrosis','ALPSS'}
for i=1:length(t.samples)
    t.samples{i}=[t.samples{i} '-' t.sex{i}];
end

ind_young=4:6;
ind_old=[2 6];
m_young=mean(t.mat_norm(:,ind_young),2);
m_old=mean(t.mat_norm(:,ind_old),2);
ratio=(1e-5+m_old)./(1e-5+m_young);
mx=max([m_young m_old],[],2);

%% perform DGE of LHD and adjacent while controlling for age

indices=1:12;
exprMatrix=t.mat_norm(:,indices);
exprMatrix=log(1+1e4*exprMatrix);
exprMatrix=(exprMatrix-nanmean(exprMatrix,2))./nanstd(exprMatrix,[],2);


sex=t.sex(indices);
age=t.age(indices);
histology=t.histology(indices);
typ=t.type(indices);

EXP_THRESH=1e-4;
ind_high_genes=find(max(t.mat_norm(:,indices),[],2)>EXP_THRESH);
length(ind_high_genes)
exprMatrix=exprMatrix(ind_high_genes,:);
nGenes = size(exprMatrix, 1);
pvals_Type = zeros(nGenes, 1);
betas_Type = zeros(nGenes, 1);
pvals_Sex = zeros(nGenes, 1);
betas_Sex = zeros(nGenes, 1);
pvals_Histology = zeros(nGenes, 1);
betas_Histology = zeros(nGenes, 1);
pvals_Age = zeros(nGenes, 1);
betas_Age = zeros(nGenes, 1);


for g=1:nGenes
    if mod(g,1000)==0
        g
    end
    y = exprMatrix(g, :)';  % gene expression for gene g
    tbl = table(categorical(typ)', categorical(sex)', age', categorical(histology)', y, 'VariableNames', {'Type','Sex', 'Age', 'Histology','Expr'});

    lm = fitlm(tbl, 'Expr ~ Type + Sex + Age + Histology');  % linear model with covariate
    %lm = fitlm(tbl, 'Expr ~ Type +  Histology');  % linear model with covariate
    coef = lm.Coefficients;
    
    rowIdx = find(contains(coef.Properties.RowNames, 'Type_'));  % Or use the exact name if you know it
    % Extract p-value for Group (second row: intercept, then GroupDisease)
     if ~isempty(rowIdx)
        pvals_Type(g) = coef.pValue(rowIdx);
        betas_Type(g) = coef.Estimate(rowIdx);
    else
        pvals_Type(g) = NaN;
        betas_Type(g) = NaN;
     end

     rowIdx = find(contains(coef.Properties.RowNames, 'Sex_'));  % Or use the exact name if you know it
    % Extract p-value for Group (second row: intercept, then GroupDisease)
     if ~isempty(rowIdx)
        pvals_Sex(g) = coef.pValue(rowIdx);
        betas_Sex(g) = coef.Estimate(rowIdx);
    else
        pvals_Sex(g) = NaN;
        betas_Sex(g) = NaN;
     end

      rowIdx = find(contains(coef.Properties.RowNames, 'Age'));  % Or use the exact name if you know it
    % Extract p-value for Group (second row: intercept, then GroupDisease)
     if ~isempty(rowIdx)
        pvals_Age(g) = coef.pValue(rowIdx);
        betas_Age(g) = coef.Estimate(rowIdx);
    else
        pvals_Age(g) = NaN;
        betas_Age(g) = NaN;
     end

       rowIdx = find(contains(coef.Properties.RowNames, 'Histology_'));  % Or use the exact name if you know it
    % Extract p-value for Group (second row: intercept, then GroupDisease)
     if ~isempty(rowIdx)
        pvals_Histology(g) = min(coef.pValue(rowIdx));
        betas_Histology(g) = max(coef.Estimate(rowIdx));
    else
        pvals_Histology(g) = NaN;
        betas_Histology(g) = NaN;
    end
end

% Multiple testing correction (FDR)
qvals_Type = mafdr(pvals_Type, 'BHFDR', true);
resultTable = table(betas_Type, pvals_Type, qvals_Type, ...
                    'VariableNames', {'Beta_Type', 'Pval_Type', 'FDR_Type'});

% Summarize results
SZ=30;
MXYLIM=5;
PVAL_THRESH=0.01;
figure;
subplot(2,3,1)
scatter(resultTable.Beta_Type,-log10(resultTable.Pval_Type),SZ,repmat(0.6,1,3),'filled')
xlabel('Effect of Type (LHD vs. Adjacent normal)')
ylabel('-log10(pval)');
%ind_signif=find(resultTable.FDR_Type<0.35);
ind_signif=find(resultTable.Pval_Type<PVAL_THRESH);
yline(-log10(PVAL_THRESH),'linestyle','--');
xline(0,'linestyle','--')
ylim([0 MXYLIM]);
A=xlim;
xlim([min(A)-0.5 max(A)+0.5])
hold on;
scatter(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),SZ,'k','filled')
text(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),...
    t.gene_name(ind_high_genes(ind_signif)),'color','k','fontsize',8)
set(gca,'fontsize',16)
box on;
title('Sample type effect','fontsize',24)

% same for histology
% Multiple testing correction (FDR)
qvals_Histology = mafdr(pvals_Histology, 'BHFDR', true);
resultTable = table(betas_Histology, pvals_Histology, qvals_Histology, ...
                    'VariableNames', {'Beta_Type', 'Pval_Type', 'FDR_Type'});

subplot(2,3,2)
scatter(resultTable.Beta_Type,-log10(resultTable.Pval_Type),SZ,repmat(0.6,1,3),'filled')
xlabel('Effect of Histology (steatotic vs. non-steatotic)')
ylabel('-log10(pval)');
%ind_signif=find(resultTable.FDR_Type<0.35);
ind_signif=find(resultTable.Pval_Type<PVAL_THRESH);
yline(-log10(PVAL_THRESH),'linestyle','--');
xline(0,'linestyle','--')
ylim([0 MXYLIM]);
A=xlim;
xlim([min(A)-0.5 max(A)+0.5])
hold on;
scatter(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),SZ,'k','filled')
text(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),...
    t.gene_name(ind_high_genes(ind_signif)),'color','k','fontsize',8)
set(gca,'fontsize',16)
box on;
title('Histology effect','fontsize',24)

% same for Age
% Multiple testing correction (FDR)
qvals_Age = mafdr(pvals_Age, 'BHFDR', true);
resultTable = table(betas_Age, pvals_Age, qvals_Age, ...
                    'VariableNames', {'Beta_Type', 'Pval_Type', 'FDR_Type'});

subplot(2,3,4)
scatter(resultTable.Beta_Type,-log10(resultTable.Pval_Type),SZ,repmat(0.6,1,3),'filled')
xlabel('Effect of Age (old vs. young)')
ylabel('-log10(pval)');
%ind_signif=find(resultTable.FDR_Type<0.35);
ind_signif=find(resultTable.Pval_Type<PVAL_THRESH);
yline(-log10(PVAL_THRESH),'linestyle','--');
xline(0,'linestyle','--')
ylim([0 MXYLIM]);
A=xlim;
xlim([min(A)-0.05 max(A)+0.05])
hold on;
scatter(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),SZ,'k','filled')
text(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),...
    t.gene_name(ind_high_genes(ind_signif)),'color','k','fontsize',8)
set(gca,'fontsize',16)
box on;
title('Age effect','fontsize',24)

% same for Sex
% Multiple testing correction (FDR)
qvals_Sex = mafdr(pvals_Sex, 'BHFDR', true);
resultTable = table(betas_Sex, pvals_Sex, qvals_Sex, ...
                    'VariableNames', {'Beta_Type', 'Pval_Type', 'FDR_Type'});

subplot(2,3,5)
scatter(resultTable.Beta_Type,-log10(resultTable.Pval_Type),SZ,repmat(0.6,1,3),'filled')
xlabel('Effect of Sex (male vs. female)')
ylabel('-log10(pval)');
%ind_signif=find(resultTable.FDR_Type<0.35);
ind_signif=find(resultTable.Pval_Type<PVAL_THRESH);
yline(-log10(PVAL_THRESH),'linestyle','--');
xline(0,'linestyle','--')
ylim([0 2.5*MXYLIM]);
hold on;
scatter(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),SZ,'k','filled')
text(resultTable.Beta_Type(ind_signif),...
    -log10(resultTable.Pval_Type(ind_signif)),...
    t.gene_name(ind_high_genes(ind_signif)),'color','k','fontsize',8)
set(gca,'fontsize',16)
box on;
title('Sex effect','fontsize',24)

% show geometric means of all covariate pvalues
% Clean up p-values: replace 0s and NaNs with small value / remove NaNs
safe_p = @(p) max(p, 1e-300); % avoids log(0)

% Clean each p-value vector
valid_Age = ~isnan(pvals_Age);
valid_Type = ~isnan(pvals_Type);
valid_Sex = ~isnan(pvals_Sex);
valid_Histology = ~isnan(pvals_Histology);

log_age = log(safe_p(pvals_Age(valid_Age)));
log_type = log(safe_p(pvals_Type(valid_Type)));
log_sex = log(safe_p(pvals_Sex(valid_Sex)));
log_histology = log(safe_p(pvals_Histology(valid_Histology)));

% Compute geometric means
geom_age = exp(mean(log_age));
geom_type = exp(mean(log_type));
geom_sex = exp(mean(log_sex));
geom_histology = exp(mean(log_histology));


% Display results
fprintf('Geometric mean p-value:\n');
fprintf('  Age:        %.4e\n', geom_age);
fprintf('  Sex:    %.4e\n', geom_sex);
fprintf('  SampleType: %.4e\n', geom_type);
fprintf('  Histology: %.4e\n', geom_histology);

tab=-log10([ geom_type geom_histology  geom_age geom_sex ]);
subplot(2,3,3);
bar(tab);
set(gca,'xtick',1:length(tab),'xticklabel',{'SampleType','Histology','Age','Sex'})
set(gca,'fontsize',16)
title('-log10(geom-mean-pval)','fontsize',24)
set(gcf,'position',[ 31         487        3258         586]);

subplot(2,3,6);
Ltab=[length(find(pvals_Type<PVAL_THRESH)) length(find(pvals_Histology<PVAL_THRESH)) ...
    length(find(pvals_Age<PVAL_THRESH)) length(find(pvals_Sex<PVAL_THRESH))]
bar(Ltab);
set(gca,'xtick',1:length(Ltab),'xticklabel',{'SampleType','Histology','Age','Sex'})
set(gca,'fontsize',16)
title('# genes with p-value<0.01','fontsize',24)
set(gcf,'position',[ 439         100        1943        1021]);