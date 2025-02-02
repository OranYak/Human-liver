function zon_struct=extract_zonation_for_github(t,EXP_THRESH,NUM_ZONES,output_plots,exclude_pval_computation)

%% extracts zonation patterns of t based on its eta
    
if nargin<2
    EXP_THRESH=0; % pvalues will only be computed for genes above this level of expression
end
if nargin<3
    NUM_ZONES=8;
end
if nargin<4
    output_plots=1;
end
if nargin<5
    exclude_pval_computation=0;
end
zon_struct.NUM_ZONES=NUM_ZONES;
% determine eta percentiles
index_spots=1:length(t.spot_name);
% discarding the fibrotic spots from zonation computaion:
if isfield(t,'ind_fib_spots')
   index_spots=setdiff(index_spots,t.ind_fib_spots);
end
if isfield(t,'ind_capsule_spots')
   index_spots=setdiff(index_spots,t.ind_capsule_spots);
end

thresh=prctile(t.eta(index_spots),(0:NUM_ZONES)*100/NUM_ZONES);
figure;
hist(t.eta,100);
for i=1:length(thresh)
    line([thresh(i) thresh(i)],ylim,'color','r');
end

% assign zones and determine zonation
zon_struct.mn=NaN(length(t.gene_name),NUM_ZONES);
zon_struct.se=NaN(length(t.gene_name),NUM_ZONES);
zon_struct.pval=NaN(length(t.gene_name),1);
zon_struct.qval=NaN(length(t.gene_name),1);
zon_struct.zone_index=zeros(1,length(t.spot_name));

display('computing zonation');
for i=1:length(thresh)-1
    index_spots=find(t.eta>=thresh(i) & t.eta<=thresh(i+1));
    % discarding the fibrotic spots from zonation computaion:
    if isfield(t,'ind_fib_spots')
        index_spots=setdiff(index_spots,t.ind_fib_spots);
    end
    if isfield(t,'ind_capsule_spots')
        index_spots=setdiff(index_spots,t.ind_capsule_spots);
    end
    zon_struct.zone_index(index_spots)=i;
    zon_struct.mn(:,i)=mean(t.mat_norm(:,index_spots),2);%%%%%%%%
    zon_struct.se(:,i)=std(t.mat_norm(:,index_spots),[],2)/sqrt(length(index_spots));
end
        
ind_genes2test=find(max(zon_struct.mn,[],2)>EXP_THRESH);
% add pvalues
display('computing pvalues');
zon_struct.pval=NaN(length(t.gene_name),1);
zon_struct.qval=NaN(length(t.gene_name),1);
if ~exclude_pval_computation
    for j=1:length(ind_genes2test)
        if mod(j,1000)==0
            j
        end
        i=ind_genes2test(j);
        index_spots=1:size(t.mat_norm,2);
        if isfield(t,'ind_fib_spots')
            index_spots=setdiff(index_spots,t.ind_fib_spots);
        end
        if isfield(t,'ind_capsule_spots')
            index_spots=setdiff(index_spots,t.ind_capsule_spots);
        end

        zon_struct.pval(i)=kruskalwallis(t.mat_norm(i,index_spots),zon_struct.zone_index(index_spots),'off');
    end
    zon_struct.qval=mafdr(zon_struct.pval,'bhfdr','true');
end
% compute center of mass
zon_struct.com=zeros(length(t.gene_name),1);
for i=1:length(t.gene_name)
    zon_struct.com(i)=sum((1:NUM_ZONES).*zon_struct.mn(i,:))./sum(zon_struct.mn(i,:));
end

end