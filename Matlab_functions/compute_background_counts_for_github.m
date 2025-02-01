function bg_vec = compute_background_counts_for_github(folder,patient,METHOD)

% this function computes the mean counts for all the spots that are not
% under the tissue and brings back this vector

% load the count table with all spots
tab = readtable([folder,patient,'_counts_ALL.csv']);
gns = tab.Var1;

spots_AL = tab.Properties.VariableNames(2:end)';
count_AL = table2array(tab(:,2:end)); 
[~,ord] = sort(spots_AL);
spots_AL = spots_AL(ord);
count_AL = count_AL(:,ord);

% find the spots that were detected under the tissue with space ranger
fid = fopen([folder,patient,'_counts_UTT.csv']);
a = textscan(fid,'%s',1);
fclose(fid);
A = split(a{1},';');
A = A(2:end);
spots_UT_ind = zeros(length(A),1);
for i=1:length(A)
        str=A{i};
        str(findstr(str,'-'))='_';
        A{i}=replace(str,'"','');
        ind = find(strcmpi(spots_AL,A{i}));
        if ~isempty(ind)
            spots_UT_ind(i) = ind;
        end
end
spots_UT = A;

% load the coordiate to plot
xy = readtable([folder,patient,'_tissue_positions_list.csv']);

indin = zeros(size(xy,1),1);
for i=1:size(xy,1)
    str=xy.Var1{i};
    str(findstr(str,'-'))='_';
    xy.Var1{i} = str;
end
[~,ord] = sort(xy.Var1);
xy = xy(ord,:);

[~,spots_NUT_ind] = setdiff(spots_AL,spots_UT);
if strcmpi(METHOD,'mean')
    bg_vec = mean(count_AL(:,spots_NUT_ind),2);
elseif strcmpi(METHOD,'median')
    bg_vec = median(count_AL(:,spots_NUT_ind),2);
else
    error('ENTER A METHOD TO USE FOR THE BG VECTOR COMPUTATION');
end

figure;
nexttile;
scatter(xy.Var6,xy.Var5,55,log10(sum(count_AL,1)),"filled");
colorbar;
hold on;
scatter(xy.Var6(spots_NUT_ind),xy.Var5(spots_NUT_ind),15,'.r');
title(patient);

[y,ord] = sort(bg_vec);
gns(ord);
nexttile;
semilogy(y);
