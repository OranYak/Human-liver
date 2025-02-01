function t = median_zone_filter_for_github(t,SHOW_PLOTS)
%% assigns the median zone over all non-zero neighbor spots
if nargin<2
   SHOW_PLOTS=0;
end
   
fprintf('calculating pyhsical distances \n');

dists = pdist(t.coor,'euclidean');
step_size = 100/min(dists(dists>0));
t.coor_um = t.coor.*step_size;
D=squareform(pdist(t.coor_um));

THRESH=150;

t.zon_struct.zone_index_med=t.zon_struct.zone_index;
ind_rel_spots=find(t.zon_struct.zone_index>0);
for j=1:length(ind_rel_spots)
    i=ind_rel_spots(j);
    ind_for_med=ind_rel_spots(find(D(i,ind_rel_spots)<THRESH));
    ind_for_med=setdiff(ind_for_med,i);
    % There are spots that are isolated and will get a NaN value in the mean function.
    % Those spots have to stay with their own original zone_index.
    if isempty(ind_for_med)
        continue
    else
        t.zon_struct.zone_index_med(i)=round(median(t.zon_struct.zone_index(ind_for_med)));
    end
end

figure;
scatter(t.coor(:,1),t.coor(:,2),25,repmat(0.7,1,3),'filled'); colorbar;
ind=find(t.zon_struct.zone_index>0);
scatter(t.coor(ind,1),t.coor(ind,2),50,t.zon_struct.zone_index(ind),'filled'); colorbar;
title('non-filtered');
set(gca,'ydir','reverse');
figure;
scatter(t.coor(:,1),t.coor(:,2),25,repmat(0.7,1,3),'filled'); colorbar;
ind=find(t.zon_struct.zone_index_med>0);
scatter(t.coor(ind,1),t.coor(ind,2),50,t.zon_struct.zone_index_med(ind),'filled'); colorbar;
title('median filtered')
set(gca,'ydir','reverse');



    
