function modularity_group= plotModularity(corrMatrix,labels,clusterApproach)

if strfind(clusterApproach,'kmeans')
    if isempty(strrep(clusterApproach,'kmeans',''))  
        %if no cluster size specified, run kmeans up to 20 and choose
        %optimal cluster size
        figure
        hold on
        for x = 1:20
        modularity_group=kmeans(corrMatrix,x,'Distance','correlation','MaxIter',100,'Replicates',10);
        subplot(2,1,1)
        [silh,h] = silhouette(corrMatrix,modularity_group,'correlation');
        silh_all(:,x)=silh;
        end
        mean(silh_all)
        subplot(2,1,2)
        plot(mean(silh_all))
        cluster_size=input('Select cluster size\n');
    else
        cluster_size=str2num(strrep(clusterApproach,'kmeans',''));
    end
    modularity_group=kmeans(corrMatrix,cluster_size,'Distance','correlation','MaxIter',100,'Replicates',10);
elseif strfind(clusterApproach,'louvain')
    modularity_group=community_louvain(corrMatrix,[],[],'negative_sym');
elseif strfind(clusterApproach,'newmansq')
    [modularity_group, Q]=modularity_und(corrMatrix);
else
    error('Could not recognize cluster approach')
end

use_rois=[1:size(labels,1)];
mod_color={'b','r','g'};
color_range=jet;
thresh=[];

cortical_corr=corrMatrix;
%[Y,eigvals]=cmdscale(1-cortical_corr);
[Y,eigvals]=cmdscale(pdist(cortical_corr, 'euclidean'));
format short g
[eigvals eigvals/max(abs(eigvals))];
% figure 
% hold on
% for x = 1: size(corrMatrix,1)
%     for y = 1:size(corrMatrix,1)
%         if cortical_corr(x,y) > 0
%             curr_color=ceil(cortical_corr(x,y)*size(color_range,1));
%             plot([Y(x,1),Y(y,1)],[Y(x,2),Y(y,2)],'-','color',[color_range(curr_color,:)],'LineWidth',(cortical_corr(x,y)*4)+2);
%         end
%     end
% end

mod_colors=distinguishable_colors(max(modularity_group));
figure; 
hold on
offset=(max(Y(:))-min(Y(:)))/100;
for x = 1: size(corrMatrix,1)
    plot([Y(x,1),Y(x,1)],[Y(x,2),Y(x,2)],'.k','MarkerSize',45);
    plot([Y(x,1),Y(x,1)],[Y(x,2),Y(x,2)],'.c','MarkerSize',30,'color',mod_colors(modularity_group(x),:));
    if iscell(labels)
        text([Y(x,1),Y(x,1)]+offset,[Y(x,2),Y(x,2)]+offset,{labels{x}},'FontSize',18,'color','k');
    else
        text([Y(x,1),Y(x,1)]+offset,[Y(x,2),Y(x,2)]+offset,{labels(x)},'FontSize',18,'color','k');   
    end
end

