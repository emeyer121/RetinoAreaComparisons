function [total_surfacearea,retino_surfacearea_smoothwm,retino_surfacearea_pial] = getSubjData( ...
    subjs,species,plot_indiv_subj,allvisual_numbers,nodearea_smoothwm_col,nodearea_pial_col)

for curr_subj = 1:size(subjs,2)
%     cd(num2str(subjs{curr_subj}))
    subj_dir = ['./',species,'/',num2str(subjs{curr_subj}),'/'];

    % File w/ various measures of cortical surface contains 7 columns (use 2 & 3 for area): 
    % 1) node index
    % 2) node area along smoothwm surface
    % 3) node area along pial surface
    % 4) avg triangle area of associated node along smoothwm surface
    % 5) avg triangle area of associated node along pial surface
    % 6) node volume between smoothwm and pial
    % 7) node thickness (Euclidean distance between smoothwm and pial)
    % To calculate total surface area, sum across all nodes
    opt.verb = 0;
    surfmeasures{1}=Read_1D([subj_dir,num2str(subjs{curr_subj}) '_surfmeasures_rh.1D.dset'],opt);
    surfmeasures{2}=Read_1D([subj_dir,num2str(subjs{curr_subj}) '_surfmeasures_lh.1D.dset'],opt);
    
    % ROI files for retinotopic maps w/ 2 columns:
    % 1) node index
    % 2) binary value identifying area assigned to each node (0 = no vis area)
    % see allvisual_numbers and allvisual_labels for identification of each area
    % this is also listed in the areal_labels file
    % IMPORTANT: file has incomplete note list. only has data for nodes in
    % visual areas
    if exist([subj_dir,'rois/allvisual-rh.1D.dset'])
        if strcmp(num2str(species),'Humans')
            allvisual{1}=Read_1D([subj_dir,'rois/allvisual-rh_fovfill.1D.dset'],opt);
            allvisual{2}=Read_1D([subj_dir,'rois/allvisual-lh_fovfill.1D.dset'],opt);
        else
            allvisual{1}=Read_1D([subj_dir,'rois/allvisual-rh.1D.dset'],opt);
            allvisual{2}=Read_1D([subj_dir,'rois/allvisual-lh.1D.dset'],opt);
        end
    else
        fprintf(['No allvisual in ' num2str(subjs{curr_subj}) '\n']);
        norois=[norois,curr_subj];
    end
    if exist([subj_dir,'rois/' num2str(subjs{curr_subj}) '_LGN.1D'])
        temp=load([subj_dir,'rois/' num2str(subjs{curr_subj}) '_LGN.1D']);
        if ~isempty(temp)
            all_LGN{1}(curr_subj)=temp(2);
            all_LGN{2}(curr_subj)=temp(4);
        else
            all_LGN{1}(curr_subj)=nan;
            all_LGN{2}(curr_subj)=nan;
        end
    end
    
    % Surface curvature files w/ 1 columns:
    % 1) surface curvature value 
    % node index = [0:length(curv_rh)-1]
    curv{1}=Read_1D([subj_dir,'rh.curv.1D.dset'],opt);
    curv{2}=Read_1D([subj_dir,'lh.curv.1D.dset'],opt);
    
    % Surface sulcal depth files w/ 1 columns:
    % 1) surface sulcal depth value     
    % node index = [0:length(sulc_rh)-1]
    sulc{1}=Read_1D([subj_dir,'rh.sulc.1D.dset'],opt);
    sulc{2}=Read_1D([subj_dir,'lh.sulc.1D.dset'],opt);

    % Surface thickness files w/ 1 columns:
    % 1) surface thickness value  
    % node index = [0:length(thickness_rh)-1]
    thickness{1}=Read_1D([subj_dir,'rh.thickness.1D.dset'],opt);
    thickness{2}=Read_1D([subj_dir,'lh.thickness.1D.dset'],opt);
    
    %neccessary to create graphs 
    clear h_combined_per_subject;
    clear subject_wm_rh;
    clear subject_wm_lh;
    clear figure
    
    for i=1:2
        V1nodes{i}=allvisual{i}(allvisual{i}(:,2)==1 | allvisual{i}(:,2)==2);
        allothernodes{i}=allvisual{i}(allvisual{i}(:,2)>2);
        
        V1allnodepairs{i}=[];
        for x = 1:ceil(length(V1nodes{i})/10):length(V1nodes{i})
            V1allnodepairs{i}=[V1allnodepairs{i};[repmat(V1nodes{i}(x),[length(allothernodes{i}) 1]),allothernodes{i}]];
        end
    end
    V1allnodepairs_rh = V1allnodepairs{1};
    V1allnodepairs_lh = V1allnodepairs{2};
    save([subj_dir,'V1nodepairs_rh.1D.dset'],'V1allnodepairs_rh','-ascii')
    save([subj_dir,'V1nodepairs_lh.1D.dset'],'V1allnodepairs_lh','-ascii')  
   
    surf_dir=subjs{curr_subj};
    
    
    V1nodepair_distance{1}=Read_1D([subj_dir,'V1_nodedistances_pial_rh.1D'],opt);
    V1nodepair_distance{2}=Read_1D([subj_dir,'V1_nodedistances_pial_lh.1D'],opt);
    for i = 1:2
        V1nodedistances_minmedian{i}=[];
        for x = 1:length(allothernodes{i})
            V1nodedistances_minmedian{i}(x,:)=[allothernodes{i}(x),min(V1nodepair_distance{i}(V1allnodepairs{i}(:,2)==allothernodes{i}(x),3)),median(V1nodepair_distance{i}(V1allnodepairs{i}(:,2)==allothernodes{i}(x),3))];
        end
        V1nodedistances_minmedian_roiavg{i}=[];
        for curr_roi = 3:length(allvisual_numbers)
            curr_avgdist=mean(V1nodedistances_minmedian{i}(ismember(V1nodedistances_minmedian{i}(:,1),allvisual{i}(allvisual{i}(:,2)==allvisual_numbers(curr_roi),1)),2));
            avgdist_allsubjs{i}(curr_roi,curr_subj)=curr_avgdist;
            curr_nodelist=V1nodedistances_minmedian{i}(ismember(V1nodedistances_minmedian{i}(:,1),allvisual{i}(allvisual{i}(:,2)==allvisual_numbers(curr_roi),1)),1);
            V1nodedistances_minmedian_roiavg{i}=[V1nodedistances_minmedian_roiavg{i};[curr_nodelist,repmat(curr_avgdist,[length(curr_nodelist) 1])]];
        end
    end
        
    dlmwrite([subj_dir,'V1nodedistances_minmedian_rh.1D.dset'],V1nodedistances_minmedian{1})
    dlmwrite([subj_dir,'V1nodedistances_minmedian_avgroi_rh.1D.dset'],V1nodedistances_minmedian_roiavg{1})
        
    dlmwrite([subj_dir,'V1nodedistances_minmedian_lh.1D.dset'],V1nodedistances_minmedian{2})
    dlmwrite([subj_dir,'V1nodedistances_minmedian_avgroi_lh.1D.dset'],V1nodedistances_minmedian_roiavg{2})

    for i = 1:2
        V1nodedist{i}{curr_subj}=V1nodedistances_minmedian{i};
        V1nodedist_avgroi{i}{curr_subj}=V1nodedistances_minmedian_roiavg{i};
        if allvisual{i}(~ismember(allvisual{i}(:,2),allvisual_numbers),2)
            allvisual{i}(~ismember(allvisual{i}(:,2),allvisual_numbers),2)
            error([num2str(subjs{curr_subj}) 'has extra ROI values'])
        end
    end    
    
    for i = 1:2
        for curr_roi=1:length(allvisual_numbers)
            curr_area_nodelist{i}=allvisual{i}(allvisual{i}(:,2)==allvisual_numbers(curr_roi),1);
           
            % calculate surface area across entire ROI
            % if code breaks due to empty ROI, add check. 
            retino_surfacearea_smoothwm{i}(curr_roi,curr_subj)=sum(surfmeasures{i}(ismember(surfmeasures{i}(:,1),curr_area_nodelist{i}),nodearea_smoothwm_col));
            retino_surfacearea_pial{i}(curr_roi,curr_subj)=sum(surfmeasures{i}(ismember(surfmeasures{i}(:,1),curr_area_nodelist{i}),nodearea_pial_col));
             
        end
        total_surfacearea{i}(curr_subj) = sum(surfmeasures{i}(:,nodearea_smoothwm_col));
    end
    
    if plot_indiv_subj
        curr_subject_combinedhemi_retino = [retino_surfacearea_smoothwm{1}(:,curr_subj), retino_surfacearea_smoothwm{2}(:,curr_subj)];  % get the last column (curr_subj) and create graph with that 

        %create bar chart comparing R & L hemis surface areas for each subject and save to folder 
        figure(curr_subj);
        bar(curr_subject_combinedhemi_retino, 'grouped');
        hold on
        title('Surface Area of Smooth White Matter, Right and Left Hemispheres');
        xlabel('Visual Areas');
        set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
        ylabel('Surface Area Mean');
        legend("Right Hemisphere", "Left Hemisphere");
%       pngFileName = sprintf('subject_%d.png',curr_subj);
        pngFileName = strcat(subjs(1,curr_subj),'_SurfaceArea.png');
        fullFileName = fullfile(save_folder, pngFileName{1,1});
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .35, .45]);
        title([num2str(subjs{curr_subj})])
        saveas(gcf, fullFileName); 
    end

end
