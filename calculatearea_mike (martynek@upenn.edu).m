addpath(genpath('afni_matlab'))
addpath(genpath('RSA_functions'))

species='Monks';
computer_name='marcelinamartynek';
%computer_name='marcaro';

plot_indiv_subj=1; % turn on/off single subject graphs

data_dir=['/Users/' num2str(computer_name) '/Box/RetinotopicTopologyProject/' num2str(species) '/'];
save_folder = ['/Users/' num2str(computer_name) '/Box/RetinotopicTopologyProject/Graphs/' num2str(species)]; %folder where graphs are saved

% Labels for allvisual-
% IMPORTANT: some subjects might not have all ROIs
if strcmp(species,'Humans')
    allvisual_numbers=[1,2,3,4,5,6,50,51,52,100,101,102,103,150,151,152,153,154];
    allvisual_labels={'V1v','V1d','V2v','V2d','V3v','V3d','V3A','V3B','V7','LO1','LO2','TO1','TO2','hV4','VO1','VO2','PHC1','PHC2'};
elseif strcmp(species,'Monks')
     allvisual_numbers=[1,2,3,4,5,6,7,8,9,10,50,51,52,53,54,55,56,100,101,102,103,150,151,152,200,201];
     allvisual_labels={'V1v','V1d','V2v','V2d','V3v','V3d','V4v','V4d','V4Av','V4Ad','V3A','DP','CIP1','CIP2','LIP1','LIP2','LIP3','MT','MST','FST','V4t','OTd','PITv','PITd','OTS1','OTS2'};       
else
   fprintf('Error. Unknown species.\n')
end

nodearea_smoothwm_col=2;
nodearea_pial_col=3;
triarea_smoothwm_col=4;
triarea_pial_col=5;
nodevol_col=6;
nodethickness_col=7;

% This should be a more efficient way of filtering the directories:
files = dir(data_dir);
names = {files.name}; 
dirFlags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..') & ~strcmp(names, 'phasemaps') & ~strcmp(names, 'VTPM') & ~strcmp(names, 'VTPM_delayedsaccade') & ~strcmp(names, 'ignore');
subjs=names(dirFlags);
norois=[];

cd(data_dir)

for curr_subj = 1:size(subjs,2)
    cd(num2str(subjs{curr_subj}))
    
    % File w/ various measures of cortical surface contains 7 columns (use 2 & 3 for area): 
    % 1) node index
    % 2) node area along smoothwm surface
    % 3) node area along pial surface
    % 4) avg triangle area of associated node along smoothwm surface
    % 5) avg triangle area of associated node along pial surface
    % 6) node volume between smoothwm and pial
    % 7) node thickness (Euclidean distance between smoothwm and pial)
    % To calculate total surface area, sum across all nodes
    surfmeasures_rh=Read_1D([num2str(subjs{curr_subj}) '_surfmeasures_rh.1D.dset']);
    surfmeasures_lh=Read_1D([num2str(subjs{curr_subj}) '_surfmeasures_lh.1D.dset']);
    
    % ROI files for retinotopic maps w/ 2 columns:
    % 1) node index
    % 2) binary value identifying area assigned to each node (0 = no vis area)
    % see allvisual_numbers and allvisual_labels for identification of each area
    % this is also listed in the areal_labels file
    % IMPORTANT: file has incomplete note list. only has data for nodes in
    % visual areas
    if exist(['rois/allvisual-rh.1D.dset'])
        allvisual_rh=Read_1D(['rois/allvisual-rh.1D.dset']);
        allvisual_lh=Read_1D(['rois/allvisual-lh.1D.dset']);
    else
        fprintf(['No allvisual in ' num2str(subjs{curr_subj}) '\n']);
        norois=[norois,curr_subj];
    end
    
    % Surface curvature files w/ 1 columns:
    % 1) surface curvature value 
    % node index = [0:length(curv_rh)-1]
    curv_rh=Read_1D(['rh.curv.1D.dset']);
    curv_lh=Read_1D(['lh.curv.1D.dset']);
    
    % Surface sulcal depth files w/ 1 columns:
    % 1) surface sulcal depth value     
    % node index = [0:length(sulc_rh)-1]
    sulc_rh=Read_1D(['rh.sulc.1D.dset']);
    sulc_lh=Read_1D(['lh.sulc.1D.dset']);

    % Surface thickness files w/ 1 columns:
    % 1) surface thickness value  
    % node index = [0:length(thickness_rh)-1]
    thickness_rh=Read_1D(['rh.thickness.1D.dset']);
    thickness_lh=Read_1D(['lh.thickness.1D.dset']);
    
    %neccessary to create graphs 
    clear h_combined_per_subject;
    clear subject_wm_rh;
    clear subject_wm_lh;
    clear figure
    
    for curr_roi=1:length(allvisual_numbers)
        curr_area_nodelist_rh=allvisual_rh(allvisual_rh(:,2)==allvisual_numbers(curr_roi),1);
        curr_area_nodelist_lh=allvisual_lh(allvisual_lh(:,2)==allvisual_numbers(curr_roi),1);
       
        % calculate surface area across entire ROI
        % if code breaks due to empty ROI, add check. 
        retino_surfacearea_smoothwm_rh(curr_roi,curr_subj)=sum(surfmeasures_rh(ismember(surfmeasures_rh(:,1),curr_area_nodelist_rh),nodearea_smoothwm_col));
        retino_surfacearea_pial_rh(curr_roi,curr_subj)=sum(surfmeasures_rh(ismember(surfmeasures_rh(:,1),curr_area_nodelist_rh),nodearea_pial_col));
              
        retino_surfacearea_smoothwm_lh(curr_roi,curr_subj)=sum(surfmeasures_lh(ismember(surfmeasures_lh(:,1),curr_area_nodelist_lh),nodearea_smoothwm_col));
        retino_surfacearea_pial_lh(curr_roi,curr_subj)=sum(surfmeasures_lh(ismember(surfmeasures_lh(:,1),curr_area_nodelist_lh),nodearea_pial_col));
        
        total_surfacearea_rh(:,curr_subj) = sum(surfmeasures_rh(:,nodearea_smoothwm_col)); 
        total_surfacearea_lh(:,curr_subj) = sum(surfmeasures_lh(:,nodearea_smoothwm_col));    
        
    end
    
    if plot_indiv_subj
        curr_subject_combinedhemi_retino = [retino_surfacearea_smoothwm_rh(:,curr_subj), retino_surfacearea_smoothwm_lh(:,curr_subj)];  % get the last column (curr_subj) and create graph with that 

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
        title([num2str(subjs{curr_subj})])
        saveas(gcf, fullFileName); 
    end

    cd ../
end

total_surfacearea_bothhemis=total_surfacearea_rh+total_surfacearea_lh;

% Make 0 sum ROIs - NaN
retino_surfacearea_smoothwm_rh(retino_surfacearea_smoothwm_rh==0)=nan;
retino_surfacearea_smoothwm_lh(retino_surfacearea_smoothwm_lh==0)=nan;
retino_surfacearea_pial_rh(retino_surfacearea_pial_rh==0)=nan;
retino_surfacearea_pial_lh(retino_surfacearea_pial_lh==0)=nan;

retino_surfacearea_smoothwm_bothhemis = retino_surfacearea_smoothwm_rh + retino_surfacearea_smoothwm_lh; %create matrix of added surface ares for left and right hemispheres 
mean_retinoareas_bothhemis = nanmean(retino_surfacearea_smoothwm_bothhemis,2); %take average for each visual area across subject 
stderror_retinoareas_bothhemis = nanstd(retino_surfacearea_smoothwm_bothhemis,0,2) / sqrt(size(retino_surfacearea_smoothwm_bothhemis,1)); %find SE of surface area size for visual areas

retino_surfacearea_smoothwm_catbothhemis=cat(2,retino_surfacearea_smoothwm_rh,retino_surfacearea_smoothwm_lh);

%% Comparisons of surface area between hemispheres
crosshemi_corrmatrix=corr(retino_surfacearea_smoothwm_rh',retino_surfacearea_smoothwm_lh','rows','complete');
figure
subplot(1,2,1)
plot(diag(crosshemi_corrmatrix))
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation (r)');
title('Across subject Correlation of Surface Area between Hemispheres')
subplot(1,2,2)
plot(nanmean(retino_surfacearea_smoothwm_rh,2),nanmean(retino_surfacearea_smoothwm_lh,2),'.k','MarkerSize',15);
%lsline
hold on
max_value=ceil(max(max([nanmean(retino_surfacearea_smoothwm_rh,2),nanmean(retino_surfacearea_smoothwm_lh,2)])));
min_value=floor(min(min([nanmean(retino_surfacearea_smoothwm_rh,2),nanmean(retino_surfacearea_smoothwm_lh,2)])));
plot([min_value, max_value],[min_value,max_value],'--k','LineWidth',2)
xlabel('Mean Surface Area RH');
ylabel('Mean Surface Area LH');
title('Comparison of Mean Surface Area between Hemispheres')
set(gcf,'Units','Normalized','OuterPosition',[0, 10, .6, .5])
pngFileName = ('CrossHemisphereComparisons.png');
fullFileName = fullfile(save_folder, pngFileName);
saveas(gcf, fullFileName); 


%% Bar graph showing across subject mean + standard error of surface area size for each visual area
figure
avgmeanbar = bar([1:length(allvisual_labels)],mean_retinoareas_bothhemis);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Surface Area Mean');
title('Across Subject Mean')

hold on
er = errorbar([1:length(allvisual_labels)],mean_retinoareas_bothhemis,stderror_retinoareas_bothhemis);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
pngFileName = ('AcrossSubjectMean.png');
fullFileName = fullfile(save_folder, pngFileName);
saveas(gcf, fullFileName); 

hold off

%% normalizing + plotting 

% For now, focus on nodearea smoothwm (add pial if you have time)
% run corr 
% to get around Nan: corr( X, Y, 'rows','complete')
% eg., corr( array of V3A sizes across subjects, array of V1 sizes)

%corr between all regions 

% hvae to create loop that takes each ROI from each subject --> nested for
% loops here --> okay 
%[rho,pval] = corr(V1_surfacearea, retino_surfacearea_smoothwm_bothhemis, 'complete');
%imagesc([r,l])

%normalize by V1 surface area 
%normalize by V1v surface area 
%normalize by V1d surface area

V1_surfacearea = (retino_surfacearea_smoothwm_bothhemis(1,:) + retino_surfacearea_smoothwm_bothhemis(2,:)); %full V1 of combined hemispheres 
V2_surfacearea = (retino_surfacearea_smoothwm_bothhemis(3,:) + retino_surfacearea_smoothwm_bothhemis(4,:)); %full V1 of combined hemispheres 
V3_surfacearea = (retino_surfacearea_smoothwm_bothhemis(5,:) + retino_surfacearea_smoothwm_bothhemis(6,:)); %full V1 of combined hemispheres 
V1v_surfacearea = retino_surfacearea_smoothwm_bothhemis(1,:);
V1d_surfacearea = retino_surfacearea_smoothwm_bothhemis(2,:);
V3d_surfacearea = retino_surfacearea_smoothwm_bothhemis(6,:);
V3v_surfacearea = retino_surfacearea_smoothwm_bothhemis(5,:);
V3A_surfacearea = retino_surfacearea_smoothwm_bothhemis(11,:); %7 if Humans

 
%human areas
hV4_surfacearea = retino_surfacearea_smoothwm_bothhemis(14,:);
TO1_surfacearea = retino_surfacearea_smoothwm_bothhemis(12,:); 

N = size(retino_surfacearea_smoothwm_bothhemis, 1);
%corr_V1 = zeros(1,N);

corr_V1 = corr(V1_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1= bar([1:length(allvisual_labels)],corr_V1);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V1 Surface Area')

corr_V1v =  corr(V1v_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1v= bar([1:length(allvisual_labels)],corr_V1v);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V1v Surface Area')


corr_V1d = corr(V1d_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1d= bar([1:length(allvisual_labels)],corr_V1d);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V1d Surface Area')


corr_V3d = corr(V3d_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1d= bar([1:length(allvisual_labels)],corr_V3d);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V3d Surface Area')

corr_V3v = corr(V3v_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1d= bar([1:length(allvisual_labels)],corr_V3v);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V3v Surface Area')

corr_hV4 = corr(hV4_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
barV1d= bar([1:length(allvisual_labels)],corr_hV4);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by hV4 Surface Area')

corr_V3A = corr(V3A_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
figure
barV1d= bar([1:length(allvisual_labels)],corr_V3A);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by V3A Surface Area')

corr_TO1 = corr(TO1_surfacearea', retino_surfacearea_smoothwm_bothhemis', 'rows','complete')
barV1d= bar([1:length(allvisual_labels)],corr_TO1);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Normalization by TO1 Surface Area')


%% Partial Correlations 


%partial corr V3d: right and left hemisphere surface area
pcorr_V3d = partialcorr(V3d_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV3d = bar([1:length(allvisual_labels)],pcorr_V3d);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')


%partial corr V3v: right and left hemisphere surface area
pcorr_V3v = partialcorr(V3v_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV3v = bar([1:length(allvisual_labels)],pcorr_V3v);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')


%partial corr hV4: right and left hemisphere surface area
pcorr_hV4 = partialcorr(hV4_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barhV4 = bar([1:length(allvisual_labels)],pcorr_hV4);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation by Right Hemisphere Surface Area')


%partial corr V3A: right and left hemisphere surface area
pcorr_V3A = partialcorr(V3A_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV3A = bar([1:length(allvisual_labels)],pcorr_V3A);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')



%partial corr TO1: right and left hemisphere surface area
pcorr_TO1 = partialcorr(TO1_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barTO1 = bar([1:length(allvisual_labels)],pcorr_TO1);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')


%partial corr V1: right and left hemisphere surface area
pcorr_V1 = partialcorr(V1_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV1 = bar([1:length(allvisual_labels)],pcorr_V1);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')


%partial corr V1v: right and left hemisphere surface area
pcorr_V1v = partialcorr(V1v_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV1v = bar([1:length(allvisual_labels)],pcorr_V1v);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')

%partial corr V1d: right and left hemisphere surface area
pcorr_V1d = partialcorr(V1d_surfacearea', retino_surfacearea_smoothwm_bothhemis', total_surfacearea_bothhemis', 'rows','complete')
barV1d = bar([1:length(allvisual_labels)],pcorr_V1d);
set(gca,'xtick',[1:length(allvisual_labels)],'xticklabel', allvisual_labels);
xlabel('Visual Areas');
ylabel('Correlation');
title('Partial Correlation Controlling for Total Hemisphere Surface Area')


% RSA 
retino_surfacearea_simmatrix=corr(retino_surfacearea_smoothwm_bothhemis','rows','complete');
plotModularity(retino_surfacearea_simmatrix,allvisual_labels,'newmansq');

retino_surfacearea_smoothwm_bothhemis_noquadrants=retino_surfacearea_smoothwm_bothhemis;
retino_surfacearea_smoothwm_bothhemis_noquadrants(1,:)=retino_surfacearea_smoothwm_bothhemis_noquadrants(1,:)+retino_surfacearea_smoothwm_bothhemis_noquadrants(2,:);
retino_surfacearea_smoothwm_bothhemis_noquadrants(2,:)=retino_surfacearea_smoothwm_bothhemis_noquadrants(3,:)+retino_surfacearea_smoothwm_bothhemis_noquadrants(4,:);
retino_surfacearea_smoothwm_bothhemis_noquadrants(3,:)=retino_surfacearea_smoothwm_bothhemis_noquadrants(5,:)+retino_surfacearea_smoothwm_bothhemis_noquadrants(6,:);
retino_surfacearea_smoothwm_bothhemis_noquadrants(4,:)=[];
retino_surfacearea_smoothwm_bothhemis_noquadrants(5,:)=[];
retino_surfacearea_smoothwm_bothhemis_noquadrants(6,:)=[];
allvisual_labels_noquadrants={'V1','V2','V3','V3A','V3B','V7','LO1','LO2','TO1','TO2','hV4','VO1','VO2','PHC1','PHC2'};
retino_surfacearea_noquadrants_simmatrix=corr(retino_surfacearea_smoothwm_bothhemis_noquadrants','rows','complete');
plotModularity(retino_surfacearea_noquadrants_simmatrix,allvisual_labels_noquadrants,'newmansq');


% eg., partialcorr( array of V3A sizes across subjects, array of V1 sizes, array total surface area)
% across subjects)

% create correlation matrix: area x area correlations
% create partial correlation matrix: area x area correlations removing
% total area
% create partial correlation matrix: area x area correlations removing
% V1 area
