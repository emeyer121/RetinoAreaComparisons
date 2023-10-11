addpath(genpath('/Applications/MATLAB_TOOLBOXES/afni_matlab'))
folder_path = '/Users/marcaro/Library/CloudStorage/Box-Box/RetinotopicTopologyProject/Humans';  % replace 'X' with the path to the directory you want to search in
mris_convert_path='/Applications/freesurfer/bin/mris_convert';

cd([folder_path])
    
files = dir(folder_path);
folders = {files([files.isdir]).name}; % get only folders

% use regular expressions to filter folders with length of lowercase 4 letters
regexp_filter = '^[a-z]{4}$'; % regular expression to match lowercase 4-letter strings
subj_folders = folders(cellfun(@(x) ~isempty(x), regexp(folders, regexp_filter)));

% display the list of folders
disp(subj_folders);

roi_value=[1,2,3,4,5,6,50,51,100,101,102,103,150,151,152,153,154];
parietal_values=[1,2,3,4,5,6,7,8,9]+200;
allrois=[0,roi_value,parietal_values];
allsubjs_alluniquerois_rh=zeros(length(allrois),size(subj_folders,2));
allsubjs_alluniquerois_lh=zeros(length(allrois),size(subj_folders,2));
alluniqueroinodepairs=allcomb(allrois,allrois);

allareadistancespial_rh=zeros(size(alluniqueroinodepairs,1),size(subj_folders,2));
allareadistancessmoothwm_rh=zeros(size(alluniqueroinodepairs,1),size(subj_folders,2));
allareadistancespial_lh=zeros(size(alluniqueroinodepairs,1),size(subj_folders,2));
allareadistancessmoothwm_lh=zeros(size(alluniqueroinodepairs,1),size(subj_folders,2));

for curr_subj = 1:size(subj_folders,2)  
    curr_surfdir=[subj_folders{curr_subj} '/surfaces/'];
    curr_roidir=[subj_folders{curr_subj} '/rois/'];

    clear datalines
    % Specify the filename
    filename = [curr_roidir 'uniquenodepairs_roicenters_nodedistances_pial_rh.1D']; % Replace 'your_file_name.txt' with the actual filename
    % Open the file for reading
    fileID = fopen(filename, 'r');
    % Read lines from the file and ignore lines starting with "#"
    dataLines = textscan(fileID, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
    fclose(fileID);
    % Remove the lines that start with "#" and the line with column names
    dataLines = dataLines{1}(3:end);
    % Preallocate a matrix to store the numeric values
    matrix = zeros(numel(dataLines), 3);
    % Extract and store the numeric values in the matrix
    for i = 1:numel(dataLines)
        values = sscanf(dataLines{i}, '%f');
        matrix(i, :) = values';
    end
    
    currentpairs=Read_1D([curr_roidir 'uniquenodepairs_roicenters_rh.1D']);
    allsubjs_alluniquerois_rh(ismember(allrois,unique(currentpairs(:,1))),curr_subj)=1;
    
    [rowsInMatrix1, locInMatrix2] = ismember(alluniqueroinodepairs, currentpairs, 'rows');
    allareadistancespial_rh(rowsInMatrix1,curr_subj)=matrix(:,3);
    if sum(sum(alluniqueroinodepairs(rowsInMatrix1,1:2)-currentpairs))>0
        warning(['There is a mismatch between the roi pairs in uniquenodepairs_roicenters and this in allpossible unique for subj' num2str(curr_subj) ' in RH'])
    end
    
    clear datalines
    % Specify the filename
    filename = [curr_roidir 'uniquenodepairs_roicenters_nodedistances_smoothwm_rh.1D']; % Replace 'your_file_name.txt' with the actual filename
    % Open the file for reading
    fileID = fopen(filename, 'r');
    % Read lines from the file and ignore lines starting with "#"
    dataLines = textscan(fileID, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
    fclose(fileID);
    % Remove the lines that start with "#" and the line with column names
    dataLines = dataLines{1}(3:end);
    % Preallocate a matrix to store the numeric values
    matrix = zeros(numel(dataLines), 3);
    % Extract and store the numeric values in the matrix
    for i = 1:numel(dataLines)
        values = sscanf(dataLines{i}, '%f');
        matrix(i, :) = values';
    end

    [rowsInMatrix1, locInMatrix2] = ismember(alluniqueroinodepairs, currentpairs, 'rows');
    allareadistancessmoothwm_rh(rowsInMatrix1,curr_subj)=matrix(:,3);
    if sum(sum(alluniqueroinodepairs(rowsInMatrix1,1:2)-currentpairs))>0
        warning(['There is a mismatch between the roi pairs in uniquenodepairs_roicenters and this in allpossible unique for subj' num2str(curr_subj) ' in RH'])
    end
    
    
     clear datalines
    % Specify the filename
    filename = [curr_roidir 'uniquenodepairs_roicenters_nodedistances_pial_lh.1D']; % Replace 'your_file_name.txt' with the actual filename
    % Open the file for reading
    fileID = fopen(filename, 'r');
    % Read lines from the file and ignore lines starting with "#"
    dataLines = textscan(fileID, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
    fclose(fileID);
    % Remove the lines that start with "#" and the line with column names
    dataLines = dataLines{1}(3:end);
    % Preallocate a matrix to store the numeric values
    matrix = zeros(numel(dataLines), 3);
    % Extract and store the numeric values in the matrix
    for i = 1:numel(dataLines)
        values = sscanf(dataLines{i}, '%f');
        matrix(i, :) = values';
    end

    currentpairs=Read_1D([curr_roidir 'uniquenodepairs_roicenters_lh.1D']);
    allsubjs_alluniquerois_lh(ismember(allrois,unique(currentpairs(:,1))),curr_subj)=1;
    
    [rowsInMatrix1, locInMatrix2] = ismember(alluniqueroinodepairs, currentpairs, 'rows');
    allareadistancespial_lh(rowsInMatrix1,curr_subj)=matrix(:,3);
    if sum(sum(alluniqueroinodepairs(rowsInMatrix1,1:2)-currentpairs))>0
        warning(['There is a mismatch between the roi pairs in uniquenodepairs_roicenters and this in allpossible unique for subj' num2str(curr_subj) ' in LH'])
    end    
    
    clear datalines
    % Specify the filename
    filename = [curr_roidir 'uniquenodepairs_roicenters_nodedistances_smoothwm_lh.1D']; % Replace 'your_file_name.txt' with the actual filename
    % Open the file for reading
    fileID = fopen(filename, 'r');
    % Read lines from the file and ignore lines starting with "#"
    dataLines = textscan(fileID, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
    fclose(fileID);
    % Remove the lines that start with "#" and the line with column names
    dataLines = dataLines{1}(3:end);
    % Preallocate a matrix to store the numeric values
    matrix = zeros(numel(dataLines), 3);
    % Extract and store the numeric values in the matrix
    for i = 1:numel(dataLines)
        values = sscanf(dataLines{i}, '%f');
        matrix(i, :) = values';
    end

    [rowsInMatrix1, locInMatrix2] = ismember(alluniqueroinodepairs, currentpairs, 'rows');
    allareadistancessmoothwm_lh(rowsInMatrix1,curr_subj)=matrix(:,3);
    if sum(sum(alluniqueroinodepairs(rowsInMatrix1,1:2)-currentpairs))>0
        warning(['There is a mismatch between the roi pairs in uniquenodepairs_roicenters and this in allpossible unique for subj' num2str(curr_subj) ' in LH'])
    end
    
end

allareadistancesavglayer_rh=[allareadistancessmoothwm_rh+allareadistancespial_rh]/2;
allareadistancesavglayer_lh=[allareadistancessmoothwm_lh+allareadistancespial_lh]/2;

allareadistancesavglayer_rh(allareadistancesavglayer_rh==0)=nan;
allareadistancesavglayer_lh(allareadistancesavglayer_lh==0)=nan;

area_labels={'V1all','V1v','V1d','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','TO1','TO2','hV4','VO1','VO2','PHC1','PHC2','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF','IFS'};

save('AreaDistances','allrois','area_labels','alluniqueroinodepairs','subj_folders','allareadistancessmoothwm_rh','allareadistancespial_rh','allareadistancesavglayer_rh','allsubjs_alluniquerois_rh','allareadistancessmoothwm_lh','allareadistancespial_lh','allareadistancesavglayer_lh','allsubjs_alluniquerois_lh')
