clear all
addpath(genpath('afni_matlab'))
species='Humans';
computer_name='marcaro';

data_dir=['/Users/' num2str(computer_name) '/Box/RetinotopicTopologyProject/Humans/'];
save_folder = ['/Users/' num2str(computer_name) '/Box/RetinotopicTopologyProject/Graphs/']; %folder where graphs are saved

files = dir(data_dir);
names = {files.name}; 
dirFlags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..') & ~strcmp(names, 'phasemaps') & ~strcmp(names, 'VTPM') & ~strcmp(names, 'VTPM_delayedsaccade') & ~strcmp(names, 'ignore');
subjs=names(dirFlags);


% Labels for allvisual-
% IMPORTANT: some subjects might not have all ROIs
allvisual_numbers=[1,2,3,4,5,6,50,51,52,100,101,102,103,150,151,152,153,154];
allvisual_labels={'V1v','V1d','V2v','V2d','V3v','V3d','V3A','V3B','V7','LO1','LO2','TO1','TO2','hV4','VO1','VO2','PHC1','PHC2'};

nodearea_smoothwm_col=2;
nodearea_pial_col=3;
triarea_smoothwm_col=4;
triarea_pial_col=5;
nodevol_col=6;
nodethickness_col=7;

cd(num2str(species))

norois=[];
wrongrois=[];
for curr_subj = 1:size(subjs,2)
    cd(data_dir)
    cd(num2str(subjs{curr_subj}))
    fprintf(['Checking ' num2str(subjs{curr_subj}) '\n'])
    
    allvisual_rh=[];
    allvisual_lh=[];
    if exist(['rois/allvisual-rh.1D.dset'])
        allvisual_rh=Read_1D(['rois/allvisual-rh.1D.dset']);
        allvisual_lh=Read_1D(['rois/allvisual-lh.1D.dset']);
        
        redflag=[];
        redflag=[sum(allvisual_rh(:,2)==9)+sum(allvisual_rh(:,2)==8)];
        if redflag>0
           fprintf(['Errors with ' num2str(subjs{curr_subj}) '\n'])
           wrongrois=[wrongrois,curr_subj]; 
        end
        
        
        
    else
        fprintf(['No allvisual in ' num2str(subjs{curr_subj}) '\n']);
        norois=[norois,curr_subj];
    end
    
    
end