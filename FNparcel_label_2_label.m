%% Set up steps
%  Be sure to copy and paste the actual freesurfer fsaverage folder into
%  your subject directory (the one it creates in the subject directoy is actually a link to the fs
%  install location). Take the label files and move them into the
%  fsaverage/labels folder. There are a few locations to make changes in the
%  code for this to run. I have indicated where to make these changes.


%% set the enviroment variable to your subjects directory

setenv('SUBJECTS_DIR','/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheckcopy'); % be sure to edit this to your path

%% cd to subject directory and generate a list of subject names
cd /data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheckcopy;
Subs = dir('**'); % change this to a base name that will pull all participants
%
%% Label names
%  These should be located in your subject directory/fsaverage/label folder.
%  Be sure you copy and pasted the fsaverage folder to your subject
%  directory and it is not just the symbolic link

labelslh = {'3_LH_V1.label' '5_LH_V1.label' '7_LH_V1.label'...
'2_LH_V1.label' '4_LH_V1.label' '6_LH_V1.label' '8_LH_V1.label' }; % all left hemi labels

labelsrh = { '3_RH_V1.label' '5_RH_V1.label' '7_RH_V1.label'...
'2_RH_V1.label' '4_RH_V1.label' '6_RH_V1.label' '8_RH_V1.label' }; % all right hemi labels