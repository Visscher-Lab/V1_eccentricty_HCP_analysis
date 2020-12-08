clear all;
addpath('/data/project/vislab/a/HCP_diff_Sara/pinar')
outputdir = '/data/project/vislab/a/HCP_diff_Sara/UP_funcanalysis/data';
fsdir='/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheckcopy';
setenv('SUBJECTS_DIR',fsdir);
parentdir= '/data/project/vislab/raw/HCP_900sub/HCP_func';
cd(outputdir);
subs = dir('*');
subs = subs(3:length(subs)-3);

for r=1:length(subs)
    if exist ([parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL']) dir
        system(['bbregister --s fsaverage --mov ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/' subs(r).name '.nii --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --init-fsl --bold']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/central/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_central.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/central/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_central.nii.gz --hemi rh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/mid/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_mid.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/mid/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_mid.nii.gz --hemi rh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/far/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_far.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/far/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_far.nii.gz --hemi rh']);
    end
    if exist ([parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL']) dir
        system(['bbregister --s fsaverage --mov ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/' subs(r).name '.nii --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --init-fsl --bold']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/central/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_central.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/central/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_central.nii.gz --hemi rh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/mid/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_mid.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/mid/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_mid.nii.gz --hemi rh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/far/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/lh_far.nii.gz --hemi lh']);
        system(['mri_vol2surf --mov ' outputdir '/' subs(r).name '/far/SCA_result.nii.gz --trgsubject fsaverage --reg ' outputdir '/' subs(r).name '/' subs(r).name '.dat --o ' outputdir '/' subs(r).name '/rh_far.nii.gz --hemi rh']);
    end
end
