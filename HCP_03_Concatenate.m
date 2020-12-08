% This script checks the denoised scans to see how many volumes remain
% after removing volumes in which there was too much motion. If there are
% enough volumes remaining (>=150 volumes or >=5 min) then all remaining
% volumes are concatenated and the new "concatenated" image is saved in a
% new folder in that subject's directory.

% This script was written for submission to Cheaha so that all processing
% could be done in parallel for time efficiency. Thus, there is a flag to
% get a subject number that is a Cheaha reference.

% March, 2019
% Pinar Demirayak

% Initialization
clear all;
addpath('/share/apps/vbc/matlab/toolbox/matlab_nifti') % add path to necessary functions (load_nii)
addpath('/data/project/vislab/a/HCP_diff_Sara/pinar'); % add path to necessary functions (fc_denoising/vl_filter)
parentdir = '/scratch/pinarde'%'/data/project/vislab/raw/HCP_900sub/HCP_func'; % directory where subject data lives
outputdir = '/data/project/vislab/a/HCP_diff_Sara/UP_funcanalysis/data';
cd(parentdir);
subs = dir('*'); % gets a list of all the subjects
subs = subs(3:length(subs));

% loop for the resting state scans
data = []; % preallocates an array to store data

% Loop to concatenate
for r = 1:1%length(subs)
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST2_RL/']) dir % navigates into the folder containing the rth resting state scan
        cd([pwd '/MNINonLinear/Results/rfMRI_REST2_RL/']);
        load fc_denoising_results.mat % load denoising results
        badvolumes = fc_denoising_results.discarded_timepoints; % get the "bad" volumes that were removed
        goodvolumes = (1:1200)'; % initialize a list of good volumes (1200 because that's how many original time points there were)
        goodvolumes(badvolumes) = []; % remove the bad volumes from the good list
        % loop to see if 1st 4 volumes have been removed (magnetic field may
        % not have necessarily been stabilized before then)
        for j = 1:4
            ind(j) = isempty(find(goodvolumes == j)); % if the jth volume was removed, then ind is 1 and if it wasnt then ind is 0
        end
        ind = 1-ind; % turns 1s into 0s and vice versa; now 1 means volume wasn't remved
        indsum = sum(ind); % takes sum of ind to see how many of the first 4 volumes need to be removed
        scan = dir('s*denoised.nii'); % gets name of scan
        image = load_nii(scan.name); % loads scan/image
        image.img(:,:,:,1:indsum) =[]; % removes correct first few volumes
        data = cat(4,data,image.img); % temporally concatenates data (basically 1st and 2nd scans but handles if there is only 1)
        % Save image
        %mkdir Concatenated_RestingState % creates new folder
        %cd Concatenated_RestingState % navigates into that folder
        image.hdr.dime.dim(5) = size(data,4); % change 4th dimension to have correct length
        image.img = data; % saves data into .img portion of the image
        save_nii(image, 'concatenated_rs.nii'); % saves the nifti
        cd(parentdir);
    else
        cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST2_LR/']) dir % navigates into the folder containing the rth resting state scan
        cd([pwd '/MNINonLinear/Results/rfMRI_REST2_LR/']);
        load fc_denoising_results.mat % load denoising results
        badvolumes = fc_denoising_results.discarded_timepoints; % get the "bad" volumes that were removed
        goodvolumes = (1:1200)'; % initialize a list of good volumes (180 because that's how many original time points there were)
        goodvolumes(badvolumes) = []; % remove the bad volumes from the good list
        % loop to see if 1st 4 volumes have been removed (magnetic field may
        % not have necessarily been stabilized before then)
        for j = 1:4
            ind(j) = isempty(find(goodvolumes == j)); % if the jth volume was removed, then ind is 1 and if it wasnt then ind is 0
        end
        ind = 1-ind; % turns 1s into 0s and vice versa; now 1 means volume wasn't removed
        indsum = sum(ind); % takes sum of ind to see how many of the first 4 volumes need to be removed
        scan = dir('s*denoised.nii'); % gets name of scan
        image = load_nii(scan.name); % loads scan/image
        image.img(:,:,:,1:indsum) =[]; % removes correct first few volumes
        data = cat(4,data,image.img); % temporally concatenates data (basically 1st and 2nd scans but handles if there is only 1)
        % Save image
        %mkdir Concatenated_RestingState % creates new folder
        %cd Concatenated_RestingState % navigates into that folder
        image.hdr.dime.dim(5) = size(data,4); % change 4th dimension to have correct length
        image.img = data; % saves data into .img portion of the image
        save_nii(image, 'concatenated_rs.nii'); % saves the nifti
        mkdir([outputdir '/' subs(r).name]);
        cd([outputdir '/' subs(r).name]);
        system(['fslmerge -t ' outputdir '/' subs(r).name '/Concatenated_RestingState ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/concatenated_rs ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_LR/concatenated_rs']);
        cd(parentdir);
    else
        cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST1_RL/']) dir % navigates into the folder containing the rth resting state scan
        cd([pwd '/MNINonLinear/Results/rfMRI_REST1_RL/']);
        load fc_denoising_results.mat % load denoising results
        badvolumes = fc_denoising_results.discarded_timepoints; % get the "bad" volumes that were removed
        goodvolumes = (1:1200)'; % initialize a list of good volumes (180 because that's how many original time points there were)
        goodvolumes(badvolumes) = []; % remove the bad volumes from the good list
        % loop to see if 1st 4 volumes have been removed (magnetic field may
        % not have necessarily been stabilized before then)
        for j = 1:4
            ind(j) = isempty(find(goodvolumes == j)); % if the jth volume was removed, then ind is 1 and if it wasnt then ind is 0
        end
        ind = 1-ind; % turns 1s into 0s and vice versa; now 1 means volume wasn't remved
        indsum = sum(ind); % takes sum of ind to see how many of the first 4 volumes need to be removed
        scan = dir('s*denoised.nii'); % gets name of scan
        image = load_nii(scan.name); % loads scan/image
        image.img(:,:,:,1:indsum) =[]; % removes correct first few volumes
        data = cat(4,data,image.img); % temporally concatenates data (basically 1st and 2nd scans but handles if there is only 1)
        % Save image
        %mkdir Concatenated_RestingState % creates new folder
        %cd Concatenated_RestingState % navigates into that folder
        image.hdr.dime.dim(5) = size(data,4); % change 4th dimension to have correct length
        image.img = data; % saves data into .img portion of the image
        save_nii(image, 'concatenated_rs.nii'); % saves the nifti
        cd(parentdir);
    else
        cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST1_LR/']) dir % navigates into the folder containing the rth resting state scan
        cd([pwd '/MNINonLinear/Results/rfMRI_REST1_LR/']);
        load fc_denoising_results.mat % load denoising results
        badvolumes = fc_denoising_results.discarded_timepoints; % get the "bad" volumes that were removed
        goodvolumes = (1:1200)'; % initialize a list of good volumes (180 because that's how many original time points there were)
        goodvolumes(badvolumes) = []; % remove the bad volumes from the good list
        % loop to see if 1st 4 volumes have been removed (magnetic field may
        % not have necessarily been stabilized before then)
        for j = 1:4
            ind(j) = isempty(find(goodvolumes == j)); % if the jth volume was removed, then ind is 1 and if it wasnt then ind is 0
        end
        ind = 1-ind; % turns 1s into 0s and vice versa; now 1 means volume wasn't remved
        indsum = sum(ind); % takes sum of ind to see how many of the first 4 volumes need to be removed
        scan = dir('s*denoised.nii'); % gets name of scan
        image = load_nii(scan.name); % loads scan/image
        image.img(:,:,:,1:indsum) =[]; % removes correct first few volumes
        data = cat(4,data,image.img); % temporally concatenates data (basically 1st and 2nd scans but handles if there is only 1)
        % Save image
        %mkdir Concatenated_RestingState % creates new folder
        %cd Concatenated_RestingState % navigates into that folder
        image.hdr.dime.dim(5) = size(data,4); % change 4th dimension to have correct length
        image.img = data; % saves data into .img portion of the image
        save_nii(image, 'concatenated_rs.nii'); % saves the nifti
        mkdir([outputdir '/' subs(r).name]);
        cd([outputdir '/' subs(r).name]);
        system(['fslmerge -t ' outputdir '/' subs(r).name '/Concatenated_RestingState ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/concatenated_rs ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_LR/concatenated_rs']);
        cd(parentdir);
    else
        cd(parentdir);
    end
    clearvars -except parentdir outputdir subs data;
    data=[];
end
