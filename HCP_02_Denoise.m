% This script runs denoising for the spatially normalized (but unsmoothed)
% resting state VINES data by calling the fc_denoising function. For more
% information on denoising, see fc_denoising. Note that fc_denoising calls
% another function called vl_filter to perform filtering. See that function
% for more information about filtering procedures

% This script was written for submission to Cheaha so that all processing
% could be done in parallel for time efficiency. Thus, there is a flag to
% get a subject number that is a Cheaha reference.

% March, 2019
% Pinar Demirayak


% Initialization

addpath('/share/apps/vbc/matlab/toolbox/spm12/'); % path where SPM lives
addpath('/share/apps/vbc/matlab/toolbox/matlab_nifti'); % add path to necessary functions (load_nii)
addpath('/data/project/vislab/a/HCP_diff_Sara/pinar'); % add path to necessary functions (fc_denoising/vl_filter)
parentdir= '/data/project/vislab/raw/HCP_900sub/HCP_func'; % directory where subject data lives
cd(parentdir);
subs = dir('*'); % gets a list of all the subjects
subs = subs(3:length(subs));
cfg = struct; % preallocate cfg structure
%cfg.mask_bold = '/share/apps/vbc/matlab/toolbox/spm8/tpm/grey.nii';
%cfg.bold_cutoff = 0.4; % cutoff threshold for grey matter mask
%cfg.mask_wm = '/share/apps/vbc/matlab/toolbox/spm8/tpm/white.nii';
%cfg.mask_csf = '/share/apps/vbc/matlab/toolbox/spm8/tpm/csf.nii';
cfg.TR = 0.7; % the repetition time
cfg.fwhm = [8 8 8]; % smoothing kernal in mm

% Loop for the resting state scans
for r = 1:length(subs)
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST2_LR']) dir
      
        cd ([pwd '/MNINonLinear/Results/rfMRI_REST2_LR']);
            importtxt = load('Movement_Regressors.txt');
            exporttxt(:,4:6) = deg2rad(importtxt(:,4:6));
            exporttxt(:,1:3)=importtxt(:,1:3);
            fid=fopen('Movement_Regressors2.txt','wt');
            fprintf(fid, '%.6f %.6f %.6f %e %e %e \n',exporttxt.');
            fclose(fid);
            cd([parentdir '/' subs(r).name]);
            cfg.prefix = subs(r).name;
            %cd(subs(r).name); % navigates into the folder containing the rth resting state scan
            disp(['Beginning scan ', num2str(r), ' of ', num2str(length(subs))]);
            cfg.path_bold = [pwd '/MNINonLinear/Results/rfMRI_REST2_LR/'];
            %mkdir(['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_LR/']);
            %cfg.path_correl = ['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_LR/'];
            cfg.path_move = [pwd '/MNINonLinear/Results/rfMRI_REST2_LR/']; % specify path to movement file as current path
            fc_denoising_new2(cfg); % perform denoising
            cd(parentdir);
   %     else
    %        cd(parentdir);
  %      end
         else
            cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST1_LR']) dir
      %  if ~isempty('rfMRI_REST1_LR')
      cd ([pwd '/MNINonLinear/Results/rfMRI_REST1_LR']);
            importtxt = load('Movement_Regressors.txt');
            exporttxt(:,4:6) = deg2rad(importtxt(:,4:6));
            exporttxt(:,1:3)=importtxt(:,1:3);
            fid=fopen('Movement_Regressors2.txt','wt');
            fprintf(fid, '%.6f %.6f %.6f %e %e %e \n',exporttxt.');
            fclose(fid);
            cd([parentdir '/' subs(r).name]);
            cfg.prefix = subs(r).name;
            %cd(subs(r).name); % navigates into the folder containing the rth resting state scan
            disp(['Beginning scan ', num2str(r), ' of ', num2str(length(subs))]);
            cfg.path_bold = [pwd '/MNINonLinear/Results/rfMRI_REST1_LR/'];
            %mkdir(['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_LR/']);
            %cfg.path_correl = ['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_LR/'];
            cfg.path_move = [pwd '/MNINonLinear/Results/rfMRI_REST1_LR/']; % specify path to movement file as current path
            fc_denoising_new2(cfg); % perform denoising
            cd(parentdir);
%         else
%             cd(parentdir);
%         end
         else
            cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST2_RL']) dir
       % if ~isempty('rfMRI_REST2_RL')
       cd ([pwd '/MNINonLinear/Results/rfMRI_REST2_RL']);
            importtxt = load('Movement_Regressors.txt');
            exporttxt(:,4:6) = deg2rad(importtxt(:,4:6));
            exporttxt(:,1:3)=importtxt(:,1:3);
            fid=fopen('Movement_Regressors2.txt','wt');
            fprintf(fid, '%.6f %.6f %.6f %e %e %e \n',exporttxt.');
            fclose(fid);
            cd([parentdir '/' subs(r).name]);
            cfg.prefix = subs(r).name;
            %cd(subs(r).name); % navigates into the folder containing the rth resting state scan
            disp(['Beginning scan ', num2str(r), ' of ', num2str(length(subs))]);
            cfg.path_bold = [pwd '/MNINonLinear/Results/rfMRI_REST2_RL/'];
            %mkdir(['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/']);
            %cfg.path_correl = ['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/'];
            cfg.path_move = [pwd '/MNINonLinear/Results/rfMRI_REST2_RL/']; % specify path to movement file as current path
            fc_denoising_new2(cfg); % perform denoising
            cd(parentdir);
%         else
%             cd(parentdir);
%         end
         else
            cd(parentdir);
    end
    cd(subs(r).name);
    if exist ([pwd '/MNINonLinear/Results/rfMRI_REST1_RL']) dir
      %  if ~isempty('rfMRI_REST1_RL')
      cd ([pwd '/MNINonLinear/Results/rfMRI_REST1_RL']);
            importtxt = load('Movement_Regressors.txt');
            exporttxt(:,4:6) = deg2rad(importtxt(:,4:6));
            exporttxt(:,1:3)=importtxt(:,1:3);
            fid=fopen('Movement_Regressors2.txt','wt');
            fprintf(fid, '%.6f %.6f %.6f %e %e %e \n',exporttxt.');
            fclose(fid);
            cd([parentdir '/' subs(r).name]);
            cfg.prefix = subs(r).name;
            %cd(subs(r).name); % navigates into the folder containing the rth resting state scan
            disp(['Beginning scan ', num2str(r), ' of ', num2str(length(subs))]);
            cfg.path_bold = [pwd '/MNINonLinear/Results/rfMRI_REST1_RL/'];
            %mkdir(['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/']);
            %cfg.path_correl = ['/data/user/pinarde/HCP_900/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/'];
            cfg.path_move = [pwd '/MNINonLinear/Results/rfMRI_REST1_RL/']; % specify path to movement file as current path
            fc_denoising_new2(cfg); % perform denoising
            cd(parentdir);
%         else
%             cd(parentdir);
%         end
         else
            cd(parentdir);
    end
end
