
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc_denoising_results = fc_denoising_new(cfg)

% FC_DENOISING Preprocesses resting state fMRI data for connectivity analysis,
% using the following steps: 
%     (1) Nuisance Regression using motion
%         parameters from rp_*.txt files as regressor (from realignment 
%         in SPM12) 
%     (2) Identifies "bad" scan time points using Framewise Displacement 
%         and DVARS (see Powers, et al. "Spurious but systemic correlations 
%         in functional connectivity MRI networks arise from subject 
%         motion". Neuroimage. 2012 Feb 1;59(3):2142-54.
%     (3) Interpolates "bad" scan time points using time points before and
%         after (see Carp, J. "Optimizing the order of operations for 
%         movement scrubbing: Comment on Power et al." Neuroimage. 2012)
%     (4) Bandpass filter data (Temporal smoothing)
%     (5) Temporally mask "bad" scans (Volumes that contains movement).
%     (6) Nuisance Regression using WM and CSF principle component scores
%         as regressors.
%     (7) Spatially smooth data using SPM toolbox. 
% --------------------------------------------------------------------
%   The image used as input should have been preprocessed with common steps
%   such as : realignement, slice-timing, coregistration,
%   Normalization. However, no spatial smoothing sould have been applied,
%   this function will take care of this step in order to optimize the
%   denoising process.
% --------------------------------------------------------------------
% The configuration should contain:
%   cfg.path_bold : folder where the bold image to denoise is. 
%   cfg.path_move : folder where the SPM mouvement file is.
%   cfg.path_correl : folder where the denoised image will be saved. If not
%                     specified, cfg.path_bold will be used.
%   cfg.prefix : Bold image prefix. This prefix will be used to identify
%                the file
%   cfg.filter_type : Type of filter to apply. 'lpfilter' for lowpass,
%                     'hpfilter' for highpass and 'bpfilter' for bandpass.
%   cfg.filter_algo : Type of filter algorithm. 'but' for butterworth,
%                     'fir' for FIR.
%   cfg.filter_pass : 'onepass' for one-pass, 'onepass-reverse' for a
%                     reverse one-pass and 'twopass' for a two-way filter.
%   cfg.filter_order : Filter order. 2 is default value.
%   cfg.TR : TR of BOLD images.
%   cfg.freq : Filter frequency cut-off. [0.01 0.08] will be used as
%              default.
%   cfg.csf_cutoff : Minimum voxel value to take into account for CSF mask.
%   cfg.wm_cutoff : Minimum voxel value to take into account for White-matter mask
%   cfg.mvt_thresh : Movement threshold in mm (A volume showing a movement
%                    component above this threshold will be corrected and deleted)
%   cfg.dvars_thresh : DVARS threshold. 
%   cfg.pca_thresh : Minimum percentage of variance explained by
%                    white-matter signal and CSF signal principal components
%   cfg.fhwm : Spatial smoothing gaussian filter values. [6 6 6] is default
%               value. Set a 0 for no spatial smoothing
% 
% OUTPUT : 
%   Creates 2 images with suffix "denoised", one denoised but non spatially smoothed, and another
%   one denoised and spatially smoothed (with prefix "s" added by SPM)
%
% Copyright (C) Rodolphe Nenert (Visscher Lab) & Jennifer Hadley (Lahti
% Lab)
% February 2012. Univeristy of Alabama at Birmingham.
%
% Version Notes :
% 28 February 2012 : First version of the function, with header (RN & JH)
% 21 March 2014 : Added some clearvars to optimize memory usage (RN)

%% Getting useful variables
if ~isfield(cfg,'path_bold'), error('Path to bold image not specified '); else path_bold = cfg.path_bold; end
if ~isfield(cfg,'path_move'), error('Path to movement file not specified '); else path_move = cfg.path_move; end
if ~isfield(cfg,'path_correl'), fprintf('Path to save denoised image not specified, same path as bold image will be used \n'); path_correl = cfg.path_bold; else path_correl = cfg.path_correl; end
if ~isfield(cfg,'prefix'), error('Bold image prefix not specified '); else prefix = cfg.prefix; end
if ~isfield(cfg,'filter_type'), fprintf('Filter type not specified, bandpass will be used as default \n'); filter_type = 'bpfilter'; else filter_type = cfg.filter_type; end
if ~isfield(cfg,'filter_algo'), fprintf('Filter algorithm not specified, Butterworth will be used as default \n'); filter_algo = 'but'; else filter_algo = cfg.filter_algo; end
if ~isfield(cfg,'filter_pass'), fprintf('Filter pass not specified, onepass will be used as default \n'); filter_pass = 'onepass'; else filter_pass = cfg.filter_pass; end
if ~isfield(cfg,'filter_order'), fprintf('Filter order not specified, 2 will be used as default \n'); N = 2; else N = cfg.filter_order; end
if ~isfield(cfg,'TR'), error('TR not specified'); else fsample = 1/(cfg.TR); end
if ~isfield(cfg,'freq'), fprintf('Filter frequency cut-off not specified, [0.01 0.08] will be used as default \n'); freq = [0.01 0.08]; else freq = cfg.freq; end
if ~isfield(cfg,'csf_cutoff'), fprintf('CSF mask cut-off not specified, 0.9 will be used as default \n'); csf_cutoff = 0.95; else csf_cutoff = cfg.csf_cutoff; end
if ~isfield(cfg,'wm_cutoff'), fprintf('White-matter mask cut-off not specified, 0.6 will be used as default \n'); wm_cutoff = 0.95; else wm_cutoff = cfg.wm_cutoff; end
if ~isfield(cfg,'grey_cutoff'), fprintf('Grey matter mask cut-off not specified, 0 will be used as default \n'); bold_cutoff = 0; else bold_cutoff = cfg.grey_cutoff; end
if ~isfield(cfg,'mvt_thresh'), fprintf('movement threshold not specified, 0.5 will be used as default \n'); mvt_thresh = 0.5; else mvt_thresh = cfg.mvt_thresh; end
if ~isfield(cfg,'pca_thresh'), fprintf('PCA threshold not specified, 0.95 will be used as default \n'); pca_thresh = 0.95; else pca_thresh = cfg.pca_thresh; end
if ~isfield(cfg,'fwhm'), fprintf('Spatial smoothing not specified, [6 6 6] will be used as default \n'); fwhm = [6 6 6]; else fwhm = cfg.fwhm; end

%- Start analysis
%--------------------------------------------------------------------------
d = spm('dir');
rTPMpath='/data/project/vislab/a/HCP_diff_Sara/pinar';
myprob = strcat(rTPMpath,filesep,'tpm',filesep,'rTPM.nii');
V = spm_vol(myprob);
wmfile = spm_read_vols(V(2));% Load White_matter mask
csffile = spm_read_vols(V(3)); % Load CSF mask
greyfile = spm_read_vols(V(1)); % Load global mask

boldfile=spm_select('ExtFPList',fullfile(path_bold), prefix); % Bold image
V = spm_vol(boldfile);
boldimg = spm_read_vols(V);
boldimg(isnan(boldimg)) = 0 ; % Replace all NaN is bold image by zeros.

fprintf(strcat('Start denoising ! \n'));

[x,y,z] = ind2sub(size(greyfile),find(greyfile > bold_cutoff)); % get coordinates of voxels with non-null value in the mask.
if isempty(x), error('Bold cutoff value is too high, list of voxels is null'); end 
[x2,y2,z2] = ind2sub(size(csffile),find(csffile > csf_cutoff)); % get coordinates of voxels with values above threshold in CSF image
if isempty(x2), error('CSF cutoff value is too high, list of voxels is null'); end 
[x3,y3,z3] = ind2sub(size(wmfile),find(wmfile > wm_cutoff)); % get coordinates of voxels with values above threshold  in WM image
if isempty(x3), error('White matter cutoff value is too high, list of voxels is null'); end 

% Preallocate arrays
imgnum = size(boldimg,4); % Get number of timepoints (fourth dimension of bold image)
voxel_timecourse = zeros(length(x),imgnum); % Will contain 2D bold timecourse (masked)
voxel_timecourse2 = zeros(length(x2),imgnum); % Will contain 2D CSF timecourse (masked)
voxel_timecourse3 = zeros(length(x3),imgnum); % Will contain 2D WM timecourse (masked)
voxel_tc_g = zeros(length(x),imgnum); % Will contain filtered voxel timecourse of the bold
voxel_tc_c = zeros(length(x2),imgnum); % Will contain filtered voxel timecourse of CSF
voxel_tc_w = zeros(length(x3),imgnum); % Will contain filtered voxel timecourse of WM

% Transformation of 4D data into a 2D matrix, easier to manipulate.
fprintf('\t Preamble : Timecourses extractions \n');

fprintf('\t\t Extracting BOLD voxel timecourses...');
for i=1:length(x)
   voxel_timecourse(i,:) = boldimg(x(i),y(i),z(i),:); 
end
fprintf('Done \n');
fprintf('\t\t Extracting CSF voxel timecourses...');
for i=1:length(x2)
   voxel_timecourse2(i,:) = boldimg(x2(i),y2(i),z2(i),:);
end
fprintf('Done \n');
fprintf('\t\t Extracting white matter voxel timecourses...');
for i=1:length(x3)
   voxel_timecourse3(i,:) = boldimg(x3(i),y3(i),z3(i),:);
end
fprintf('Done \n');

%- Regressing out movement
%--------------------------------------------------------------------------

fprintf('\t STEP 1 : Regressing out movement \n');
cd(path_move)
temp = dir('Movement_Regressors2.txt');
filemvt = load(temp(1).name); % Load files that contains movement parameters
filemvt_deriv= vertcat(zeros(1,size(filemvt,2)),diff(filemvt)); % Pull out first derivative from movement
pred = [ones(length(filemvt_deriv),1) filemvt filemvt_deriv]; % Constructing predictors matrix

fprintf('\t\t On BOLD voxels:  ');
temp='';
parfor i=1:length(x)
    resp = voxel_timecourse(i,:)';
    [~,~,r] = regress(resp,pred); % Calculate multiple regression, residuals = new voxel timecourse;
    voxel_tc_g(i,:) = r;
%     fprintf(repmat('\b',1,length(temp)));
%     temp = sprintf('%3.4f %% ',i/length(x)*100);
%     fprintf('%s',temp);
end
fprintf('...Done \n');

fprintf('\t\t On CSF voxels:  ');
temp='';
parfor i=1:length(x2)
    resp = voxel_timecourse2(i,:)';
    [~,~,r] = regress(resp,pred); % Calculate multiple regression, residuals = new voxel timecourse;
    voxel_tc_c(i,:) = r;
%     fprintf(repmat('\b',1,length(temp)));
%     temp = sprintf('%3.4f %% ',i/length(x2)*100);
%     fprintf('%s',temp);
end
fprintf('...Done \n');

fprintf('\t\t On White-matter voxels:  ');
temp='';
parfor i=1:length(x3)
    resp = voxel_timecourse3(i,:)';
    [~,~,r] = regress(resp,pred); % Calculate multiple regression, residuals = new voxel timecourse;
    voxel_tc_w(i,:) = r;
%     fprintf(repmat('\b',1,length(temp)));
%     temp = sprintf('%3.4f %% ',i/length(x3)*100);
%     fprintf('%s',temp);
end
fprintf('...Done \n');
clearvars voxel_timecourse voxel_timecourse2 voxel_timecourse3; % Let's free some RAM


%- STEP 2 Find bad volumes with high movement effect on BOLD
%--------------------------------------------------------------------------

fprintf('\t STEP 2 : Search for bad volumes \n');
% First, find find bad volumes based on spm movement values.
fprintf('\t\t On movement values...');
bad_trans = abs(filemvt_deriv(:,1)) + abs(filemvt_deriv(:,2)) + abs(filemvt_deriv(:,3));
bad_rotat = (abs(filemvt_deriv(:,4)) + abs(filemvt_deriv(:,5)) + abs(filemvt_deriv(:,6))) .*50;
bad_all = bad_trans + bad_rotat;
bad_id = find(bad_all > mvt_thresh);
fprintf('...Done \n');
num_bad = round(((length(bad_id)/imgnum)*2)*imgnum); % Number of bad volumes to take into account for DVARS
if (num_bad >= imgnum)
    num_bad = imgnum - 1;
end

% Second, DVARS-like analysis
fprintf('\t\t DVARS analysis:  ');
temp_tc = vertcat(zeros(1,length(x)),diff(voxel_tc_g')); % DVARS calculation (Powers et al.)
temp_tc2 = temp_tc.*temp_tc;
temp_tc_avg = sqrt(mean(temp_tc2'))*10;
sort_avg = sort(temp_tc_avg,'descend');
dvars_thresh = sort_avg(1,num_bad + 1);
bad_id2 = find(temp_tc_avg > dvars_thresh);   % Find volumes above threshold. 
bad_volumes = intersect(bad_id,bad_id2);
fprintf('...Done \n');
clearvars temp_tc temp_tc2; % Let's free some RAM again


%- STEP 3 Interpolation of bad volumes.
%--------------------------------------------------------------------------

fprintf('\t STEP 3 : Bad volumes interpolation ');
for j=1:length(bad_volumes)
    if (bad_volumes(j) + 1 < imgnum)
        voxel_tc_g(:,bad_volumes(j)) = (voxel_tc_g(:,bad_volumes(j)-1) + voxel_tc_g(:,bad_volumes(j) + 1))./2;
    else
        voxel_tc_g(:,bad_volumes(j)) = voxel_tc_g(:,bad_volumes(j)-1);
    end
end
fprintf('...Done \n');


%- STEP 4 Temporal smoothing
%--------------------------------------------------------------------------

fprintf('\t STEP 4 : Temporal smoothing \n');
fprintf('\t\t Filtering timecourses...');
voxel_timecoursef = vl_filter(voxel_tc_g,filter_type,filter_algo,filter_pass,freq,fsample,N); % filter all voxels timecourse
voxel_timecourse2f = vl_filter(voxel_tc_c,filter_type,filter_algo,filter_pass,freq,fsample,N); % filter voxels from CSF timecourse
voxel_timecourse3f = vl_filter(voxel_tc_w,filter_type,filter_algo,filter_pass,freq,fsample,N); % filter voxels from WM timecourse
fprintf('Done \n');
clearvars voxel_tc_g voxel_tc_c voxel_tc_w; % Free some memory


%- STEP 5 Removing bad volumes
%--------------------------------------------------------------------------

fprintf('\t STEP 5 : Bad volumes removal \n');
voxel_timecoursef(:,bad_volumes) = [];
voxel_timecourse2f(:,bad_volumes) = [];
voxel_timecourse3f(:,bad_volumes) = [];


%- STEP 6 Regressing out WM and CSF
%--------------------------------------------------------------------------

[~,score1,latent1] = princomp(voxel_timecourse2f', 'econ');  % PCA over CSF timecourses
var1 = cumsum(latent1)./sum(latent1);
var_id = find(var1 > pca_thresh); % Find components that explains the minimum of variance specified in cfg
csf_comp = score1(:, 1:var_id(1,1));
clear voxel_timecourse2f

[~,score2,latent2] = princomp(voxel_timecourse3f', 'econ'); % PCA over WM timecourses
var2 = cumsum(latent2)./sum(latent2);
var_id = find(var2 > pca_thresh); % Find components that explains the minimum of variance specified in cfg
wm_comp = score2(:, 1:var_id(1,1));
clear voxel_timecourse3f

pred = [ones(length(csf_comp),1) csf_comp wm_comp]; % Build predictors matrix with components

fprintf('\t\t Getting residuals:  ');
voxel_timecourseff = zeros(size(voxel_timecoursef));
temp='';
parfor i=1:length(voxel_timecoursef) % iterate over columns of matrix
    resp = voxel_timecoursef(i,:)';
    [~,~,r] = regress(resp,pred); % Calculate multiple regression, residuals = new voxels timecourses;
    voxel_timecourseff(i,:) = r;
%     fprintf(repmat('\b',1,length(temp)));
%     temp = sprintf('%3.4f %% ',i/length(x)*100);
%     fprintf('%s',temp);
end
fprintf('Done \n');
clearvars score1 latent1 score2 latent2 csf_comp wm_comp voxel_timecoursef voxel_timecourse2f voxel_timecourse3f;


%- STEP 7 SAVE DENOISED IMAGE + SPATIAL SMOOTHING
%--------------------------------------------------------------------------

denoised_blank = zeros(size(boldimg,1),size(boldimg,2),size(boldimg,3),size(voxel_timecourseff,2)); % Filled with zeros, to be refilled with correlation coeff later
for i=1:length(x)
    denoised_blank(x(i),y(i),z(i),:) = voxel_timecourseff(i,:); % Filling blank image with timecourses values
end

newV = spm_create_vol(V);
newV = newV(1:size(voxel_timecourseff,2));
newV = rmfield(newV,'pinfo');
[path,file,ext] = fileparts(V(1).fname);
for i=1:size(voxel_timecourseff,2)
    newV(i).fname = strcat(path,filesep,file,'_denoised',ext);
    newV(i).dt = [64,0];
    spm_write_vol(newV(i),denoised_blank(:,:,:,i));
end

% Last step, spatial smoothing 
fprintf('\t STEP 7 : Spatial smoothing using SPM function \n');
mynorm=cellstr(spm_select('ExtFPList',fullfile(path), '_denoised')); % ref = functional image
matlabbatch=[];
matlabbatch{1}.spm.spatial.smooth.data = mynorm;
matlabbatch{1}.spm.spatial.smooth.fwhm = fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);

% Save variables of interest
fc_denoising_results = [];
fc_denoising_results.fd = bad_all;
fc_denoising_results.fd_threshold = mvt_thresh;
fc_denoising_results.dvars = temp_tc_avg;
fc_denoising_results.discarded_timepoints = bad_volumes;
save fc_denoising_results.mat fc_denoising_results;
% Done!
fprintf('\t Hurray! Image successfully denoised \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtered = vl_filter(timecourse,filter_type, filter_algo,filter_pass,freq,freq_sample,N)

% vl_filter applies a filter on one or several timecourses in a matrix
% (Nsample x Ndatapoints)
%
% INPUTS:
%   - timecourse: matrix containing each timecourse on a line
%   - filter_type: 'lpfilter' for lowpass, 'hpfilter' for highpass and 'bpfilter' for bandpass
%   - filter_algo: 'but' for butterworth, 'fir' for FIR
%   - filter_pass: 'onepass' for one-pass, 'onepass-reverse' for a reverse one-pass and 'twopass' for a two-way filter
%   - freq: cut-off frequency, put a table for bandpass filter
%   - freq_sample: sample frequency
%   - N: filter order (put default values if non-specified), specify the
%   amount of signal reduction as following : Nx6 db per octave
%
% Copyright (c) 2011
% AUTHOR: Rodolphe Nenert

% If order non-specified, put an empty value
if nargin < 7
    N = [];
end

Fn = freq_sample/2; % Nyquist frequency

%% Calculate coefficients according to the filter type
switch filter_type
    
    case 'lpfilter' % low-pass filter    
        switch filter_algo
            case 'but'
                 if isempty(N)
                    N = 6;
                 end
                [B, A] = butter(N, max(freq)/Fn); % Butterworth coefficients calculation;
            case 'fir'
                if isempty(N)
                    N = 25; 
                end
                [B, A] = fir1(N, max(freq)/Fn); % FIR digital filter design calculation;
        end
    
    case 'hpfilter' % high-pass filter
        switch filter_algo
          case 'but'
            if isempty(N)
              N = 6;
            end
            [B, A] = butter(N, max(freq)/Fn, 'high');
          case 'fir'
            if isempty(N)
              N = 25;
            end
            [B, A] = fir1(N, max(freq)/Fn, 'high');
        end  
        
    case 'bpfilter' % bandpass filter
        switch filter_algo
          case 'but'
            if isempty(N)
              N = 4;
            end
            [B, A] = butter(N, [min(freq)/Fn max(freq)/Fn]);
          case 'fir'
            if isempty(N)
              N = 25;
            end
            [B, A] = fir1(N, [min(freq)/Fn max(freq)/Fn]);
        end
end

%% Applies a filter to the data and corrects edge-artifacts for one-pass filtering.
dcGain = sum(B)/sum(A);

[~,Nb] = size(timecourse);

switch filter_pass
  case 'onepass'
	offset = timecourse(:,1);
	timecourse = timecourse - repmat(offset,1,Nb);
    filtered = filter(B, A, timecourse')' + repmat(dcGain*offset, 1, Nb);
  case 'onepass-reverse'
  	offset = timecourse(:,end);
    timecourse  = fliplr(timecourse) - repmat(offset,1,Nb);
    filt = filter(B, A, timecourse')';
    filtered = fliplr(filt) + repmat(dcGain*offset, 1, Nb);
  case 'twopass'
	% filtfilt does the correction for us
    filtered = filtfilt(B, A, timecourse')';
end



