clear all;
% Set up Freesurfer
addpath(genpath('/share/apps/rc/software/FreeSurfer/6.0.0-centos6_x86_64/matlab/'))
% Set up FSL environment
setenv( 'FSLDIR', '/share/apps/rc/software/FSL/5.0.9-centos6_64');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;
%Necessary paths
addpath('/data/project/vislab/a/HCP_diff_Sara/pinar')
outputdir = '/data/project/vislab/a/HCP_diff_Sara/UP_funcanalysis/data'; % directory where scans are
datadir = '/data/project/vislab/a/HCP_diff_Sara/FNlabelvol';
parentdir= '/data/project/vislab/raw/HCP_900sub/HCP_func';
cd(outputdir);
subs = dir('*'); % gets a list of all the subjects
subs = subs(3:length(subs)-3);

%Open text files
fid=fopen('seedtovoxel_correlations771end.txt','wt');

for r=771:length(subs)
    
    cd([datadir '/' subs(r).name '/']);
    system(['fslmaths ' datadir '/' subs(r).name '/net4_funcROI.nii -bin ' datadir '/' subs(r).name '/net4_funcROI_bin.nii.gz']);
    system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/net4ROI_ts.txt -m ' datadir '/' subs(r).name '/net4_funcROI_bin.nii.gz']);
    system(['fslmaths ' datadir '/' subs(r).name '/net6_funcROI.nii -bin ' datadir '/' subs(r).name '/net6_funcROI_bin.nii.gz']);
    system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/net6ROI_ts.txt -m ' datadir '/' subs(r).name '/net6_funcROI_bin.nii.gz']);
    system(['fslmaths ' datadir '/' subs(r).name '/net7_funcROI.nii -bin ' datadir '/' subs(r).name '/net7_funcROI_bin.nii.gz']);
    system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/net7ROI_ts.txt -m ' datadir '/' subs(r).name '/net7_funcROI_bin.nii.gz']);
    ts1 = 'centralROI_ts.txt';
    ts2 = 'midROI_ts.txt';
    ts3 = 'farROI_ts.txt';
    ts4 = 'net4ROI_ts.txt';
    ts5 = 'net6ROI_ts.txt';
    ts6 = 'net7ROI_ts.txt';
    
    % Load seed ROI timeseries:
    T1 = load(ts1);
    T2 = load(ts2);
    T3 = load(ts3);
    T4 = load(ts4);
    T5 = load(ts5);
    T6 = load(ts6);
    
    % Perform correlation:
    out(:,1) = corr(T1,T4);
    out(:,2) = corr(T1,T5);
    out(:,3) = corr(T1,T6);
    out(:,4) = corr(T2,T4);
    out(:,5) = corr(T2,T5);
    out(:,6) = corr(T2,T6);
    out(:,7) = corr(T3,T4);
    out(:,8) = corr(T3,T5);
    out(:,9) = corr(T3,T6);
    out(isnan(out)==1) = 0;
    
    % Perform r to z transform:
    out = 0.5*log((1+out)./(1-out));
    
    % Write results in a text file
    fprintf(fid, '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n',out.');
    
    %clearvars -except parentdir outputdir subs datadir r fid;
end
fclose(fid);
