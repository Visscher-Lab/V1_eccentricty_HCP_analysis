clear all
addpath('/data/project/vislab/a/HCP_diff_Sara/pinar')
outputdir = '/data/project/vislab/a/HCP_diff_Sara/UP_funcanalysis/data'; % directory where scans are
datadir = '/data/project/vislab/a/HCP_diff_Sara/FNlabelvol';
parentdir= '/data/project/vislab/raw/HCP_900sub/HCP_func';
cd(parentdir);
subs = dir('*'); % gets a list of all the subjects
subs = subs(3:length(subs));

for r=1:length(subs)
    if exist ([parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/' subs(r).name '.nii'], 'file' ) ==2
        system(['fslmaths ' datadir '/' subs(r).name '/newcentral_funcROI.nii -bin ' datadir '/' subs(r).name '/central_funcROI_bin.nii.gz']);
        system(['fslmaths ' outputdir '/' subs(r).name '/Concatenated_RestingState.nii.gz -mas ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST2_RL/' subs(r).name '.nii ' datadir '/' subs(r).name '/concatenated_masked.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/centralROI_ts.txt -m ' datadir '/' subs(r).name '/central_funcROI_bin.nii.gz']);
        cd([datadir '/' subs(r).name '/']);
        system(['gzip -d concatenated_masked.nii.gz']);
        ts = 'centralROI_ts.txt';
        D='concatenated_masked.nii';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/central']);
        clear ts;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/central']);
        
        system(['fslmaths ' datadir '/' subs(r).name '/mid_funcROI.nii -bin ' datadir '/' subs(r).name '/mid_funcROI_bin.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/midROI_ts.txt -m ' datadir '/' subs(r).name '/mid_funcROI_bin.nii.gz']);
        ts = 'midROI_ts.txt';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/mid']);
        clear ts;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/mid']);
        
        system(['fslmaths ' datadir '/' subs(r).name '/newfar_funcROI.nii -bin ' datadir '/' subs(r).name '/far_funcROI_bin.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/farROI_ts.txt -m ' datadir '/' subs(r).name '/far_funcROI_bin.nii.gz']);
        ts = 'farROI_ts.txt';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/far']);
        clear ts D;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/far']);
    end
    
    
    if exist ([parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/' subs(r).name '.nii'], 'file' ) ==2
        system(['fslmaths ' datadir '/' subs(r).name '/newcentral_funcROI.nii -bin ' datadir '/' subs(r).name '/central_funcROI_bin.nii.gz']);
        system(['fslmaths ' outputdir  '/' subs(r).name '/Concatenated_RestingState.nii.gz -mas ' parentdir '/' subs(r).name '/MNINonLinear/Results/rfMRI_REST1_RL/' subs(r).name '.nii ' datadir '/' subs(r).name '/concatenated_masked.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/centralROI_ts.txt -m ' datadir '/' subs(r).name '/central_funcROI_bin.nii.gz']);
        cd([datadir '/' subs(r).name]);
        system(['gzip -d concatenated_masked.nii.gz']);
        ts = 'centralROI_ts.txt';
        D = 'concatenated_masked.nii';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/central']);
        clear ts;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/central']);
        
        system(['fslmaths ' datadir '/' subs(r).name '/mid_funcROI.nii -bin ' datadir '/' subs(r).name '/mid_funcROI_bin.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/midROI_ts.txt -m ' datadir '/' subs(r).name '/mid_funcROI_bin.nii.gz']);
        ts = 'midROI_ts.txt';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/mid']);
        clear ts;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/mid']);
        
        system(['fslmaths ' datadir '/' subs(r).name '/newfar_funcROI.nii -bin ' datadir '/' subs(r).name '/far_funcROI_bin.nii.gz']);
        system(['fslmeants -i ' datadir '/' subs(r).name '/concatenated_masked.nii.gz -o ' datadir '/' subs(r).name '/farROI_ts.txt -m ' datadir '/' subs(r).name '/far_funcROI_bin.nii.gz']);
        mkdir([outputdir '/' subs(r).name '/far']);
        ts = 'farROI_ts.txt';
        SCA(ts,D);
        mkdir([outputdir '/' subs(r).name '/far']);
        clear ts D;
        movefile('SCA_result.nii.gz',[outputdir '/' subs(r).name '/far']);
    end
    clearvars -except parentdir outputdir subs datadir r;
end
