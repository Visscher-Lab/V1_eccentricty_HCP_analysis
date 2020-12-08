%FP central versus far
clear
clc
diary_dir ='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis'; %set location for output of diary
 diary([diary_dir '/preprocscriptlistFPcflh.txt']); %diary location and filename
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['**']);
cfiles(1:2) = [];
cfiles(787:end) = [];
str='';
for ii = 1:length(cfiles)
    sub = cfiles(ii).name;
    flag= [' --isp $in/FP/', sub, '.central/lhsurf.mgh --isp $in/FP/', sub, '.far/lhsurf.mgh'];
    str = [str sprintf(flag)]
end
diary OFF;
%copy and paste the string output in the diary to the mris_preproc script


