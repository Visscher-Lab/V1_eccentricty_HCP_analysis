%sublist for paired
%central and far
clear
clc
outpdir='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis';
mkdir (outpdir)
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['*']);
cfiles(1) = [];
cfiles(1) = [];
cfiles(787:end) = [];
cfiles = rmfield(cfiles, {'folder', 'date', 'bytes', 'isdir', 'datenum'});
u=repelem(cfiles,2)
[u(1:1572).unnamed] = deal([]); 
[u(1:2:end).unnamed] = deal('.central');
[u(2:2:end).unnamed] = deal('.far');
cd(outpdir)
writetable(struct2table(u), 'sublistcf.txt')
% delete first row and remove commas

%sublist mid and central
clear
clc
outpdir='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis'; 
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['*']);
cfiles(1) = [];
cfiles(1) = [];
cfiles(787:end) = [];
cfiles = rmfield(cfiles, {'folder', 'date', 'bytes', 'isdir', 'datenum'});
u=repelem(cfiles,2)
[u(1:1572).unnamed] = deal([]); 
[u(1:2:end).unnamed] = deal('.mid');
[u(2:2:end).unnamed] = deal('.central');
cd(outpdir)
writetable(struct2table(u), 'sublistmc.txt')
%then copy from reconall folder, delete first row and remove commas

%sublist mid and far
clear
clc
outpdir='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis'; 
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['*']);
cfiles(1) = [];
cfiles(1) = [];
cfiles(787:end) = [];
cfiles = rmfield(cfiles, {'folder', 'date', 'bytes', 'isdir', 'datenum'});
u=repelem(cfiles,2)
[u(1:1572).unnamed] = deal([]); 
[u(1:2:end).unnamed] = deal('.mid');
[u(2:2:end).unnamed] = deal('.far');
cd(outpdir)
writetable(struct2table(u), 'sublistmf.txt')
%then copy from reconall folder, delete first row and remove commas

%sublist f v c
clear
clc
outpdir='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis'; 
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['*']);
cfiles(1) = [];
cfiles(1) = [];
cfiles(787:end) = [];
cfiles = rmfield(cfiles, {'folder', 'date', 'bytes', 'isdir', 'datenum'});
u=repelem(cfiles,2)
[u(1:1572).unnamed] = deal([]); 
[u(1:2:end).unnamed] = deal('.far');
[u(2:2:end).unnamed] = deal('.central');
cd(outpdir)
writetable(struct2table(u), 'sublistfc.txt')
%then copy from reconall folder, delete first row and remove commas

%fsgd for paired
clear
clc
outpdir='/data/project/vislab/a/HCP_diff_Sara/V1stoFNtprob/surface_analysis'; 
parentdir = '/data/project/vislab/a/HCP_diff_Sara/subjects_reconall_postcheck';
cd(parentdir);
cfiles = dir(['*']);
cfiles(1) = [];
cfiles(1) = [];
cfiles(787:end) = [];
[cfiles(1:786).folder] = deal([]); 
[cfiles(1:786).date] = deal([]); 
[cfiles(1:786).folder] = deal(['paired']);
[cfiles(1:786).date] = deal(['Main']);
[cfiles(1:786).unnamed] = deal(['Input']);
cfiles = rmfield(cfiles, {'bytes', 'isdir', 'datenum'});
cd(outpdir)
writetable(struct2table(cfiles), 'FSGD.txt')
%delete first row and remove commas and moved paired over a space
% so that the lines look like this: Input 100206paired Main