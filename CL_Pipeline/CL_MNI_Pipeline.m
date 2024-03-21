%% Centiloid pipeline in MNI sopace for all possible combinations
%% Using subject-base RR and cortical target ROI and GAAIN RR

clear all
clc
spm_get_defaults;
spm_jobman('initcfg')
pathFolder=uigetdir('' , 'Choose PET Folder');d = dir(pathFolder);isub = [d(:).isdir]; %# returns logical vector
dirlist = {d(isub).name}';
dirlist(1:2) = [];
MpathFolder=uigetdir('','Chosse MRI Folder');md = dir(MpathFolder);misub = [md(:).isdir]; %# returns logical vector
mdirlist = {md(misub).name}';
mdirlist(1:2) = [];
[AtlasName,Path]=uigetfile('*.nii', 'Choose hammer atlas');
Atlaspath=[Path,AtlasName];
[CLcompositeName,CL_compositePath]=uigetfile('*.nii', 'Choose CL composite');
CLTargetpath=[CL_compositePath,CLcompositeName];
[GAAIN_Name,GAAIN_Path]=uigetfile('*.nii', 'Choose GAAIN cotical composite');
GAAINTargetpath=[GAAIN_Path,GAAIN_Name];
for kk=1:length(dirlist)
PETfile=dir([pathFolder '\' dirlist{kk} '\*Static' ,'*.nii']);
PETpath=[pathFolder '\' dirlist{kk} '\' PETfile(1).name];
    st.vol = spm_vol(PETpath);
    vs = st.vol.mat\eye(4);
    vs(1:3,4) = (st.vol.dim+1)/2;
    spm_get_space(st.vol.fname,inv(vs));
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(spm('dir'),'toolbox','DARTEL','icbm152.nii,1')};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {PETpath};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run',matlabbatch);
    clear matlabbatch
copyfile(Atlaspath,[pathFolder '\' dirlist{kk}])
copyfile(CLTargetpath,[pathFolder '\' dirlist{kk}])
%% RR
Sub_atlas=[pathFolder '\' dirlist{kk} '\Hammer_atlas.nii'];
Sub_Targetp=[pathFolder '\' dirlist{kk} '\CL_Composite.nii'];
%% 
MRfile=dir([MpathFolder '\' mdirlist{kk} '\' ,'*.nii']);%Reading MRI
if ~isempty(MRfile)
MRpath=[MpathFolder '\' mdirlist{kk} '\' MRfile(1).name];
invDffile=dir([MpathFolder '\' mdirlist{kk} '\iy' ,'*.nii']);%Reading MRI
Dfpath=[MpathFolder '\' mdirlist{kk} '\' invDffile.name];
%% 
atlasvol={[Sub_atlas,',',int2str(1)]};
Targetvol={[Sub_Targetp,',',int2str(1)]};
%% 
matlabbatch{1}.spm.spatial.normalise.write.subj.def = ({Dfpath});
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [atlasvol;Targetvol];
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
spm_jobman('run', matlabbatch);
clear matlabbatch
delete(Sub_atlas);
delete(Sub_Targetp);
%% Norm files
files=dir([pathFolder '\' dirlist{kk} '\w' ,'*.nii']);
for nn=1:length(files)
    s{nn}=[files(nn).folder '\' files(nn).name  ,',',int2str(1)];
    path{nn}=[files(nn).folder '\' files(nn).name ];
end
GMpath=[MpathFolder '\' mdirlist{kk} '\c1' MRfile(1).name];
%% Coregister PET with MRI
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[MRpath,',',int2str(1)]};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[PETpath ,',',int2str(1)]};
matlabbatch{1}.spm.spatial.coreg.estimate.other = s';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
clear matlabbatch
%% Reslice
flags = struct('mask',1,'interp',0,'mean',0,'which',1,'wrap',[0 0 0]);
Image{1} = MRpath;
for nn=2:length(files)+1
Image{nn} = [files(nn-1).folder '\' files(nn-1).name ];
end
P = {Image{1:3}}';
spm_reslice(P,flags)
clear Image
files2=dir([pathFolder '\' dirlist{kk} '\w' ,'*.nii']);
for nn=1:length(files2)
    delete([files2(nn).folder '\' files2(nn).name]);
end
%% Defining RR and CL composite mask
rNSub_Targetp=[pathFolder '\' dirlist{kk} '\rwCL_Composite.nii'];
rNSub_atlas=[pathFolder '\' dirlist{kk} '\rwHammer_atlas.nii'];
%% Binarize GM and Creating subject specific Composite
GMpath=[MpathFolder '\' mdirlist{kk} '\c1' MRfile(1).name];
GM=spm_read_vols(spm_vol((GMpath)));GM(isnan(GM))=0;
GM_nii=load_untouch_nii(GMpath);
WMfile=dir([MpathFolder '\' mdirlist{kk} '\c2' ,'*.nii']);%Reading MRI
WMpath=[MpathFolder '\' mdirlist{kk} '\' WMfile.name];
CSFfile=dir([MpathFolder '\' mdirlist{kk} '\c3' ,'*.nii']);%Reading MRI
CSFpath=[MpathFolder '\' mdirlist{kk} '\' CSFfile.name];
WM=spm_read_vols(spm_vol((WMpath)));WM(isnan(WM))=0;
CSF=spm_read_vols(spm_vol((CSFpath)));CSF(isnan(CSF))=0;
GMmask=zeros(size(GM));
GMmask(GM>WM & GM>CSF)=1;
WMmask=zeros(size(GM));
WMmask(WM>GM & WM>CSF)=1;
brain=WMmask+GMmask;
hammer=double(spm_read_vols(spm_vol(rNSub_atlas)));hammer(isnan(hammer))=0;
hammer(brain<1)=0;
hammer (hammer==17 & GMmask==1)=84;
hammer (hammer==18 & GMmask==1)=84;
VOI=double(spm_read_vols(spm_vol(rNSub_Targetp)));VOI(isnan(VOI))=0;
VOI(GMmask<1)=0;
nii=make_nii(VOI);
nii.hdr.hist = GM_nii.hdr.hist;
nii.hdr.dime.pixdim = GM_nii.hdr.dime.pixdim;
newname = [ 'N_CL_Atlas_' dirlist{kk} ,'.nii'];
save_nii(nii,newname);
movefile(newname,dirlist{kk})
clear nii
nii=make_nii(hammer);
nii.hdr.hist = GM_nii.hdr.hist;
nii.hdr.dime.pixdim = GM_nii.hdr.dime.pixdim;
hnewname = [ 'N_hammer_' dirlist{kk} ,'.nii'];
save_nii(nii,hnewname);
movefile(hnewname,dirlist{kk}) 
clear nii
delete(rNSub_Targetp);
delete(rNSub_atlas);
%% Move atlases to MNI space
Targetp=[pathFolder '\' dirlist{kk} '\N_CL_Atlas_' dirlist{kk} '.nii'];
atlas=[pathFolder '\' dirlist{kk} '\N_hammer_' dirlist{kk} '.nii'];
atlasvol={[atlas,',',int2str(1)]};
targetvol={[Targetp,',',int2str(1)]};
PETvol={[PETpath,',',int2str(1)]};
MRvol={[MRpath,',',int2str(1)]};
Dffile=dir([MpathFolder '\' mdirlist{kk} '\y_' ,'*.nii']);%Reading MRI
Dfpath=[MpathFolder '\' mdirlist{kk} '\' Dffile.name];
%% Moving the cortical ROIs to the MNI space
matlabbatch{1}.spm.spatial.normalise.write.subj.def = ({Dfpath});
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [MRvol;PETvol;atlasvol;targetvol];
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
spm_jobman('run', matlabbatch);
clear matlabbatch
delete(atlas);delete(Targetp);
fprintf(1,['subject ' num2str(kk) '  done\n'])
end
end

