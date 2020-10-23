clear
clc

%% Prepare paths and regexp

maindir = '/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor';

% try
% do_delete(gdir(pwd,'spm_'),0)
% do_delete(fullfile(pwd,'do_workflow_qsub.sh'),0)
% catch err
% end

par.redo= 0;
par.run = 0;
par.pct = 0;
par.sge = 1;


%% Get files paths

% e = exam(maindir,'nifti',{'INISGHTEC','INSIGHTEC','ULTRABRAIN'});
e = exam(maindir,'nifti','ULTRABRAIN'); 


% T1
e.addSerie('INV1$'      ,'anat_T1_INV1',1)
e.addSerie('INV2$'      ,'anat_T1_INV2',1)
e.addSerie('UNI_Images$','anat_T1_UNI' ,1)
e.getExam( '002_FJ_INCLUSION'   ).getSerie('anat').addVolume('^s.*nii','s')
e.getExam('0003_PP_V1_INCLUSION').getSerie('anat').addVolume('^v.*nii','v')
e.getExam('0002_FJ_V4'          ).getSerie('anat').addVolume('^v.*nii','v')

% Run : Exec
e.getExam( '002_FJ_INCLUSION'   ).addSerie('fMRI_tremor$'                 , 'run_nm', 1)
e.getExam( '002_FJ_INCLUSION'   ).addSerie('fMRI_tremor_refBLIP$'         , 'run_bp', 1) % refAP
e.getExam( '002_FJ_INCLUSION'   ).getSerie('run').addVolume('^f.*nii','f',1)

e.getExam('0003_PP_V1_INCLUSION').addSerie('fMRI_tremor_run\d{2}$'        , 'run_nm',2)
e.getExam('0003_PP_V1_INCLUSION').addSerie('fMRI_tremor_run\d{2}_refBLIP$', 'run_bp',2) % refAP
e.getExam('0003_PP_V1_INCLUSION').getSerie('run').addVolume('^v.*nii','v',1)

e.getExam('0002_FJ_V4'          ).addSerie('fMRI_tremor$'                 , 'run_nm', 1)
e.getExam('0002_FJ_V4'          ).addSerie('fMRI_tremor_refBLIP$'         , 'run_bp', 1) % refAP
e.getExam('0002_FJ_V4'          ).getSerie('run').addVolume('^v.*nii','v',1)

% Unzip if necessary (with PCT ?)
e.unzipVolume(par);

e.reorderSeries('name'); % mostly useful for topup, that requires pairs of (AP,PA)/(PA,AP) scans

e.explore


%% denoise mp2rage

INV1 = e.getSerie('anat_T1_INV1').getVolume('^v|^s');
INV2 = e.getSerie('anat_T1_INV2').getVolume('^v|^s');
UNI  = e.getSerie('anat_T1_UNI' ).getVolume('^v|^s');

par.regularization = 50;
par.prefix         = 'c';
job_denoise_mp2rage(INV1,INV2,UNI,par);
par = rmfield(par,'prefix');


%% Segment anat with cat12

par.subfolder = 0;         % 0 means "do not write in subfolder"
par.biasstr   = 0.5;
par.accstr    = 0.5;
par.WM        = [1 0 1 0]; %                          (wp2*)     /                        (mwp2*)     /              (p2*)     /                            (rp2*)
par.CSF       = [1 0 1 0]; %                          (wp3*)     /                        (mwp3*)     /              (p3*)     /                            (rp3*)
par.TPMC      = [1 0 1 0]; %                          (wp[456]*) /                        (mwp[456]*) /              (p[456]*) /                            (rp[456]*)
par.label     = [1 0 0] ;  % native (p0*)  / normalize (wp0*)  / dartel (rp0*)       This will create a label map : p0 = (1 x p1) + (3 x p2) + (1 x p3)
par.bias      = [1 1 0] ;  % native (ms*)  / normalize (wms*)  / dartel (rms*)       This will save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % native (mis*) / normalize (wmis*) / dartel (rmis*)      This will save the bias field corrected  + SANLM (local) T1
par.warp      = [1 1];     % Warp fields  : native->template (y_*) / native<-template (iy_*)
par.doSurface = 0;
par.doROI     = 0;         % Will compute the volume in each atlas region
par.jacobian  = 0;         % Write jacobian determinant in normalize space

anat = e.gser('anat_T1_UNI').gvol('^c(s|v)');
job_do_segmentCAT12(anat,par);

par.jobname = 'zipWMCSF';
e.gser('anat_T1_UNI').gvol('^wp[23]').zip_and_keep(par);
par = rmfield(par,'jobname');


%% Preprocess fMRI runs

%realign and reslice
par.type = 'estimate_and_reslice';
ffunc_nm = e.getSerie('run_nm').getVolume('^(f|v)');
j_realign_reslice_nm = job_realign(ffunc_nm,par);

%realign and reslice opposite phase
par.type = 'estimate_and_reslice';
ffunc_bp = e.getSerie('run_bp').getVolume('^(f|v)');
j_realign_reslice_op = job_realign(ffunc_bp,par);

%topup and unwarp
ffunc_all = e.getSerie('run').getVolume('^r(f|v)');
do_topup_unwarp_4D(ffunc_all,par);
par.jobname = 'zipUT';
e.gser('run').gvol('^ut').zip_and_keep(par);
par = rmfield(par,'jobname');

%coregister mean fonc on brain_anat
fanat = e.getSerie('anat_T1_UNI').getVolume('^p0');
fmean = e.getSerie('run_nm').getVolume('^utmean(f|v)'); fmean = fmean(:,1); % use the mean of the run1 to estimate the coreg
fo    = e.getSerie('run_nm').getVolume('^utr(f|v)');
par.type = 'estimate';
par.jobname = 'spm_coreg_epi2anat';
j_coregister=job_coregister(fmean,fanat,fo,par);
par = rmfield(par,'jobname');

%apply normalize
par.vox    = [2.01923  2.01923  2];
to_warp    = fo.removeEmpty + fmean.removeEmpty;
warp_field = to_warp.getExam.getSerie('anat_T1_UNI').getVolume('^y');
j_apply_normalize=job_apply_normalize(warp_field,to_warp,par);

%smooth the data
ffonc = e.getSerie('run_nm').getVolume('wutr(f|v)').removeEmpty;
par.smooth = [4 4 4];
j_smooth=job_smooth(ffonc,par);

% coregister WM & CSF on functionnal (using the warped mean)
if isfield(par,'prefix'), par = rmfield(par,'prefix'); end
ref = e.getSerie('run_nm');
ref = ref(:,1).getVolume('^wutmean(f|v)'); % first acquired run (time)
src = e.getSerie('anat_T1_UNI').getVolume('^wp2');
oth = e.getSerie('anat_T1_UNI').getVolume('^wp3');
par.type = 'estimate_and_write';
par.jobname = 'spm_coreg_wmcsf2epi';
job_coregister(src,ref,oth,par);
par = rmfield(par,'jobname');

% robust EPI mask for firstlevel
epi = e.getSerie('run_nm').getVolume('^wutr(f|v)').removeEmpty;
par.fsl_output_format = 'NIFTI';
do_fsl_robust_mask_epi(epi,par);

save('e','e')

