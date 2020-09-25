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
par.run = 1;
par.pct = 0;
par.sge = 0;


%% Get files paths

% e = exam(maindir,'nifti',{'INISGHTEC','INSIGHTEC','ULTRABRAIN'});
e = exam(maindir,'nifti','ULTRABRAIN');

% T1
e.addSerie('INV1$'      ,'anat_T1_INV1',1)
e.addSerie('INV2$'      ,'anat_T1_INV2',1)
e.addSerie('UNI_Images$','anat_T1_UNI' ,1)
e.getSerie('anat').addVolume('^s.*nii','s',1)

% Run : Exec
e.addSerie('RS$'        , 'rs_nm', 1)
% e.addSerie('RS_refBLIP$', 'rs_bp', 1) % refAP

e.getSerie('rs').addVolume('^f.*nii','f',3)
e.getSerie().addJson('^dic','j')

% Unzip if necessary (with PCT ?)

e.reorderSeries('name'); % mostly useful for topup, that requires pairs of (AP,PA)/(PA,AP) scans

e.explore

%% denoise mp2rage

% done graphically

e.getSerie('anat_T1_UNI').addVolume('^cs.*nii','cs',1)


%% Segment anat with cat12

clear par
par.redo= 0;
par.run = 1;
par.pct = 0;
par.sge = 0;

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

anat = e.gser('anat_T1_UNI').gvol('^cs');
job_do_segmentCAT12(anat,par);

% temporary, until the next PR :
e.getSerie('anat_T1_UNI').addVolume('^p0.*nii','p0',1)


%% Sort echos

clear par
par.run  = 1;
par.fake = 0;
par.sge  = 0;

par.redo = 0;
meinfo = job_sort_echos( e.getSerie('rs') , par );


%% job_afni_proc_multi_echo

clear par
par.run  = 1;
par.fake = 0;
par.sge  = 0;
par.redo = 0;

par.seperate = 1;
par.write_nifti = 1;

par.subdir = 'afni_vtd';
par.blocks  = {'despike', 'tshift', 'volreg'};
job_afni_proc_multi_echo( meinfo, par );


%% do_fsl_robust_mask_epi

fin  = e.getSerie('rs').getVolume('^vtde1');

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.redo  = 0;
[ job, fmask ] = do_fsl_robust_mask_epi( fin, par );


%% tedana

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.redo  = 0;
par.pct   = 0;

job_tedana_009a1( meinfo, 'vtd', 'tedana_vtd', 'bet_Tmean_vtde1_mask.nii.gz ', par );

e.getSerie('rs').getVolume({'tedana','bet'}).unzip_and_keep(par);


%% Coregister

clear par
par.run   = 1;
par.sge   = 0;
par.redo  = 0;
par.type  = 'estimate';

src = e.getSerie('rs').removeEmpty.getVolume('^bet_Tmean_vtde1$');
oth = e.getSerie('rs').removeEmpty.getVolume('^tedana.*ts_OC$');
tmp_exam = [e.getSerie('rs').removeEmpty.exam];
ref = tmp_exam.getSerie('anat_T1').getVolume('^p0');

% Save checkpoint
to_zip_and_keep = src(:)+oth(:);
to_zip_and_keep.zip_and_keep(par);

job_coregister(src,ref,oth,par);


%% Normalize

clear par
par.run  = 1;
par.sge  = 0;
par.redo = 0;
par.vox = [2.5 2.5 2.5]; % IMPORTANT keep original EPI voxel size
img = e.getSerie('rs').getVolume('^tedana.*ts_OC$').removeEmpty;
tmp_exam = [img.exam];
y   = tmp_exam.getSerie('anat_T1_UNI').getVolume('^y');

par.auto_add_obj = 0;
job_apply_normalize(y,img,par);

% auto_add not yet proprely coded (17/02/2020)
e.getSerie('rs').addVolume('tedana',   '^wts_OC.nii$',   'wts_OC',1)
e.getSerie('rs').addVolume('tedana','^wdn_ts_OC.nii$','wdn_ts_OC',1)
