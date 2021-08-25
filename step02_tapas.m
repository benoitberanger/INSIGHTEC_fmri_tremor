clear
clc

load e.mat

outdir = e.getSerie('run_nm').removeEmpty.getPath;
func   = e.getSerie('run_nm').removeEmpty.getVolume('^wutr(f|v)');
mask   = e.getSerie('run_nm').removeEmpty.getExam.getSerie('anat_T1_UNI').getVolume('^rwp[23]').squeeze();

rp = e.getSerie('run_nm').removeEmpty.getRP('rp');


%%

%----------------------------------------------------------------------------------------------------------------------------------------------------
% ALWAYS MANDATORY
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.TR     = 1.000;
par.nSlice = 72;
par.volume = func;

par.outdir = outdir;

par.physio   = 0;
par.noiseROI = 1;
par.rp       = 1;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Physio
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio_RETROICOR        = 1;
par.physio_HRV              = 1;
par.physio_RVT              = 1;
par.physio_logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
par.logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
par.slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.


%----------------------------------------------------------------------------------------------------------------------------------------------------
% noiseROI
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.noiseROI_mask   = mask;
par.noiseROI_volume = func;

par.noiseROI_thresholds   = [0.95 0.95];     % keep voxels with tissu probabilty >= 95%
par.noiseROI_n_voxel_crop = [2 1];           % crop n voxels in each direction, to avoid partial volume
par.noiseROI_n_components = 10;              % keep n PCA componenets


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Realignment Parameters
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.rp_file      = rp;

par.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
par.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'
par.rp_threshold = 1;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Other
%----------------------------------------------------------------------------------------------------------------------------------------------------
par.print_figures = 0; % 0 , 1 , 2 , 3

par.run           = 1;
par.display       = 0;

par.redo = 0;
par.pct  = 0;

job_physio_tapas( par );
