clear
clc

load e.mat


%% Cluster ?

CLUSTER = 0;


%% dirs & files

model_name = 'boxcar';

dirStats = e.mkdir('model',model_name);

dirFonc = e.getSerie('run_nm').toJob;

e.getExam( '002_FJ').getSerie('run_nm'    ).addStim('onsets',     '.mat','stim'   ,1)
e.getExam('0003_PP').getSerie('run_nm_001').addStim('onsets','run01.mat','stim_01',1)
e.getExam('0003_PP').getSerie('run_nm_002').addStim('onsets','run02.mat','stim_02',1)

onsetFile = e.getSerie('run_nm').getStim('stim').toJob;
run = e.getSerie('run_nm');
mask = run(:,1).getVolume('mask').toJob; % 1 mask per model, whatever the number of runs

par.rp       = 1;
par.file_reg = '^sw.*nii';
par.TR       = 1.000;


%% Specify

if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo    = 0;
par.display = 0;
par.jobname  = 'spm_glm_def';
par.mask_thr = 0.1;

% Boxcar + TAPAS
par.rp_regex = 'multiple_regressors.txt';
job_first_level_specify(dirFonc,dirStats,onsetFile,par);


%% Estimate

e.addModel('model','^boxcar$','boxcar');

save('e','e')

fspm = e.getModel('boxcar').removeEmpty.toJob;

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.display = 0;
job_first_level_estimate(fspm,par);


%% Contrast

Instructions = [1 0 0 0];
Relax        = [0 1 0 0];
Posture      = [0 0 1 0];
EndText      = [0 0 0 1];

contrast_T.values = {
    
Instructions
Relax
Posture
EndText

Posture - Relax

}';

contrast_T.names = {
    
'Instructions'
'Relax'
'Posture'
'EndText'

'Posture - Relax'

}';

contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
contrast_F.values = {eye(4)};
contrast_F.types  = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

contrast.names  = [contrast_F.names  contrast_T.names ];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types ];


%% Contrast : write

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.display = 0;
par.jobname  = 'spm_glm_con';

par.sessrep = 'repl';

par.delete_previous = 1;
par.report=0;
job_first_level_contrast(fspm,contrast,par);

