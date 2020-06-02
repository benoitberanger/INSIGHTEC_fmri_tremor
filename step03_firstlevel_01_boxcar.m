clear
clc

load e.mat


%% Prepare paths and regexp

par.display = 0;
par.run     = 1;
par.pct     = 0;
par.verbose = 2;


%% dirs & files

model_name = {'classic','tapas'};

dirStats_1 = e.mkdir('model',model_name{1});
dirStats_2 = e.mkdir('model',model_name{2});

dirFonc = e.getSerie('run_nm').toJob;
e.getSerie('run_nm').addStim('onsets','.mat','stim',1)
onsetFile = e.getSerie('run_nm').getStim('stim').toJob;

par.rp       = 1;
par.file_reg = '^sw.*nii';


%% Specify boxcar

par.redo    = 0;
par.sge     = 0;
par.run     = 0;
par.display = 0;

% Boxcar + 6rp (classic)
par.rp_regex = '^rp.*txt';
job1 = job_first_level_specify(dirFonc,dirStats_1,onsetFile,par);

% Boxcar + TAPAS
par.rp_regex = 'multiple_regressors.txt';
job2 = job_first_level_specify(dirFonc,dirStats_2,onsetFile,par);

jobs = [job1 job2 ];

par.sge     = 1;
par.run     = 0;
par.display = 0;
par.jobname  = 'spm_glm_def';

job_ending_rountines(jobs,[],par);

return

%% Estimate

e.addModel('model','^classic$','classic');
e.addModel('model','^tapas$','tapas');

save('e','e')

fspm = e.getModel(cellstr2regex({'classic','tapas'},1)).removeEmpty.toJob;

clear par
par.sge     = 1;
par.run     = 0;
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
par.sge     = 1;
par.run     = 0;
par.display = 0;
par.jobname  = 'spm_glm_con';

par.sessrep = 'none';

par.delete_previous = 1;
par.report=0;
fspm = e.getModel(cellstr2regex({'classic','tapas'},1)).removeEmpty.toJob;
job_first_level_contrast(fspm,contrast,par);


