clear
clc

load e.mat


%% Prepare paths and regexp

par.display = 0;
par.run     = 0;
par.pct     = 0;
par.verbose = 2;


%% dirs & files

dirFonc = e.getSerie('run_nm').toJob;
e.getSerie('run_nm').addStim('onsets','.mat','stim',1)
onsetFile = e.getSerie('run_nm').getStim('stim').toJob;
e.getSerie('run_nm').addVolume('mask.nii','mask',1)

par.rp       = 1;
par.file_reg = '^sw.*nii';
par.mask_thr = 0.1;


%% Specify user regressors

list = {
    'FCR_L'
    'FCR_R'
    'ECR_L'
    'ECR_R'
    'EMG'
    'ACC_X'
    'ACC_Y'
    'ACC_Z'
    'ACC_XYZ_mean'
    'ACC_XYZ_pca'
    
    'power_FCR_L'
    'power_FCR_R'
    'power_ECR_L'
    'power_ECR_R'
    'power_ACC_X'
    'power_ACC_Y'
    'power_ACC_Z'
    
    'logpower_FCR_L'
    'logpower_FCR_R'
    'logpower_ECR_L'
    'logpower_ECR_R'
    'logpower_ACC_X'
    'logpower_ACC_Y'
    'logpower_ACC_Z'
    
    };

dir_regressor = gdir(e.gpath,'electrophy');

for l = 1 : length(list)
    e.getSerie('run_nm').addStim('electrophy',['^' list{l} '.mat'],list{l},1)
end

jobs = cell(0);
j = 0;
for l = 1 : length(list)
    j = j + 1;
    
    par.mask = e.getSerie('run_nm').getVolume('mask').toJob;
    
    % FARM + 6rp (classic)
    par.rp_regex = '^rp.*txt';
    par.file_regressor = e.getSerie('run_nm').getStim(list{l}).toJob(1);
    dirStats = e.mkdir('model',['classic_' list{l}]);
    jobtmp = job_first_level_specify(dirFonc,dirStats,repmat(onsetFile,size(e)),par);
    jobs = [jobs jobtmp];
    
    % FARM + TAPAS
    par.rp_regex = 'multiple_regressors.txt';
    par.file_regressor = e.getSerie('run_nm').getStim(list{l}).toJob(1);
    dirStats = e.mkdir('model',['tapas_' list{l}]);
    jobtmp = job_first_level_specify(dirFonc,dirStats,repmat(onsetFile,size(e)),par);
    jobs = [jobs jobtmp];
    
end

par.sge     = 1;
par.run     = 0;
par.display = 0;
par.jobname  = 'spm_glm_def';

job_ending_rountines(jobs,[],par);

return

%% Estimate

for l = 1 : length(list)
    e.addModel('model',['classic_' list{l} '$'],['classic_' list{l}]);
    e.addModel('model',[  'tapas_' list{l} '$'],[  'tapas_' list{l}]);
end

save('e','e')

fspm = e.getModel(cellstr2regex({'classic_','tapas_'})).removeEmpty.toJob;

clear par
par.sge     = 1;
par.run     = 0;
par.display = 0;
job_first_level_estimate(fspm,par);


%% Contrast

Instructions = [1 0 0 0   0 0];
Relax        = [0 1 0 0   0 0];
Posture      = [0 0 1 0   0 0];
EndText      = [0 0 0 1   0 0];
reg          = [0 0 0 0   1 0];
dreg         = [0 0 0 0   0 1];

contrast_T.values = {
    
Instructions
Relax
Posture
EndText
reg
dreg

Posture-Relax

}';

contrast_T.names = {
    
'Instructions'
'Relax'
'Posture'
'EndText'
'reg'
'dreg'

'Posture-Relax'

}';

% reg  = [1 0];
% dreg = [0 1];
% 
% contrast_T.values = {
%     
% reg
% dreg
% 
% }';
% 
% contrast_T.names = {
%     
% 'reg'
% 'dreg'
% 
% }';

contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
% contrast_F.values = {eye(2)};
contrast_F.values = {eye(6)};
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
job_first_level_contrast(fspm,contrast,par);

