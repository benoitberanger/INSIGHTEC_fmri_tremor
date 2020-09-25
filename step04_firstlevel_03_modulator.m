clear
clc

load e.mat

%% dirs & files


dirFonc = e.getSerie('run_nm').toJob;
e.getSerie('run_nm').addStim('onsets','.mat','stim',1)
e.getSerie('run_nm').addVolume('mask.nii','mask',1)
stim = e.getSerie('run_nm').getStim('stim').load; stim = stim{1};

par.rp       = 1;
par.file_reg = '^sw.*nii';
TR = 1.000; % seconds


%% Make Posture regressor modular
% We need to modulate the regressor for each value of power (sampled @ TR)

new_onset = [];
for block = 1 : length(stim.onsets{3})
    next = stim.onsets{3}(block) : TR : stim.onsets{3}(block)+stim.durations{3}(block);
    new_onset = [new_onset next]; %#ok<AGROW>
end
new_duration = TR*ones(size(new_onset));


%% Specify user regressors

list = {
    'modulator_FCR_L'
    'modulator_FCR_R'
    'modulator_ECR_L'
    'modulator_ECR_R'
    'modulator_ACC_X'
    'modulator_ACC_Y'
    'modulator_ACC_Z'
    };

dir_regressor = gdir(e.gpath,'electrophy');

onsets = cell(0);

for l = 1 : length(list)
    
    e.getSerie('run_nm').addStim('electrophy',['^' list{l} '.mat'],list{l},1)
    emg = e.getSerie('run_nm').getStim(list{l}).load; emg = emg{1};
    
%     hold on
%     plot(emg.R(:,1))
    
    bin = false(1,size(emg.R,1));
    bin(round(new_onset)) = 1;
    
    onsets{l}(1).name     = stim.names    {1};
    onsets{l}(1).onset    = stim.onsets   {1};
    onsets{l}(1).duration = stim.durations{1};
    onsets{l}(1).tmod     = 0;
    onsets{l}(1).orth     = 1;
    
    onsets{l}(2).name     = stim.names    {2};
    onsets{l}(2).onset    = stim.onsets   {2};
    onsets{l}(2).duration = stim.durations{2};
    onsets{l}(2).tmod     = 0;
    onsets{l}(2).orth     = 1;
    
    onsets{l}(3).name     = stim.names    {3};
    onsets{l}(3).onset    = new_onset;
    onsets{l}(3).duration = new_duration;
    onsets{l}(3).tmod     = 0;
    onsets{l}(3).orth     = 1;
    
    onsets{l}(3).pmod(1).name  = 'reg';
    onsets{l}(3).pmod(1).param = emg.R(bin,1);
    onsets{l}(3).pmod(1).poly  = 1;
    onsets{l}(3).pmod(2).name  = 'dreg';
    onsets{l}(3).pmod(2).param = emg.R(bin,2);
    onsets{l}(3).pmod(2).poly  = 1;
    
    onsets{l}(4).name     = stim.names    {4};
    onsets{l}(4).onset    = stim.onsets   {4};
    onsets{l}(4).duration = stim.durations{4};
    onsets{l}(4).tmod     = 0;
    onsets{l}(4).orth     = 1;
    
end

% legend(list)

% return


%% Specify

par.sge     = 1;
par.run     = 0;
par.display = 0;
par.mask_thr = 0.1;
par.jobname  = 'spm_glm_def';
par.rp_regex = 'multiple_regressors.txt';
mask = e.getSerie('run_nm').getVolume('mask').toJob;
par.mask = repmat(mask,size(list));
dirStats = r_mkdir( fullfile( e.getPath, 'model') ,  list );

job_first_level_specify( repmat(dirFonc,size(list)), dirStats, onsets,par);

return


%% Estimation

for l = 1 : length(list)
    e.addModel('model',list{l},list{l});
end

save('e','e')

fspm = e.getModel('modulator').removeEmpty.toJob;

clear par
par.sge     = 1;
par.run     = 0;
par.display = 0;
job_first_level_estimate(fspm,par);


%% Contrast

Instructions        = [1 0 0 0 0 0];
Relax               = [0 1 0 0 0 0];
Posture             = [0 0 1 0 0 0];
Posture__power      = [0 0 0 1 0 0];
Posture__power_diff = [0 0 0 0 1 0];
EndText             = [0 0 0 0 0 1];

contrast_T.values = {
    
Instructions
Relax
Posture
Posture__power
Posture__power_diff
EndText
Posture - Relax

}';

contrast_T.names = {
    
'Instructions'
'Relax'
'Posture'
'Posture__power'
'Posture__power_diff'
'EndText'

'Posture - Relax'

}';

contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
contrast_F.values = {eye(6)};
contrast_F.types  = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

contrast.names  = [contrast_F.names  contrast_T.names ];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types ];


%% Contrast : write

clear par
par.sessrep = 'none';

par.sge     = 1;
par.run     = 0;
par.display = 0;
par.jobname  = 'spm_glm_con';

par.delete_previous = 1;
par.report          = 0;
job = job_first_level_contrast(fspm,contrast,par);

