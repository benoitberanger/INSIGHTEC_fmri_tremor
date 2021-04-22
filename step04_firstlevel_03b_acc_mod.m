clear
clc

load e.mat


%% Cluster ?

CLUSTER = 1;


%% dirs & files

volume = e.getSerie('run_nm').getVolume('^sw').toJob(1);

e.getSerie('run_nm').addJson('multiple_regressors.txt','rp')
rp = e.getSerie('run_nm').getJson('rp').toJob(1);

e_2run = e.getExam('0003_PP_V1_INCLUSION');
e_1run = e - e_2run;

e_1run.getSerie('run_nm'    ).addStim('onsets','run01.mat$','stim'   ,1)
e_2run.getSerie('run_nm_001').addStim('onsets','run01.mat$','stim_01',1)
e_2run.getSerie('run_nm_002').addStim('onsets','run02.mat$','stim_02',1)
onsetFile = e.getSerie('run_nm').getStim('stim').toJob(1);

run = e.getSerie('run_nm');
mask = run(:,1).getVolume('mask').toJob; % 1 mask per model, whatever the number of runs

TR = 1.000;


%% Models

model_list = {'mod', 'log_mod', 'dmod', 'dlog_mod'};

for l = 1 : length(model_list)
    a=e_1run.getSerie('run_nm'    ).addStim('electrophy',[     '__euclidiannorm@ACCXYZ__'  model_list{l} '.mat$'], ['acc_' model_list{l} '_01']);
    a=e_2run.getSerie('run_nm_001').addStim('electrophy',['run01__euclidiannorm@ACCXYZ__'  model_list{l} '.mat$'], ['acc_' model_list{l} '_01']);
    a=e_2run.getSerie('run_nm_001').addStim('electrophy',['run02__euclidiannorm@ACCXYZ__'  model_list{l} '.mat$'], ['acc_' model_list{l} '_02']);
end


model{1} = e.getSerie('run_nm').getStim('acc_.*_01').toJob(1);
model{2} = e.getSerie('run_nm').getStim('acc_.*_02').toJob(1);

fpath = model{1}{1}{1}; % pick one, to extract name
[~,fname,~] = fileparts(fpath);
split1 = strsplit( fname, '__' );
name = split1{2};
model_name = strcat(name,'__', model_list)';

modeldir = e.mkdir('model',model_name);
modeldir = cellfun(@cellstr, modeldir,'UniformOutput',0);

fspm = [modeldir{:}]; fspm = fspm(:);
fspm = addsuffixtofilenames(fspm, 'SPM.mat');


%% Specify user regressors

nSubj = length(e);

matlabbatch = cell(0); j = 0;
for iSubj = 1 : nSubj
    
    for m = 1 : length(model_list)
        
        j = j + 1;
        matlabbatch{j}.spm.stats.fmri_spec.dir = modeldir{iSubj}(m);
        matlabbatch{j}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{j}.spm.stats.fmri_spec.timing.RT = 1.000;
        matlabbatch{j}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{j}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        nRun = length(volume{iSubj});
        for iRun = 1 : nRun
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).scans = spm_select('expand',volume{iSubj}(iRun));
            O = load(onsetFile{iSubj}{iRun});
            M = load(model{iRun}{iSubj}{m});
            for c = 1 : length(O.names)
                matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).name     = O.names    {c};
                matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).onset    = O.onsets   {c};
                matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).duration = O.durations{c};
                if strcmp(O.names{c},'Posture')
                    posture_onset = [];
                    for block = 1 : length(O.onsets{c})
                        posture_onset = [posture_onset  O.onsets{c}(block) : TR : O.onsets{c}(block)+O.durations{c}(block)]; %#ok<AGROW>
                    end
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).onset    = posture_onset;
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).duration = ones(size(posture_onset))*TR;
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod.name  = 'TARGET';
                    bin = false(size(M.R));
                    bin(round(posture_onset)) = true;
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod.param = M.R(bin);
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod.poly  = 1;
                else
                    matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                end
                matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).tmod = 0;
                matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond(c).orth = 1;
            end
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).multi = {''};
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).multi_reg = rp{iSubj}(iRun);
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
        end % iRun
        
        matlabbatch{j}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{j}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{j}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{j}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{j}.spm.stats.fmri_spec.mthresh = 0.1;
        matlabbatch{j}.spm.stats.fmri_spec.mask = mask(iSubj);
        matlabbatch{j}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
    end % m
    
end % iSubj

return
switch CLUSTER
    case 1
        par.run = 0;
        par.sge = 1;
        par.sge_queu = 'normal,bigmem';
    case 0
        par.run = 1;
        par.sge = 0;
    case -1
        par.display = 1;
end
par.jobname  = 'spm_glm_def';

job_ending_rountines(matlabbatch,[],par);


%% Estimate

clear par
switch CLUSTER
    case 1
        par.run = 0;
        par.sge = 1;
        par.sge_queu = 'normal,bigmem';
    case 0
        par.run = 1;
        par.sge = 0;
    case -1
        par.display = 1;
end
par.display = 0;
job_first_level_estimate(fspm,par);


%% Contrast

Instructions = [1 0 0 0   0];
Relax        = [0 1 0 0   0];
Posture      = [0 0 1 0   0];
TARGET       = [0 0 0 1   0];  % <=== modulator
EndText      = [0 0 0 0   1];

contrast_T.values = {
    
Instructions
Relax
Posture
TARGET
EndText

Posture-Relax

}';

contrast_T.names = {
    
'Instructions'
'Relax'
'Posture'
'TARGET'
'EndText'

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
contrast_F.values = {eye(5)};
contrast_F.types  = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

contrast.names  = [contrast_F.names  contrast_T.names ];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types ];


%% Contrast : write

clear par
switch CLUSTER
    case 1
        par.run = 0;
        par.sge = 1;
        par.sge_queu = 'normal,bigmem';
    case 0
        par.run = 1;
        par.sge = 0;
    case -1
        par.display = 1;
end
par.display = 0;
par.jobname  = 'spm_glm_con';

par.sessrep = 'repl';

par.delete_previous = 1;
par.report=0;
job_first_level_contrast(fspm,contrast,par);

