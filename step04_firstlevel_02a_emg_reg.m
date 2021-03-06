clear
clc

load e.mat


%% Cluster ?

CLUSTER = 1;


%% dirs & files

volume = e.getSerie('run_nm').getVolume('^sw').toJob(1);

e.getSerie('run_nm').addJson('multiple_regressors.txt','rp')
rp = e.getSerie('run_nm').getJson('rp').toJob(1);

e.getExam( '002_FJ').getSerie('run_nm'    ).addStim('onsets',     '.mat','stim'   ,1)
e.getExam('0003_PP').getSerie('run_nm_001').addStim('onsets','run01.mat','stim_01',1)
e.getExam('0003_PP').getSerie('run_nm_002').addStim('onsets','run02.mat','stim_02',1)
onsetFile = e.getSerie('run_nm').getStim('stim').toJob(1);

run = e.getSerie('run_nm');
mask = run(:,1).getVolume('mask').toJob; % 1 mask per model, whatever the number of runs



%% Models

model_list = {'reg', 'dreg', 'log_reg', 'dlog_reg'};

for l = 1 : length(model_list)
    a=e.getExam( '002_FJ').getSerie('run_nm'    ).addStim('electrophy',[     '__peakpower@bestemg==.*__' model_list{l} '.mat$'], ['emg_' model_list{l} '_01']);
    a=e.getExam('0003_PP').getSerie('run_nm_001').addStim('electrophy',['run01__peakpower@bestemg==.*__' model_list{l} '.mat$'], ['emg_' model_list{l} '_01']);
    a=e.getExam('0003_PP').getSerie('run_nm_001').addStim('electrophy',['run02__peakpower@bestemg==.*__' model_list{l} '.mat$'], ['emg_' model_list{l} '_02']);
end

model{1} = e.getSerie('run_nm').getStim('emg_.*_01').toJob(1);
model{2} = e.getSerie('run_nm').getStim('emg_.*_02').toJob(1);

fpath = model{1}{1}{1}; % pick one, to extract name
[~,fname,~] = fileparts(fpath);
split1 = strsplit( fname, '__' );
name = split1{2};
splot2 = strsplit( name, '==' );
model_name = strcat(splot2{1},'__', model_list)';

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
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).multi = onsetFile{iSubj}(iRun);
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
            matlabbatch{j}.spm.stats.fmri_spec.sess(iRun).multi_reg = [model{iRun}{iSubj}(m) ; rp{iSubj}(iRun)];
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
EndText      = [0 0 0 1   0];
TARGET       = [0 0 0 0   1];

contrast_T.values = {
    
Instructions
Relax
Posture
EndText
TARGET


Posture-Relax

}';

contrast_T.names = {
    
'Instructions'
'Relax'
'Posture'
'EndText'
'TARGET'

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

