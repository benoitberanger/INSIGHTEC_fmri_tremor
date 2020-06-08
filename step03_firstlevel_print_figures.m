clear
clc

load e.mat

model = e.getModel(cellstr2regex({'tapas_power','tapas_logpower','modulator'})).getPath';

anat = e.getSerie('anat_T1_UNI').getVolume('^wmcs').getPath;
anat = repmat(anat,size(model));

coord = {
    [ -1  -9  58]  'SMA_L'
    [-36 -25  57]   'M1_L'
    [-17 -23   6]  'VIM_L'
    [  3 -68 -12] 'CB_V_R'
    };

contrast = {
%     'Posture', 'Posture - Relax',...                       % boxcar
    'reg', 'dreg', %'Posture-Relax',...                     % EMG_envelope_power_logpower
    'Posture__power', 'Posture__power_diff',...            % modulator
    };

print_spm_figure( ...
    'model',model, ...
    'anat',anat, ...
    'coord',coord ,...
    'contrast',contrast,...
    'correction','none',...
    'threshold',0.001,...
    'nvoxel', 3,...
    'outdir',fullfile(pwd,'png'))

